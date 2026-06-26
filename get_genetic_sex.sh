#!/bin/bash
# get_genetic_sex.sh - Build AoU sex covariate and sex/ploidy QC summaries.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${WORKSPACE_CDR:?WORKSPACE_CDR not set}"
: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${WORKSPACE_BUCKET_URI:?WORKSPACE_BUCKET_URI not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_GENETIC_SEX_DIR:?DX_GENETIC_SEX_DIR not set}"
: "${AOU_GENOMIC_METRICS_FILE:?AOU_GENOMIC_METRICS_FILE not set}"
: "${LOCAL_REGENIE_DIR:?LOCAL_REGENIE_DIR not set}"

GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE="${GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE:-1}"
DX_GENETIC_SEX_URI="${DX_GENETIC_SEX_URI:-${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/genetic_sex}"

local_scrap="${LOCAL_REGENIE_DIR}/genetic_sex_scrap"
mkdir -p "${DX_GENETIC_SEX_DIR}/scrap" "${LOCAL_REGENIE_DIR}" "${local_scrap}"

run_wb() {
    local timeout_seconds="${WB_RESOURCE_TIMEOUT_SECONDS:-300}"
    if command -v timeout >/dev/null 2>&1; then
        timeout "${timeout_seconds}" wb "$@"
    else
        wb "$@"
    fi
}

choose_bq_tmp_dataset() {
    local requested="${1:-}" candidate
    if [[ -n "${requested}" ]]; then
        printf '%s\n' "${requested}"
        return 0
    fi
    for candidate in sbayesrc_tmp high_quality_cohort; do
        if bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${candidate}" >/dev/null 2>&1; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done
    candidate=$(bq --project_id="${GOOGLE_PROJECT}" ls --max_results=100 2>/dev/null |
        awk 'NR > 2 && $1 !~ /^-/ {print $1; exit}')
    if [[ -n "${candidate}" ]]; then
        printf '%s\n' "${candidate}"
        return 0
    fi

    candidate="sbayesrc_tmp"
    echo "  No existing BigQuery dataset found in ${GOOGLE_PROJECT}; attempting to create ${candidate} ..." >&2
    if command -v wb >/dev/null 2>&1; then
        if run_wb resource create bq-dataset \
            --id="${candidate}" \
            --dataset-id="${candidate}" \
            --location=us-central1 >/dev/null 2>&1; then
            if bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${candidate}" >/dev/null 2>&1; then
                printf '%s\n' "${candidate}"
                return 0
            fi
        fi
    fi

    if bq --project_id="${GOOGLE_PROJECT}" mk --dataset "${GOOGLE_PROJECT}:${candidate}" >/dev/null 2>&1; then
        if bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${candidate}" >/dev/null 2>&1; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    fi

    echo "ERROR: no existing BigQuery dataset found in ${GOOGLE_PROJECT}, and automatic creation of ${candidate} failed." >&2
    echo "  Create a writable BigQuery dataset in this AoU workspace or set SBAYESRC_BQ_TMP_DATASET/GENETIC_SEX_BQ_TMP_DATASET." >&2
    return 1
}

fam="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.fam"
if [[ ! -s "${fam}" ]]; then
    echo "ERROR: missing sample universe fam file ${fam}" >&2
    exit 1
fi
if [[ ! -s "${AOU_GENOMIC_METRICS_FILE}" ]]; then
    echo "ERROR: missing AoU genomic metrics file ${AOU_GENOMIC_METRICS_FILE}" >&2
    exit 1
fi

desired_params="${LOCAL_REGENIE_DIR}/genetic_sex.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'workspace_cdr\t%s\n' "${WORKSPACE_CDR}"
    printf 'sample_universe\t%s\n' "direct_bfile_hq/chr1_22_merged_hq.fam"
    printf 'sample_universe_size\t%s\n' "$(stat -c%s "${fam}")"
    printf 'genomic_metrics_file\t%s\n' "v8/wgs/short_read/snpindel/aux/qc/genomics_metrics_Dec142023_1859_02_tz0000.tsv"
    printf 'genomic_metrics_size\t%s\n' "$(stat -c%s "${AOU_GENOMIC_METRICS_FILE}")"
    printf 'require_ploidy_concordance\t%s\n' "${GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE}"
    printf 'get_genetic_sex_py_sha256\t%s\n' "$(sha256sum "${SCRIPT_DIR}/get_genetic_sex.py" | awk '{print $1}')"
} > "${desired_params}"

params="${DX_GENETIC_SEX_DIR}/genetic_sex.params.tsv"
summary="${DX_GENETIC_SEX_DIR}/genetic_sex_summary.tsv"
sex_covar="${DX_GENETIC_SEX_DIR}/sex_covar.txt"
crosstab="${DX_GENETIC_SEX_DIR}/sex_ploidy_crosstab.tsv"
qc_rows="${DX_GENETIC_SEX_DIR}/sex_ploidy_qc.tsv"
log_file="${DX_GENETIC_SEX_DIR}/genetic_sex_log.txt"

if [[ -s "${params}" && -s "${summary}" && -s "${sex_covar}" &&
      -s "${crosstab}" && -s "${qc_rows}" && -s "${log_file}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' '$1 == "confident_sex_samples" {print $2; exit}' "${summary}")
        observed=$(( $(wc -l < "${sex_covar}") - 1 ))
        if [[ -n "${expected}" && "${observed}" -eq "${expected}" ]]; then
            pct=$(awk -F'\t' '$1 == "confident_sex_percent" {print $2; exit}' "${summary}")
            echo "  Genetic sex covariate already exists (${observed} confident samples; ${pct}%) — skipping"
            exit 0
        fi
    fi
    echo "  Genetic sex outputs exist but params/counts do not match — rebuilding"
fi

sex_query="${local_scrap}/sex_at_birth_query.csv"
query_sql="${local_scrap}/sex_at_birth_query.sql"
query_log="${local_scrap}/sex_at_birth_query.log"
sample_csv="${local_scrap}/sample_universe.csv"
local_out="${local_scrap}/outputs"
tmp_dataset="$(choose_bq_tmp_dataset "${GENETIC_SEX_BQ_TMP_DATASET:-${SBAYESRC_BQ_TMP_DATASET:-}}")"
tmp_table="genetic_sex_samples_$(date +%Y%m%d_%H%M%S)_$$"
tmp_result_table="genetic_sex_result_$(date +%Y%m%d_%H%M%S)_$$"
tmp_table_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${tmp_table}"
tmp_result_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${tmp_result_table}"
tmp_table_sql="${GOOGLE_PROJECT}.${tmp_dataset}.${tmp_table}"
sex_query_gs="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/genetic_sex/scrap/sex_at_birth_query.csv"

cleanup_tmp_table() {
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_table_ref}" >/dev/null 2>&1 || true
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_result_ref}" >/dev/null 2>&1 || true
}
trap cleanup_tmp_table EXIT

{
    printf 'IID\n'
    awk '{print $2}' "${fam}"
} > "${sample_csv}"

echo "  Loading $(($(wc -l < "${sample_csv}") - 1)) sample IDs to temporary BigQuery table ${tmp_table_ref} ..."
if ! bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${tmp_dataset}" >/dev/null 2>&1; then
    echo "ERROR: BigQuery dataset ${GOOGLE_PROJECT}:${tmp_dataset} does not exist or is not readable." >&2
    echo "  Set SBAYESRC_BQ_TMP_DATASET or GENETIC_SEX_BQ_TMP_DATASET to an existing writable dataset." >&2
    exit 1
fi
bq --project_id="${GOOGLE_PROJECT}" load \
    --replace \
    --source_format=CSV \
    --skip_leading_rows=1 \
    "${tmp_table_ref}" \
    "${sample_csv}" \
    IID:STRING \
    >/dev/null

cat > "${query_sql}" <<SQL
SELECT
  CAST(p.person_id AS STRING) AS IID,
  p.sex_at_birth_concept_id,
  COALESCE(c.concept_name, '') AS sex_at_birth_concept_name,
  p.sex_at_birth_source_concept_id,
  COALESCE(sc.concept_name, '') AS sex_at_birth_source_concept_name,
  COALESCE(p.sex_at_birth_source_value, '') AS sex_at_birth_source_value
FROM \`${WORKSPACE_CDR}.person\` p
JOIN \`${tmp_table_sql}\` s
  ON s.IID = CAST(p.person_id AS STRING)
LEFT JOIN \`${WORKSPACE_CDR}.concept\` c
  ON c.concept_id = p.sex_at_birth_concept_id
LEFT JOIN \`${WORKSPACE_CDR}.concept\` sc
  ON sc.concept_id = p.sex_at_birth_source_concept_id
ORDER BY CAST(s.IID AS INT64);
SQL

echo "  Querying AoU CDR sex-at-birth fields for genotyped sample universe ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${tmp_result_ref}" \
    < "${query_sql}" \
    > "${query_log}" 2>&1

echo "  Exporting sex-at-birth query result to workspace bucket ..."
gcloud storage rm "${sex_query_gs}" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true
bq --project_id="${GOOGLE_PROJECT}" extract \
    --destination_format=CSV \
    --print_header=true \
    "${tmp_result_ref}" \
    "${sex_query_gs}" \
    >> "${query_log}" 2>&1
gcloud storage cp "${sex_query_gs}" "${sex_query}" \
    --billing-project="${GOOGLE_PROJECT}" >> "${query_log}" 2>&1

query_rows=$(( $(wc -l < "${sex_query}") - 1 ))
echo "  BigQuery sex-at-birth rows: ${query_rows}"
if [[ "${query_rows}" -le 0 ]]; then
    echo "ERROR: sex-at-birth query returned no rows" >&2
    exit 1
fi

rm -rf "${local_out}"
mkdir -p "${local_out}"

python3 "${SCRIPT_DIR}/get_genetic_sex.py" \
    --fam "${fam}" \
    --sex-query "${sex_query}" \
    --genomic-metrics "${AOU_GENOMIC_METRICS_FILE}" \
    --out-dir "${local_out}" \
    --require-ploidy-concordance "${GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE}"

cp "${desired_params}" "${local_out}/genetic_sex.params.tsv"

echo "  Copying genetic sex outputs to workspace bucket ..."
gcloud storage cp \
    "${local_out}/sex_covar.txt" \
    "${local_out}/sex_ploidy_qc.tsv" \
    "${local_out}/sex_ploidy_crosstab.tsv" \
    "${local_out}/genetic_sex_summary.tsv" \
    "${local_out}/genetic_sex_log.txt" \
    "${local_out}/genetic_sex.params.tsv" \
    "${DX_GENETIC_SEX_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null

echo "  Genetic sex summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
