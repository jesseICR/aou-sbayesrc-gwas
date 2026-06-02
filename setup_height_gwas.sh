#!/bin/bash
# setup_height_gwas.sh - Build AoU height GWAS phenotype/covariate files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${WORKSPACE_CDR:?WORKSPACE_CDR not set}"
: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${WORKSPACE_BUCKET_URI:?WORKSPACE_BUCKET_URI not set}"
: "${DX_HEIGHT_REGENIE_INPUT_DIR:?DX_HEIGHT_REGENIE_INPUT_DIR not set}"
: "${DX_HEIGHT_REGENIE_INPUT_URI:?DX_HEIGHT_REGENIE_INPUT_URI not set}"
: "${DX_EUROPEANS_DIR:?DX_EUROPEANS_DIR not set}"
: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_GENETIC_SEX_DIR:?DX_GENETIC_SEX_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${LOCAL_REGENIE_DIR:?LOCAL_REGENIE_DIR not set}"

HEIGHT_MEASUREMENT_CONCEPT_ID="${HEIGHT_MEASUREMENT_CONCEPT_ID:-3036277}"
HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID="${HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID:-903133}"
HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID="${HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID:-44818701}"
HEIGHT_UNIT_CONCEPT_ID="${HEIGHT_UNIT_CONCEPT_ID:-8582}"
HEIGHT_MIN_CM="${HEIGHT_MIN_CM:-140}"
HEIGHT_N_PCS="${HEIGHT_N_PCS:-10}"

local_scrap="${LOCAL_REGENIE_DIR}/height_example_scrap"
local_out="${local_scrap}/outputs"
mkdir -p "${DX_HEIGHT_REGENIE_INPUT_DIR}/scrap" "${LOCAL_REGENIE_DIR}" "${local_scrap}"

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
    echo "ERROR: no existing BigQuery dataset found in ${GOOGLE_PROJECT} for temporary tables." >&2
    echo "  Set SBAYESRC_BQ_TMP_DATASET or HEIGHT_BQ_TMP_DATASET to an existing writable dataset." >&2
    return 1
}

europeans="${DX_EUROPEANS_DIR}/classified_european_iids.txt"
sscore="${DX_PCA_EUR_DIR}/aou_projected.sscore"
sex_covar_input="${DX_GENETIC_SEX_DIR}/sex_covar.txt"
sex_params="${DX_GENETIC_SEX_DIR}/genetic_sex.params.tsv"
sex_summary="${DX_GENETIC_SEX_DIR}/genetic_sex_summary.tsv"
fam="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.fam"
for f in "${europeans}" "${sscore}" "${sex_covar_input}" "${sex_params}" "${sex_summary}" "${fam}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing required input ${f}" >&2
        exit 1
    fi
done

desired_params="${LOCAL_REGENIE_DIR}/height_example.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'workspace_cdr\t%s\n' "${WORKSPACE_CDR}"
    printf 'height_measurement_concept_id\t%s\n' "${HEIGHT_MEASUREMENT_CONCEPT_ID}"
    printf 'height_measurement_source_concept_id\t%s\n' "${HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID}"
    printf 'height_measurement_type_concept_id\t%s\n' "${HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID}"
    printf 'height_unit_concept_id\t%s\n' "${HEIGHT_UNIT_CONCEPT_ID}"
    printf 'height_min_cm\t%s\n' "${HEIGHT_MIN_CM}"
    printf 'height_n_pcs\t%s\n' "${HEIGHT_N_PCS}"
    printf 'europeans_size\t%s\n' "$(stat -c%s "${europeans}")"
    printf 'sscore_size\t%s\n' "$(stat -c%s "${sscore}")"
    printf 'sex_covar_size\t%s\n' "$(stat -c%s "${sex_covar_input}")"
    printf 'sex_params_size\t%s\n' "$(stat -c%s "${sex_params}")"
    printf 'sex_summary_size\t%s\n' "$(stat -c%s "${sex_summary}")"
    printf 'fam_size\t%s\n' "$(stat -c%s "${fam}")"
    printf 'setup_height_gwas_py_sha256\t%s\n' "$(sha256sum "${SCRIPT_DIR}/setup_height_gwas.py" | awk '{print $1}')"
} > "${desired_params}"

params="${DX_HEIGHT_REGENIE_INPUT_DIR}/height_gwas.params.tsv"
summary="${DX_HEIGHT_REGENIE_INPUT_DIR}/height_gwas.summary.tsv"
phen="${DX_HEIGHT_REGENIE_INPUT_DIR}/phen.txt"
base_covar="${DX_HEIGHT_REGENIE_INPUT_DIR}/base_covar.txt"
covar="${DX_HEIGHT_REGENIE_INPUT_DIR}/covar.txt"
training_iids="${DX_HEIGHT_REGENIE_INPUT_DIR}/training_iids.txt"
log_file="${DX_HEIGHT_REGENIE_INPUT_DIR}/height_gwas_log.txt"

if [[ -s "${params}" && -s "${summary}" && -s "${phen}" && -s "${base_covar}" &&
      -s "${covar}" && -s "${training_iids}" && -s "${log_file}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' '$1 == "gwas_samples" {print $2; exit}' "${summary}")
        phen_rows=$(( $(wc -l < "${phen}") - 1 ))
        covar_rows=$(( $(wc -l < "${covar}") - 1 ))
        keep_rows=$(wc -l < "${training_iids}")
        if [[ -n "${expected}" && "${phen_rows}" -eq "${expected}" &&
              "${covar_rows}" -eq "${expected}" && "${keep_rows}" -eq "${expected}" ]]; then
            echo "  Height GWAS inputs already exist (${expected} samples) — skipping"
            exit 0
        fi
    fi
    echo "  Height GWAS inputs exist but params/counts do not match — rebuilding"
fi

height_query="${local_scrap}/height_query.csv"
query_sql="${local_scrap}/height_query.sql"
query_log="${local_scrap}/height_query.log"
tmp_dataset="$(choose_bq_tmp_dataset "${HEIGHT_BQ_TMP_DATASET:-${GENETIC_SEX_BQ_TMP_DATASET:-${SBAYESRC_BQ_TMP_DATASET:-}}}")"
tmp_result_table="height_gwas_result_$(date +%Y%m%d_%H%M%S)_$$"
tmp_result_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${tmp_result_table}"
height_query_gs="${DX_HEIGHT_REGENIE_INPUT_URI}/scrap/height_query.csv"

cleanup_tmp_table() {
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_result_ref}" >/dev/null 2>&1 || true
}
trap cleanup_tmp_table EXIT

if ! bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${tmp_dataset}" >/dev/null 2>&1; then
    echo "ERROR: BigQuery dataset ${GOOGLE_PROJECT}:${tmp_dataset} does not exist or is not readable." >&2
    echo "  Set SBAYESRC_BQ_TMP_DATASET or HEIGHT_BQ_TMP_DATASET to an existing writable dataset." >&2
    exit 1
fi

cat > "${query_sql}" <<SQL
WITH height_rows AS (
  SELECT
    CAST(m.person_id AS STRING) AS IID,
    CAST(m.value_as_number AS FLOAT64) AS height_cm,
    DATE_DIFF(m.measurement_date, DATE(p.birth_datetime), DAY) / 365.25 AS age_at_height
  FROM \`${WORKSPACE_CDR}.measurement\` m
  JOIN \`${WORKSPACE_CDR}.person\` p
    ON p.person_id = m.person_id
  WHERE m.measurement_concept_id = ${HEIGHT_MEASUREMENT_CONCEPT_ID}
    AND m.measurement_source_concept_id = ${HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID}
    AND m.measurement_type_concept_id = ${HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID}
    AND m.unit_concept_id = ${HEIGHT_UNIT_CONCEPT_ID}
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number >= ${HEIGHT_MIN_CM}
    AND p.birth_datetime IS NOT NULL
)
SELECT DISTINCT
  IID,
  PERCENTILE_CONT(height_cm, 0.5) OVER (PARTITION BY IID) AS height,
  AVG(age_at_height) OVER (PARTITION BY IID) AS age_at_height,
  COUNT(*) OVER (PARTITION BY IID) AS n_height_records
FROM height_rows
ORDER BY CAST(IID AS INT64);
SQL

echo "  Querying AoU CDR for program-collected height into temporary BigQuery table ${tmp_result_ref} ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${tmp_result_ref}" \
    < "${query_sql}" \
    > "${query_log}" 2>&1

echo "  Exporting height query result to workspace bucket ..."
gcloud storage rm "${height_query_gs}" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true
bq --project_id="${GOOGLE_PROJECT}" extract \
    --destination_format=CSV \
    --print_header=true \
    "${tmp_result_ref}" \
    "${height_query_gs}" \
    >> "${query_log}" 2>&1
gcloud storage cp "${height_query_gs}" "${height_query}" \
    --billing-project="${GOOGLE_PROJECT}" >> "${query_log}" 2>&1

height_query_rows=$(( $(wc -l < "${height_query}") - 1 ))
echo "  BigQuery height rows after source/min-height filters: ${height_query_rows}"
if [[ "${height_query_rows}" -le 0 ]]; then
    echo "ERROR: height query returned no rows" >&2
    exit 1
fi

rm -rf "${local_out}"
mkdir -p "${local_out}"

python3 "${SCRIPT_DIR}/setup_height_gwas.py" \
    --height-query "${height_query}" \
    --europeans "${europeans}" \
    --sex-covar "${sex_covar_input}" \
    --fam "${fam}" \
    --sscore "${sscore}" \
    --out-dir "${local_out}" \
    --height-min-cm "${HEIGHT_MIN_CM}" \
    --n-pcs "${HEIGHT_N_PCS}"

cp "${desired_params}" "${local_out}/height_gwas.params.tsv"

echo "  Copying height GWAS inputs to workspace bucket ..."
gcloud storage cp \
    "${local_out}/training_iids.txt" \
    "${local_out}/phen.txt" \
    "${local_out}/base_covar.txt" \
    "${local_out}/covar.txt" \
    "${local_out}/height_gwas.summary.tsv" \
    "${local_out}/height_gwas_log.txt" \
    "${local_out}/height_gwas.params.tsv" \
    "${DX_HEIGHT_REGENIE_INPUT_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null

echo "  Height GWAS setup summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
