#!/bin/bash
# setup_income_gwas.sh - Build AoU household-income GWAS inputs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${WORKSPACE_CDR:?WORKSPACE_CDR not set}"
: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${WORKSPACE_BUCKET_URI:?WORKSPACE_BUCKET_URI not set}"
: "${DX_INCOME_REGENIE_INPUT_DIR:?DX_INCOME_REGENIE_INPUT_DIR not set}"
: "${DX_INCOME_REGENIE_INPUT_URI:?DX_INCOME_REGENIE_INPUT_URI not set}"
: "${DX_EUROPEANS_DIR:?DX_EUROPEANS_DIR not set}"
: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_GENETIC_SEX_DIR:?DX_GENETIC_SEX_DIR not set}"
: "${DX_SAMPLE_QC_DIR:?DX_SAMPLE_QC_DIR not set}"
: "${DX_GWAS_STEP1_BFILE_DIR:?DX_GWAS_STEP1_BFILE_DIR not set}"
: "${LOCAL_REGENIE_DIR:?LOCAL_REGENIE_DIR not set}"

INCOME_QUESTION_CONCEPT_ID="${INCOME_QUESTION_CONCEPT_ID:-1585375}"
INCOME_N_PCS="${INCOME_N_PCS:-10}"
IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE="${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE:-3}"

local_scrap="${LOCAL_REGENIE_DIR}/income_gwas_scrap"
local_out="${local_scrap}/outputs"
mkdir -p "${DX_INCOME_REGENIE_INPUT_DIR}/scrap" "${LOCAL_REGENIE_DIR}" "${local_scrap}"

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
    echo "  Set SBAYESRC_BQ_TMP_DATASET or INCOME_BQ_TMP_DATASET to an existing writable dataset." >&2
    return 1
}

europeans="${DX_EUROPEANS_DIR}/classified_european_iids.txt"
sscore="${DX_PCA_EUR_DIR}/aou_projected.sscore"
sex_covar_input="${DX_GENETIC_SEX_DIR}/sex_covar.txt"
sex_params="${DX_GENETIC_SEX_DIR}/genetic_sex.params.tsv"
sex_summary="${DX_GENETIC_SEX_DIR}/genetic_sex_summary.tsv"
sample_qc_exclude="${DX_SAMPLE_QC_DIR}/exclude_identical_component_size_ge${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}_iids.txt"
sample_qc_params="${DX_SAMPLE_QC_DIR}/identical_component_sample_qc.params.tsv"
sample_qc_summary="${DX_SAMPLE_QC_DIR}/identical_component_sample_qc.summary.tsv"
fam="${DX_GWAS_STEP1_BFILE_DIR}/chr1_22_merged_gwas_step1.fam"
for f in "${europeans}" "${sscore}" "${sex_covar_input}" "${sex_params}" "${sex_summary}" \
         "${sample_qc_exclude}" "${sample_qc_params}" "${sample_qc_summary}" "${fam}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing required input ${f}" >&2
        exit 1
    fi
done

desired_params="${LOCAL_REGENIE_DIR}/income_gwas.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'workspace_cdr\t%s\n' "${WORKSPACE_CDR}"
    printf 'income_question_concept_id\t%s\n' "${INCOME_QUESTION_CONCEPT_ID}"
    printf 'income_mapping\t%s\n' "1585376=5;1585377=17.5;1585378=30;1585379=42.5;1585380=62.5;1585381=87.5;1585382=125;1585383=175;1585384=250"
    printf 'income_n_pcs\t%s\n' "${INCOME_N_PCS}"
    printf 'europeans_size\t%s\n' "$(stat -c%s "${europeans}")"
    printf 'sscore_size\t%s\n' "$(stat -c%s "${sscore}")"
    printf 'sex_covar_size\t%s\n' "$(stat -c%s "${sex_covar_input}")"
    printf 'sex_params_size\t%s\n' "$(stat -c%s "${sex_params}")"
    printf 'sex_summary_size\t%s\n' "$(stat -c%s "${sex_summary}")"
    printf 'sample_qc_exclude_size\t%s\n' "$(stat -c%s "${sample_qc_exclude}")"
    printf 'sample_qc_exclude_sha256\t%s\n' "$(sha256sum "${sample_qc_exclude}" | awk '{print $1}')"
    printf 'sample_qc_params_size\t%s\n' "$(stat -c%s "${sample_qc_params}")"
    printf 'sample_qc_params_sha256\t%s\n' "$(sha256sum "${sample_qc_params}" | awk '{print $1}')"
    printf 'sample_qc_summary_size\t%s\n' "$(stat -c%s "${sample_qc_summary}")"
    printf 'sample_qc_summary_sha256\t%s\n' "$(sha256sum "${sample_qc_summary}" | awk '{print $1}')"
    printf 'identical_component_exclude_min_size\t%s\n' "${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}"
    printf 'fam_size\t%s\n' "$(stat -c%s "${fam}")"
    printf 'setup_income_gwas_py_sha256\t%s\n' "$(sha256sum "${SCRIPT_DIR}/setup_income_gwas.py" | awk '{print $1}')"
} > "${desired_params}"

params="${DX_INCOME_REGENIE_INPUT_DIR}/income_gwas.params.tsv"
summary="${DX_INCOME_REGENIE_INPUT_DIR}/income_gwas.summary.tsv"
phen="${DX_INCOME_REGENIE_INPUT_DIR}/phen.txt"
base_covar="${DX_INCOME_REGENIE_INPUT_DIR}/base_covar.txt"
covar="${DX_INCOME_REGENIE_INPUT_DIR}/covar.txt"
training_iids="${DX_INCOME_REGENIE_INPUT_DIR}/training_iids.txt"
answer_counts="${DX_INCOME_REGENIE_INPUT_DIR}/income_answer_counts.tsv"
log_file="${DX_INCOME_REGENIE_INPUT_DIR}/income_gwas_log.txt"

if [[ -s "${params}" && -s "${summary}" && -s "${phen}" && -s "${base_covar}" &&
      -s "${covar}" && -s "${training_iids}" && -s "${answer_counts}" && -s "${log_file}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' '$1 == "gwas_samples" {print $2; exit}' "${summary}")
        phen_rows=$(( $(wc -l < "${phen}") - 1 ))
        covar_rows=$(( $(wc -l < "${covar}") - 1 ))
        keep_rows=$(wc -l < "${training_iids}")
        if [[ -n "${expected}" && "${phen_rows}" -eq "${expected}" &&
              "${covar_rows}" -eq "${expected}" && "${keep_rows}" -eq "${expected}" ]]; then
            echo "  Income GWAS inputs already exist (${expected} samples) — skipping"
            exit 0
        fi
    fi
    echo "  Income GWAS inputs exist but params/counts do not match — rebuilding"
fi

income_query="${local_scrap}/income_query.csv"
query_sql="${local_scrap}/income_query.sql"
query_log="${local_scrap}/income_query.log"
tmp_dataset="$(choose_bq_tmp_dataset "${INCOME_BQ_TMP_DATASET:-${SBAYESRC_BQ_TMP_DATASET:-}}")"
tmp_result_table="income_gwas_result_$(date +%Y%m%d_%H%M%S)_$$"
tmp_result_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${tmp_result_table}"
income_query_gs="${DX_INCOME_REGENIE_INPUT_URI}/scrap/income_query.csv"

cleanup_tmp_table() {
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_result_ref}" >/dev/null 2>&1 || true
}
trap cleanup_tmp_table EXIT

if ! bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${tmp_dataset}" >/dev/null 2>&1; then
    echo "ERROR: BigQuery dataset ${GOOGLE_PROJECT}:${tmp_dataset} does not exist or is not readable." >&2
    echo "  Set SBAYESRC_BQ_TMP_DATASET or INCOME_BQ_TMP_DATASET to an existing writable dataset." >&2
    exit 1
fi

cat > "${query_sql}" <<SQL
WITH codeable AS (
  SELECT
    CAST(s.person_id AS STRING) AS IID,
    CAST(CASE s.answer_concept_id
      WHEN 1585376 THEN 5.0
      WHEN 1585377 THEN 17.5
      WHEN 1585378 THEN 30.0
      WHEN 1585379 THEN 42.5
      WHEN 1585380 THEN 62.5
      WHEN 1585381 THEN 87.5
      WHEN 1585382 THEN 125.0
      WHEN 1585383 THEN 175.0
      WHEN 1585384 THEN 250.0
      ELSE NULL
    END AS FLOAT64) AS income_k,
    DATE_DIFF(DATE(s.survey_datetime), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_survey,
    CAST(s.answer_concept_id AS INT64) AS answer_concept_id,
    s.answer AS answer,
    s.survey_datetime,
    COUNT(*) OVER (PARTITION BY s.person_id) AS n_income_records,
    ROW_NUMBER() OVER (
      PARTITION BY s.person_id
      ORDER BY s.survey_datetime DESC, s.answer_concept_id DESC
    ) AS rn
  FROM \`${WORKSPACE_CDR}.ds_survey\` s
  JOIN \`${WORKSPACE_CDR}.person\` p
    ON p.person_id = s.person_id
  WHERE s.question_concept_id = ${INCOME_QUESTION_CONCEPT_ID}
    AND s.answer_concept_id IN (1585376,1585377,1585378,1585379,1585380,1585381,1585382,1585383,1585384)
    AND s.survey_datetime IS NOT NULL
    AND p.birth_datetime IS NOT NULL
)
SELECT
  IID,
  income_k,
  age_at_survey,
  answer_concept_id,
  answer,
  n_income_records
FROM codeable
WHERE rn = 1
ORDER BY CAST(IID AS INT64);
SQL

echo "  Querying AoU CDR for income phenotype into temporary BigQuery table ${tmp_result_ref} ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${tmp_result_ref}" \
    < "${query_sql}" \
    > "${query_log}" 2>&1

echo "  Exporting income query result to workspace bucket ..."
gcloud storage rm "${income_query_gs}" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true
bq --project_id="${GOOGLE_PROJECT}" extract \
    --destination_format=CSV \
    --print_header=true \
    "${tmp_result_ref}" \
    "${income_query_gs}" \
    >> "${query_log}" 2>&1
gcloud storage cp "${income_query_gs}" "${income_query}" \
    --billing-project="${GOOGLE_PROJECT}" >> "${query_log}" 2>&1

income_query_rows=$(( $(wc -l < "${income_query}") - 1 ))
echo "  BigQuery income rows after codeable-answer/birth-date filters: ${income_query_rows}"
if [[ "${income_query_rows}" -le 0 ]]; then
    echo "ERROR: income query returned no rows" >&2
    exit 1
fi

rm -rf "${local_out}"
mkdir -p "${local_out}"

python3 "${SCRIPT_DIR}/setup_income_gwas.py" \
    --income-query "${income_query}" \
    --europeans "${europeans}" \
    --sex-covar "${sex_covar_input}" \
    --exclude-iids "${sample_qc_exclude}" \
    --fam "${fam}" \
    --sscore "${sscore}" \
    --out-dir "${local_out}" \
    --n-pcs "${INCOME_N_PCS}"

cp "${desired_params}" "${local_out}/income_gwas.params.tsv"

echo "  Copying income GWAS inputs to workspace bucket ..."
gcloud storage cp \
    "${local_out}/training_iids.txt" \
    "${local_out}/phen.txt" \
    "${local_out}/base_covar.txt" \
    "${local_out}/covar.txt" \
    "${local_out}/income_answer_counts.tsv" \
    "${local_out}/income_gwas.summary.tsv" \
    "${local_out}/income_gwas_log.txt" \
    "${local_out}/income_gwas.params.tsv" \
    "${DX_INCOME_REGENIE_INPUT_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null

echo "  Income GWAS setup summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
