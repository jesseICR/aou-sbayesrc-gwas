#!/bin/bash
# setup_ses_ea_proxy_gwas.sh - Build primary ses_ea_proxy scores and REGENIE inputs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${WORKSPACE_CDR:?WORKSPACE_CDR not set}"
: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${WORKSPACE_BUCKET_URI:?WORKSPACE_BUCKET_URI not set}"
: "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR:?DX_SES_EA_PROXY_REGENIE_INPUT_DIR not set}"
: "${DX_SES_EA_PROXY_REGENIE_INPUT_URI:?DX_SES_EA_PROXY_REGENIE_INPUT_URI not set}"
: "${DX_EUROPEANS_DIR:?DX_EUROPEANS_DIR not set}"
: "${DX_KINSHIP_DIR:?DX_KINSHIP_DIR not set}"
: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_GENETIC_SEX_DIR:?DX_GENETIC_SEX_DIR not set}"
: "${DX_SAMPLE_QC_DIR:?DX_SAMPLE_QC_DIR not set}"
: "${DX_GWAS_STEP1_BFILE_DIR:?DX_GWAS_STEP1_BFILE_DIR not set}"
: "${LOCAL_REGENIE_DIR:?LOCAL_REGENIE_DIR not set}"

SES_EA_PROXY_N_PCS="${SES_EA_PROXY_N_PCS:-10}"
SES_EA_PROXY_OUTER_FOLDS="${SES_EA_PROXY_OUTER_FOLDS:-5}"
SES_EA_PROXY_SEED="${SES_EA_PROXY_SEED:-2026}"
SES_EA_PROXY_THREADS="${SES_EA_PROXY_THREADS:-$(nproc)}"
SES_EA_PROXY_NUM_BOOST_ROUND="${SES_EA_PROXY_NUM_BOOST_ROUND:-2000}"
SES_EA_PROXY_EARLY_STOPPING_ROUNDS="${SES_EA_PROXY_EARLY_STOPPING_ROUNDS:-50}"
SES_EA_PROXY_CV_FOLDS="${SES_EA_PROXY_CV_FOLDS:-4}"
SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD="${SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD:-0.0441941}"
GWAS_MIN_AGE_AT_SURVEY="${GWAS_MIN_AGE_AT_SURVEY:-26}"
IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE="${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE:-3}"
WORKSPACE_MHWB_CDR="${WORKSPACE_MHWB_CDR:-${WORKSPACE_CDR%%.*}.C_V8_R2_offcycle_mhwb}"

local_scrap="${LOCAL_REGENIE_DIR}/ses_ea_proxy_scrap"
local_out="${local_scrap}/outputs"
mkdir -p "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap" "${LOCAL_REGENIE_DIR}" "${local_scrap}"

choose_bq_tmp_dataset() {
    local requested="${1:-}" candidate
    if [[ -n "${requested}" ]]; then
        printf '%s\n' "${requested}"
        return 0
    fi
    for candidate in sbayesrc_tmp high_quality_cohort dataset_test2; do
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
    echo "  Set SBAYESRC_BQ_TMP_DATASET or SES_EA_PROXY_BQ_TMP_DATASET to an existing writable dataset." >&2
    return 1
}

sql_table_ref_to_bq_show_ref() {
    local ref="$1"
    printf '%s:%s\n' "${ref%%.*}" "${ref#*.}"
}

europeans="${DX_EUROPEANS_DIR}/classified_european_iids.txt"
kin0="${DX_KINSHIP_DIR}/aou_hq_direct_rel.kin0"
fit_pca_iids="${DX_PCA_EUR_DIR}/fit_pca_iids.txt"
sscore="${DX_PCA_EUR_DIR}/aou_projected.sscore"
sex_covar_input="${DX_GENETIC_SEX_DIR}/sex_covar.txt"
sex_params="${DX_GENETIC_SEX_DIR}/genetic_sex.params.tsv"
sex_summary="${DX_GENETIC_SEX_DIR}/genetic_sex_summary.tsv"
sample_qc_exclude="${DX_SAMPLE_QC_DIR}/exclude_identical_component_size_ge${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}_iids.txt"
sample_qc_params="${DX_SAMPLE_QC_DIR}/identical_component_sample_qc.params.tsv"
sample_qc_summary="${DX_SAMPLE_QC_DIR}/identical_component_sample_qc.summary.tsv"
fam="${DX_GWAS_STEP1_BFILE_DIR}/chr1_22_merged_gwas_step1.fam"
metadata="${SCRIPT_DIR}/data/aou_metadata/aou_ds_survey_question_concepts.tsv"
for f in "${europeans}" "${kin0}" "${fit_pca_iids}" "${sscore}" "${sex_covar_input}" "${sex_params}" \
         "${sex_summary}" "${sample_qc_exclude}" "${sample_qc_params}" "${sample_qc_summary}" \
         "${fam}" "${metadata}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing required input ${f}" >&2
        exit 1
    fi
done

main_cdr_show="$(sql_table_ref_to_bq_show_ref "${WORKSPACE_CDR}")"
mhwb_cdr_show="$(sql_table_ref_to_bq_show_ref "${WORKSPACE_MHWB_CDR}")"
if ! bq show "${main_cdr_show}.ds_survey" >/dev/null 2>&1; then
    echo "ERROR: cannot read ${WORKSPACE_CDR}.ds_survey" >&2
    exit 1
fi
if ! bq show "${mhwb_cdr_show}.survey_conduct" >/dev/null 2>&1; then
    echo "ERROR: cannot read ${WORKSPACE_MHWB_CDR}.survey_conduct" >&2
    echo "  Set WORKSPACE_MHWB_CDR to the Mental Health / Well-Being off-cycle dataset." >&2
    exit 1
fi

main_ids="$(PYTHONPATH="${SCRIPT_DIR}" python3 - <<'PY'
import setup_ses_ea_proxy_gwas as m
print(",".join(str(x) for x in sorted(m.MAIN_PRIMARY_IDS)))
PY
)"
bhp_codes_sql="$(PYTHONPATH="${SCRIPT_DIR}" python3 - <<'PY'
import setup_ses_ea_proxy_gwas as m
print(",".join("'" + x + "'" for x in m.BHP_CODES))
PY
)"
main_ids_hash="$(printf '%s' "${main_ids}" | sha256sum | awk '{print $1}')"
bhp_codes_hash="$(printf '%s' "${bhp_codes_sql}" | sha256sum | awk '{print $1}')"

desired_params="${LOCAL_REGENIE_DIR}/ses_ea_proxy_gwas.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'workspace_cdr\t%s\n' "${WORKSPACE_CDR}"
    printf 'workspace_mhwb_cdr\t%s\n' "${WORKSPACE_MHWB_CDR}"
    printf 'gwas_min_age_at_survey\t%s\n' "${GWAS_MIN_AGE_AT_SURVEY}"
    printf 'ses_ea_proxy_n_pcs\t%s\n' "${SES_EA_PROXY_N_PCS}"
    printf 'ses_ea_proxy_outer_folds\t%s\n' "${SES_EA_PROXY_OUTER_FOLDS}"
    printf 'ses_ea_proxy_seed\t%s\n' "${SES_EA_PROXY_SEED}"
    printf 'ses_ea_proxy_num_boost_round\t%s\n' "${SES_EA_PROXY_NUM_BOOST_ROUND}"
    printf 'ses_ea_proxy_early_stopping_rounds\t%s\n' "${SES_EA_PROXY_EARLY_STOPPING_ROUNDS}"
    printf 'ses_ea_proxy_cv_folds\t%s\n' "${SES_EA_PROXY_CV_FOLDS}"
    printf 'ses_ea_proxy_final_kinship_threshold\t%s\n' "${SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD}"
    printf 'main_feature_question_ids_sha256\t%s\n' "${main_ids_hash}"
    printf 'bhp_question_codes_sha256\t%s\n' "${bhp_codes_hash}"
    printf 'europeans_size\t%s\n' "$(stat -c%s "${europeans}")"
    printf 'kin0_size\t%s\n' "$(stat -c%s "${kin0}")"
    printf 'fit_pca_iids_size\t%s\n' "$(stat -c%s "${fit_pca_iids}")"
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
    printf 'metadata_size\t%s\n' "$(stat -c%s "${metadata}")"
    printf 'metadata_sha256\t%s\n' "$(sha256sum "${metadata}" | awk '{print $1}')"
    printf 'setup_ses_ea_proxy_gwas_sh_sha256\t%s\n' "$(sha256sum "${SCRIPT_DIR}/setup_ses_ea_proxy_gwas.sh" | awk '{print $1}')"
    printf 'setup_ses_ea_proxy_gwas_py_sha256\t%s\n' "$(sha256sum "${SCRIPT_DIR}/setup_ses_ea_proxy_gwas.py" | awk '{print $1}')"
} > "${desired_params}"

params="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/ses_ea_proxy_gwas.params.tsv"
summary="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/ses_ea_proxy_gwas.summary.tsv"
phen="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/phen.txt"
base_covar="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/base_covar.txt"
covar="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/covar.txt"
training_iids="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/training_iids.txt"
oof_scores="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/oof_scores.tsv"
applied_scores="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/applied_scores.tsv"
all_scores="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/all_scores.tsv"
fold_metrics="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/fold_metrics.tsv"
applied_metrics="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/applied_metrics.tsv"
covar_corr="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/proxy_covariate_correlations.tsv"
feature_manifest="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/feature_manifest.resolved.tsv"
feature_counts="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/feature_counts.tsv"
feature_missingness="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/feature_missingness.tsv"
pmi_missingness_counts="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/pmi_missingness_counts.tsv"
branch_recoding_summary="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/branch_recoding_summary.tsv"
missing_data_handling="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/missing_data_handling.tsv"
xgboost_model_manifest="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_model_manifest.tsv"
xgboost_feature_columns_json="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_feature_columns.json"
xgboost_feature_columns_tsv="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_feature_columns.tsv"
xgboost_final_model="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_models/final_model.json"
final_model_train_iids="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/final_model_train_iids.txt"
final_model_excluded_iids="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/final_model_excluded_related_to_applied_iids.txt"
final_model_kinholdout_summary="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/final_model_kinholdout_summary.tsv"
log_file="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/ses_ea_proxy_gwas_log.txt"

xgboost_fold_models_ok=1
for fold in $(seq 0 $((SES_EA_PROXY_OUTER_FOLDS - 1))); do
    if [[ ! -s "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_models/fold_${fold}.json" ]]; then
        xgboost_fold_models_ok=0
        break
    fi
done

if [[ -s "${params}" && -s "${summary}" && -s "${phen}" && -s "${base_covar}" &&
      -s "${covar}" && -s "${training_iids}" && -s "${oof_scores}" &&
      -s "${applied_scores}" && -s "${all_scores}" && -s "${fold_metrics}" &&
      -s "${applied_metrics}" && -s "${covar_corr}" && -s "${feature_manifest}" &&
      -s "${feature_counts}" && -s "${feature_missingness}" &&
      -s "${pmi_missingness_counts}" && -s "${branch_recoding_summary}" &&
      -s "${missing_data_handling}" && -s "${xgboost_model_manifest}" &&
      -s "${xgboost_feature_columns_json}" && -s "${xgboost_feature_columns_tsv}" &&
      -s "${xgboost_final_model}" && -s "${final_model_train_iids}" &&
      -s "${final_model_excluded_iids}" && -s "${final_model_kinholdout_summary}" &&
      "${xgboost_fold_models_ok}" -eq 1 &&
      -s "${log_file}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' '$1 == "eligible_classified_eur_samples" {print $2; exit}' "${summary}")
        phen_rows=$(( $(wc -l < "${phen}") - 1 ))
        covar_rows=$(( $(wc -l < "${covar}") - 1 ))
        keep_rows=$(wc -l < "${training_iids}")
        if [[ -n "${expected}" && "${phen_rows}" -eq "${expected}" &&
              "${covar_rows}" -eq "${expected}" && "${keep_rows}" -eq "${expected}" ]]; then
            echo "  ses_ea_proxy inputs/scores already exist (${expected} samples) — skipping"
            exit 0
        fi
    fi
    echo "  ses_ea_proxy inputs exist but params/counts do not match — rebuilding"
fi

if ! python3 - <<'PY' >/dev/null 2>&1
import xgboost
PY
then
    echo "  xgboost not found; installing with pip as approved ..."
    python3 -m pip install --user xgboost
fi

tmp_dataset="$(choose_bq_tmp_dataset "${SES_EA_PROXY_BQ_TMP_DATASET:-${SBAYESRC_BQ_TMP_DATASET:-}}")"
if ! bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${tmp_dataset}" >/dev/null 2>&1; then
    echo "ERROR: BigQuery dataset ${GOOGLE_PROJECT}:${tmp_dataset} does not exist or is not readable." >&2
    exit 1
fi

run_tag="$(date +%Y%m%d_%H%M%S)_$$"
candidate_table="ses_ea_proxy_candidates_${run_tag}"
ea_table="ses_ea_proxy_ea_${run_tag}"
main_table="ses_ea_proxy_main_${run_tag}"
bhp_table="ses_ea_proxy_bhp_${run_tag}"
area_table="ses_ea_proxy_area_${run_tag}"
candidate_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${candidate_table}"
ea_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${ea_table}"
main_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${main_table}"
bhp_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${bhp_table}"
area_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${area_table}"
candidate_sql_ref="${GOOGLE_PROJECT}.${tmp_dataset}.${candidate_table}"

cleanup_tmp_tables() {
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${candidate_ref}" >/dev/null 2>&1 || true
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${ea_ref}" >/dev/null 2>&1 || true
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${main_ref}" >/dev/null 2>&1 || true
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${bhp_ref}" >/dev/null 2>&1 || true
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${area_ref}" >/dev/null 2>&1 || true
}
trap cleanup_tmp_tables EXIT

candidate_csv="${local_scrap}/candidate_iids.csv"
{
    printf 'IID\n'
    awk '{print $NF}' "${europeans}"
} > "${candidate_csv}"

echo "  Loading classified-European candidate IID table ${candidate_ref} ..."
bq --project_id="${GOOGLE_PROJECT}" load \
    --source_format=CSV \
    --skip_leading_rows=1 \
    --replace \
    "${candidate_ref}" \
    "${candidate_csv}" \
    "IID:STRING" >/dev/null

ea_query="${local_scrap}/ea_query.csv"
main_survey="${local_scrap}/main_survey_features.csv"
bhp_survey="${local_scrap}/bhp_survey_features.csv"
area_ses="${local_scrap}/area_ses.csv"
query_log="${local_scrap}/bq_extract.log"
: > "${query_log}"

ea_query_gs="${DX_SES_EA_PROXY_REGENIE_INPUT_URI}/scrap/ea_query.csv"
main_survey_gs="${DX_SES_EA_PROXY_REGENIE_INPUT_URI}/scrap/main_survey_features.csv"
bhp_survey_gs="${DX_SES_EA_PROXY_REGENIE_INPUT_URI}/scrap/bhp_survey_features.csv"
area_ses_gs="${DX_SES_EA_PROXY_REGENIE_INPUT_URI}/scrap/area_ses.csv"

if [[ "${SES_EA_PROXY_REUSE_LOCAL_EXTRACTS:-0}" == "1" &&
      -s "${ea_query}" && -s "${main_survey}" && -s "${bhp_survey}" && -s "${area_ses}" ]]; then
    echo "  Reusing existing local survey extracts because SES_EA_PROXY_REUSE_LOCAL_EXTRACTS=1"
else
echo "  Querying EA teacher label ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${ea_ref}" \
    >> "${query_log}" 2>&1 <<SQL
WITH codeable AS (
  SELECT
    CAST(s.person_id AS STRING) AS IID,
    CAST(CASE s.answer_concept_id
      WHEN 1585941 THEN 9.0
      WHEN 1585942 THEN 9.0
      WHEN 1585943 THEN 9.0
      WHEN 1585944 THEN 10.0
      WHEN 1585945 THEN 13.0
      WHEN 1585946 THEN 15.0
      WHEN 1585947 THEN 18.0
      WHEN 1585948 THEN 20.0
      ELSE NULL
    END AS FLOAT64) AS ea_years,
    EXTRACT(YEAR FROM DATE(p.birth_datetime)) +
      SAFE_DIVIDE(
        DATE_DIFF(DATE(p.birth_datetime), DATE(EXTRACT(YEAR FROM DATE(p.birth_datetime)), 1, 1), DAY),
        DATE_DIFF(DATE(EXTRACT(YEAR FROM DATE(p.birth_datetime)) + 1, 1, 1), DATE(EXTRACT(YEAR FROM DATE(p.birth_datetime)), 1, 1), DAY)
      ) AS yob,
    DATE_DIFF(DATE(s.survey_datetime), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_basics,
    CAST(s.answer_concept_id AS INT64) AS answer_concept_id,
    s.answer AS answer,
    s.survey_datetime,
    ROW_NUMBER() OVER (
      PARTITION BY s.person_id
      ORDER BY s.survey_datetime DESC, s.answer_concept_id DESC
    ) AS rn
  FROM \`${WORKSPACE_CDR}.ds_survey\` s
  JOIN \`${WORKSPACE_CDR}.person\` p
    ON p.person_id = s.person_id
  JOIN \`${candidate_sql_ref}\` c
    ON c.IID = CAST(s.person_id AS STRING)
  WHERE s.question_concept_id = 1585940
    AND s.answer_concept_id IN (1585941,1585942,1585943,1585944,1585945,1585946,1585947,1585948)
    AND s.survey_datetime IS NOT NULL
    AND p.birth_datetime IS NOT NULL
)
SELECT IID, ea_years, yob, age_at_basics, answer_concept_id, answer
FROM codeable
WHERE rn = 1
ORDER BY CAST(IID AS INT64)
SQL

echo "  Querying primary main-CDR survey features ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${main_ref}" \
    >> "${query_log}" 2>&1 <<SQL
WITH latest_q AS (
  SELECT
    s.person_id,
    s.question_concept_id,
    MAX(s.survey_datetime) AS latest_dt
  FROM \`${WORKSPACE_CDR}.ds_survey\` s
  JOIN \`${candidate_sql_ref}\` c
    ON c.IID = CAST(s.person_id AS STRING)
  WHERE s.question_concept_id IN (${main_ids})
    AND s.survey_datetime IS NOT NULL
  GROUP BY s.person_id, s.question_concept_id
)
SELECT
  CAST(s.person_id AS STRING) AS IID,
  s.survey,
  CAST(s.question_concept_id AS INT64) AS question_concept_id,
  s.question,
  CAST(s.answer_concept_id AS INT64) AS answer_concept_id,
  s.answer,
  s.survey_datetime,
  DATE_DIFF(DATE(s.survey_datetime), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_survey
FROM \`${WORKSPACE_CDR}.ds_survey\` s
JOIN latest_q l
  ON s.person_id = l.person_id
 AND s.question_concept_id = l.question_concept_id
 AND s.survey_datetime = l.latest_dt
JOIN \`${WORKSPACE_CDR}.person\` p
  ON p.person_id = s.person_id
ORDER BY CAST(s.person_id AS INT64), s.survey, s.question_concept_id, s.answer_concept_id
SQL

echo "  Querying BHP primary features from ${WORKSPACE_MHWB_CDR} ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${bhp_ref}" \
    >> "${query_log}" 2>&1 <<SQL
WITH bhp_rows AS (
  SELECT
    o.person_id,
    o.observation_concept_id AS question_concept_id,
    q.concept_code AS question_code,
    q.concept_name AS question,
    o.value_as_concept_id AS answer_concept_id,
    COALESCE(a.concept_name, o.value_source_value, o.value_as_string, CAST(o.value_as_number AS STRING)) AS answer,
    o.observation_datetime AS survey_datetime,
    DATE_DIFF(DATE(o.observation_datetime), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_survey
  FROM \`${WORKSPACE_MHWB_CDR}.observation\` o
  JOIN \`${WORKSPACE_MHWB_CDR}.survey_conduct\` sc
    ON o.questionnaire_response_id = sc.survey_conduct_id
  JOIN \`${WORKSPACE_MHWB_CDR}.concept\` q
    ON q.concept_id = o.observation_concept_id
  LEFT JOIN \`${WORKSPACE_MHWB_CDR}.concept\` a
    ON a.concept_id = o.value_as_concept_id
  JOIN \`${WORKSPACE_CDR}.person\` p
    ON p.person_id = o.person_id
  JOIN \`${candidate_sql_ref}\` c
    ON c.IID = CAST(o.person_id AS STRING)
  WHERE sc.survey_concept_id = 1703870
    AND sc.survey_source_value = 'bhp'
    AND q.concept_code IN (${bhp_codes_sql})
    AND o.observation_datetime IS NOT NULL
    AND p.birth_datetime IS NOT NULL
),
latest_q AS (
  SELECT person_id, question_concept_id, MAX(survey_datetime) AS latest_dt
  FROM bhp_rows
  GROUP BY person_id, question_concept_id
)
SELECT
  CAST(b.person_id AS STRING) AS IID,
  CAST(b.question_concept_id AS INT64) AS question_concept_id,
  b.question_code,
  b.question,
  CAST(b.answer_concept_id AS INT64) AS answer_concept_id,
  b.answer,
  b.survey_datetime,
  b.age_at_survey
FROM bhp_rows b
JOIN latest_q l
  ON b.person_id = l.person_id
 AND b.question_concept_id = l.question_concept_id
 AND b.survey_datetime = l.latest_dt
ORDER BY CAST(b.person_id AS INT64), b.question_code, b.answer_concept_id
SQL

echo "  Querying ZIP3-derived area SES ..."
bq --project_id="${GOOGLE_PROJECT}" query \
    --use_legacy_sql=false \
    --replace \
    --destination_table="${area_ref}" \
    >> "${query_log}" 2>&1 <<SQL
WITH ranked AS (
  SELECT
    CAST(z.PERSON_ID AS STRING) AS IID,
    CAST(z.DEPRIVATION_INDEX AS FLOAT64) AS deprivation_index,
    CAST(z.MEDIAN_INCOME AS FLOAT64) AS median_income,
    CAST(z.FRACTION_POVERTY AS FLOAT64) AS fraction_poverty,
    CAST(z.FRACTION_ASSISTED_INCOME AS FLOAT64) AS fraction_assisted_income,
    CAST(z.FRACTION_NO_HEALTH_INS AS FLOAT64) AS fraction_no_health_ins,
    CAST(z.FRACTION_VACANT_HOUSING AS FLOAT64) AS fraction_vacant_housing,
    ROW_NUMBER() OVER (
      PARTITION BY z.PERSON_ID
      ORDER BY z.OBSERVATION_DATETIME DESC, z.ACS DESC
    ) AS rn
  FROM \`${WORKSPACE_CDR}.ds_zip_code_socioeconomic\` z
  JOIN \`${candidate_sql_ref}\` c
    ON c.IID = CAST(z.PERSON_ID AS STRING)
)
SELECT
  IID,
  deprivation_index,
  median_income,
  fraction_poverty,
  fraction_assisted_income,
  fraction_no_health_ins,
  fraction_vacant_housing
FROM ranked
WHERE rn = 1
ORDER BY CAST(IID AS INT64)
SQL

copy_extract() {
    local table_ref="$1" gs_uri="$2" local_path="$3" shard_uri part_dir
    shard_uri="${gs_uri%.csv}-*.csv"
    part_dir="${local_path}.parts"
    rm -rf "${part_dir}"
    mkdir -p "${part_dir}"
    gcloud storage rm "${shard_uri}" --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true
    bq --project_id="${GOOGLE_PROJECT}" extract \
        --destination_format=CSV \
        --print_header=true \
        "${table_ref}" \
        "${shard_uri}" >> "${query_log}" 2>&1
    gcloud storage cp "${shard_uri}" "${part_dir}/" \
        --billing-project="${GOOGLE_PROJECT}" >> "${query_log}" 2>&1
    mapfile -t parts < <(find "${part_dir}" -maxdepth 1 -type f -name '*.csv' | sort)
    if [[ "${#parts[@]}" -eq 0 ]]; then
        echo "ERROR: no export shards found for ${table_ref}" >&2
        return 1
    fi
    awk 'FNR == 1 && NR != 1 {next} {print}' "${parts[@]}" > "${local_path}"
    rm -rf "${part_dir}"
}

echo "  Exporting query results to workspace-local files ..."
copy_extract "${ea_ref}" "${ea_query_gs}" "${ea_query}"
copy_extract "${main_ref}" "${main_survey_gs}" "${main_survey}"
copy_extract "${bhp_ref}" "${bhp_survey_gs}" "${bhp_survey}"
copy_extract "${area_ref}" "${area_ses_gs}" "${area_ses}"
fi

echo "  Extract row counts:"
printf '    EA rows:          %s\n' "$(( $(wc -l < "${ea_query}") - 1 ))"
printf '    Main survey rows: %s\n' "$(( $(wc -l < "${main_survey}") - 1 ))"
printf '    BHP rows:         %s\n' "$(( $(wc -l < "${bhp_survey}") - 1 ))"
printf '    Area SES rows:    %s\n' "$(( $(wc -l < "${area_ses}") - 1 ))"

rm -rf "${local_out}"
mkdir -p "${local_out}"

python3 "${SCRIPT_DIR}/setup_ses_ea_proxy_gwas.py" \
    --ea-query "${ea_query}" \
    --main-survey "${main_survey}" \
    --bhp-survey "${bhp_survey}" \
    --area-ses "${area_ses}" \
    --metadata "${metadata}" \
    --europeans "${europeans}" \
    --final-kinship-holdout-kin0 "${kin0}" \
    --final-kinship-holdout-threshold "${SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD}" \
    --fit-pca-iids "${fit_pca_iids}" \
    --sex-covar "${sex_covar_input}" \
    --exclude-iids "${sample_qc_exclude}" \
    --fam "${fam}" \
    --sscore "${sscore}" \
    --out-dir "${local_out}" \
    --n-pcs "${SES_EA_PROXY_N_PCS}" \
    --min-age-at-basics "${GWAS_MIN_AGE_AT_SURVEY}" \
    --outer-folds "${SES_EA_PROXY_OUTER_FOLDS}" \
    --seed "${SES_EA_PROXY_SEED}" \
    --threads "${SES_EA_PROXY_THREADS}" \
    --num-boost-round "${SES_EA_PROXY_NUM_BOOST_ROUND}" \
    --early-stopping-rounds "${SES_EA_PROXY_EARLY_STOPPING_ROUNDS}" \
    --cv-folds "${SES_EA_PROXY_CV_FOLDS}"

cp "${desired_params}" "${local_out}/ses_ea_proxy_gwas.params.tsv"

echo "  Copying ses_ea_proxy outputs to workspace bucket ..."
gcloud storage cp \
    "${local_out}/training_iids.txt" \
    "${local_out}/phen.txt" \
    "${local_out}/base_covar.txt" \
    "${local_out}/covar.txt" \
    "${local_out}/oof_scores.tsv" \
    "${local_out}/applied_scores.tsv" \
    "${local_out}/all_scores.tsv" \
    "${local_out}/fold_assignment.tsv" \
    "${local_out}/fold_metrics.tsv" \
    "${local_out}/applied_metrics.tsv" \
    "${local_out}/proxy_covariate_correlations.tsv" \
    "${local_out}/feature_importance.tsv" \
    "${local_out}/final_model_train_iids.txt" \
    "${local_out}/final_model_excluded_related_to_applied_iids.txt" \
    "${local_out}/final_model_excluded_related_to_applied_edges.tsv" \
    "${local_out}/final_model_kinholdout_summary.tsv" \
    "${local_out}/feature_manifest.resolved.tsv" \
    "${local_out}/feature_counts.tsv" \
    "${local_out}/feature_missingness.tsv" \
    "${local_out}/pmi_missingness_counts.tsv" \
    "${local_out}/branch_recoding_summary.tsv" \
    "${local_out}/missing_data_handling.tsv" \
    "${local_out}/xgboost_model_manifest.tsv" \
    "${local_out}/xgboost_feature_columns.json" \
    "${local_out}/xgboost_feature_columns.tsv" \
    "${local_out}/leakage_denylist_hits.tsv" \
    "${local_out}/runtime_manifest.json" \
    "${local_out}/ses_ea_proxy_gwas.summary.tsv" \
    "${local_out}/ses_ea_proxy_gwas_log.txt" \
    "${local_out}/ses_ea_proxy_gwas.params.tsv" \
    "${DX_SES_EA_PROXY_REGENIE_INPUT_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null

gcloud storage cp -r \
    "${local_out}/xgboost_models" \
    "${DX_SES_EA_PROXY_REGENIE_INPUT_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null

echo "  ses_ea_proxy setup summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
