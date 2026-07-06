#!/usr/bin/env bash
# run_pan_aou_gwas.sh - Pan-UKB-style All of Us HapMap3 residualize-first GWAS.
#
# Runs the whole phenotype-wide GWAS from inside an AoU Verily Jupyter terminal:
#   1. resolve the pre-built HapMap3 HQ bfile and the sample-QC / PC / sex inputs;
#   2. build the unrelated-European GWAS keep-list;
#   3. extract survey responses (ds_survey [+ MHWB off-cycle for BHP/EHW]) and
#      physical measurements to local CSVs via `bq`;
#   4. hand off to pan_aou_gwas.py, which builds residualized phenotypes and runs
#      the covariate-free PLINK2 linear GWAS per phenotype.
#
# See SPECSHEET.md for the method. This script is idempotent: existing extracts
# and GWAS outputs are reused unless --force is given.
#
# Usage:
#   bash run_pan_aou_gwas.sh --setup-only   # extract + build phenotypes, no GWAS
#   bash run_pan_aou_gwas.sh --smoke        # a few phenotypes end-to-end
#   bash run_pan_aou_gwas.sh                 # full run
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Shared AoU environment (WORKSPACE_CDR, workspace bucket, DX_OUTPUT_DIR, ...).
# aou_downstream_env.sh lives one level up in the repo root.
# shellcheck source=/dev/null
source "${REPO_DIR}/aou_downstream_env.sh"

SETUP_ONLY=0
SMOKE=0
FORCE=0
SMOKE_PHENOS="${SMOKE_PHENOS:-height_cm,bmi_kg_m2,ord_educationlevel_highestgrade,ord_overallhealth_generalhealth,bin_overallhealth_generalhealth__excellent}"
for arg in "$@"; do
  case "${arg}" in
    --setup-only) SETUP_ONLY=1 ;;
    --smoke) SMOKE=1 ;;
    --force) FORCE=1 ;;
    *) echo "Unknown arg: ${arg}" >&2; exit 2 ;;
  esac
done

# --- inputs from the main pipeline ----------------------------------------- #
HM3_BFILE="${HM3_BFILE:-${DX_OUTPUT_DIR}/hapmap3_bfile_hq/hapmap3_bfile_hq}"
EUR_IIDS="${DX_OUTPUT_DIR}/europeans/classified_european_iids.txt"
FIT_PCA_IIDS="${DX_OUTPUT_DIR}/pca_eur/fit_pca_iids.txt"
PROJECTED_PCS="${DX_OUTPUT_DIR}/pca_eur/aou_projected.sscore"
SEX_COVAR="${DX_OUTPUT_DIR}/genetic_sex/sex_covar.txt"
EXCLUDE_IDENTICAL="${DX_OUTPUT_DIR}/sample_qc/exclude_identical_component_size_ge3_iids.txt"

for f in "${HM3_BFILE}.fam" "${FIT_PCA_IIDS}" "${PROJECTED_PCS}" "${SEX_COVAR}"; do
  [[ -s "${f}" ]] || { echo "ERROR: required input missing: ${f}" >&2; exit 1; }
done

# --- local working dirs ---------------------------------------------------- #
WORK="${SCRIPT_DIR}/work"
EXTRACT_DIR="${WORK}/extract"
KEEP_DIR="${WORK}/sample_qc"
mkdir -p "${EXTRACT_DIR}" "${KEEP_DIR}"

# --- 1. unrelated-EUR keep-list -------------------------------------------- #
# fit_pca_iids is the third-degree-unrelated European PCA-fit set. Intersect
# with confident genetic sex and subtract the identical-component exclusion.
KEEP="${KEEP_DIR}/unrelated_eur.keep"
if [[ ! -s "${KEEP}" || "${FORCE}" == 1 ]]; then
  echo "Building unrelated-EUR keep-list ..."
  python3 - "${FIT_PCA_IIDS}" "${SEX_COVAR}" "${EXCLUDE_IDENTICAL}" "${KEEP}" <<'PY'
import sys
fit, sexp, excl, out = sys.argv[1:5]
def col(path, idx=-1):
    s=set()
    try:
        with open(path) as f:
            for line in f:
                p=line.split()
                if p: s.add(p[idx])
    except FileNotFoundError:
        pass
    return s
fit_iids=col(fit)
# sex_covar has a header IID/sex_01; take the IID column
sex=set()
with open(sexp) as f:
    header=f.readline().split()
    iidx=header.index("IID") if "IID" in header else 1
    for line in f:
        p=line.split()
        if len(p)>iidx: sex.add(p[iidx])
excl_iids=col(excl)
keep=sorted((fit_iids & sex) - excl_iids, key=lambda x:(0,int(x)) if x.isdigit() else (1,x))
with open(out,"w") as g:
    for iid in keep:
        g.write(f"{iid}\t{iid}\n")
print(f"  keep={len(keep)}  (fit_pca={len(fit_iids)}, sex={len(sex)}, excluded={len(excl_iids)})")
PY
fi

# --- 2. extract survey responses ------------------------------------------- #
SURVEY_CSV="${EXTRACT_DIR}/survey_responses.csv"
if [[ ! -s "${SURVEY_CSV}" || "${FORCE}" == 1 ]]; then
  echo "Extracting survey responses from ${WORKSPACE_CDR}.ds_survey ..."
  bq --project_id="${GOOGLE_PROJECT}" query --nouse_legacy_sql --format=csv --max_rows=100000000 "
    SELECT
      CAST(s.person_id AS STRING)     AS person_id,
      s.survey                        AS survey,
      s.question_concept_id           AS question_concept_id,
      s.question                      AS question,
      s.answer_concept_id             AS answer_concept_id,
      s.answer                        AS answer,
      s.survey_datetime               AS survey_datetime,
      DATE_DIFF(DATE(s.survey_datetime), DATE(p.birth_datetime), DAY)/365.25 AS age_at_survey
    FROM \`${WORKSPACE_CDR}.ds_survey\` s
    JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
    WHERE s.survey_datetime IS NOT NULL
  " > "${SURVEY_CSV}"
  echo "  wrote ${SURVEY_CSV} ($(wc -l < "${SURVEY_CSV}") lines)"
fi

# Behavioral Health / Emotional Health may live in the off-cycle MHWB CDR.
# Extract it too and let the worker union it (deduped by person/question/datetime).
BHP_CSV="${EXTRACT_DIR}/bhp_ehw_responses.csv"
if [[ "${PAN_AOU_SKIP_MHWB:-0}" != 1 && ( ! -s "${BHP_CSV}" || "${FORCE}" == 1 ) ]]; then
  if bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_MHWB_CDR/./:}.survey_conduct" >/dev/null 2>&1; then
    echo "Extracting BHP/EHW responses from ${WORKSPACE_MHWB_CDR} ..."
    bq --project_id="${GOOGLE_PROJECT}" query --nouse_legacy_sql --format=csv --max_rows=100000000 "
      SELECT
        CAST(o.person_id AS STRING)                  AS person_id,
        sc.survey_concept_id                         AS survey,
        o.observation_concept_id                     AS question_concept_id,
        qc.concept_name                              AS question,
        o.value_as_concept_id                        AS answer_concept_id,
        COALESCE(a.concept_name, o.value_source_value, o.value_as_string,
                 CAST(o.value_as_number AS STRING))  AS answer,
        o.observation_datetime                       AS survey_datetime,
        DATE_DIFF(DATE(o.observation_datetime), DATE(p.birth_datetime), DAY)/365.25 AS age_at_survey
      FROM \`${WORKSPACE_MHWB_CDR}.observation\` o
      JOIN \`${WORKSPACE_MHWB_CDR}.survey_conduct\` sc
        ON o.questionnaire_response_id = sc.survey_conduct_id
      JOIN \`${WORKSPACE_MHWB_CDR}.person\` p USING (person_id)
      LEFT JOIN \`${WORKSPACE_MHWB_CDR}.concept\` qc ON o.observation_concept_id = qc.concept_id
      LEFT JOIN \`${WORKSPACE_MHWB_CDR}.concept\` a  ON o.value_as_concept_id   = a.concept_id
      WHERE o.observation_datetime IS NOT NULL
    " > "${BHP_CSV}" || echo "  WARN: MHWB extract failed; continuing without it."
  else
    echo "  MHWB CDR ${WORKSPACE_MHWB_CDR} not readable; skipping (set PAN_AOU_SKIP_MHWB=1 to silence)."
    : > "${BHP_CSV}"
  fi
fi

# --- 3. extract physical measurements -------------------------------------- #
MEAS_CSV="${EXTRACT_DIR}/measurements.csv"
if [[ ! -s "${MEAS_CSV}" || "${FORCE}" == 1 ]]; then
  echo "Extracting physical measurements from ${WORKSPACE_CDR}.measurement ..."
  bq --project_id="${GOOGLE_PROJECT}" query --nouse_legacy_sql --format=csv --max_rows=100000000 "
    SELECT
      CAST(m.person_id AS STRING) AS person_id,
      m.measurement_concept_id    AS measurement_concept_id,
      c.concept_name              AS measurement_name,
      m.measurement_datetime      AS measurement_datetime,
      m.value_as_number           AS value_as_number,
      m.unit_source_value         AS unit,
      DATE_DIFF(DATE(m.measurement_datetime), DATE(p.birth_datetime), DAY)/365.25 AS age_at_measurement
    FROM \`${WORKSPACE_CDR}.measurement\` m
    JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
    LEFT JOIN \`${WORKSPACE_CDR}.concept\` c ON m.measurement_concept_id = c.concept_id
    WHERE m.measurement_concept_id IN (
        903133,  -- Height
        903121,  -- Weight
        903124,  -- BMI
        903118,  -- Systolic BP
        903115,  -- Diastolic BP
        903126   -- Heart rate / pulse
      )
      AND m.value_as_number IS NOT NULL
  " > "${MEAS_CSV}"
  echo "  wrote ${MEAS_CSV} ($(wc -l < "${MEAS_CSV}") lines)"
fi

# --- 3b. extract Fitbit activity + sleep (optional) ------------------------ #
FITBIT_ACT_CSV="${EXTRACT_DIR}/fitbit_activity.csv"
FITBIT_SLEEP_CSV="${EXTRACT_DIR}/fitbit_sleep.csv"
if [[ "${PAN_AOU_SKIP_FITBIT:-0}" != 1 ]]; then
  if [[ ( ! -s "${FITBIT_ACT_CSV}" || "${FORCE}" == 1 ) ]] && \
     bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_CDR/./:}.activity_summary" >/dev/null 2>&1; then
    echo "Extracting Fitbit daily activity ..."
    bq --project_id="${GOOGLE_PROJECT}" query --nouse_legacy_sql --format=csv --max_rows=1000000000 "
      SELECT
        CAST(a.person_id AS STRING) AS person_id,
        a.steps                      AS steps,
        a.sedentary_minutes          AS sedentary_minutes,
        (a.fairly_active_minutes + a.very_active_minutes) AS active_minutes,
        DATE_DIFF(a.date, DATE(p.birth_datetime), DAY)/365.25 AS age
      FROM \`${WORKSPACE_CDR}.activity_summary\` a
      JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
      WHERE a.steps IS NOT NULL AND a.steps > 0
    " > "${FITBIT_ACT_CSV}" || echo "  WARN: Fitbit activity extract failed."
  fi
  if [[ ( ! -s "${FITBIT_SLEEP_CSV}" || "${FORCE}" == 1 ) ]] && \
     bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_CDR/./:}.sleep_daily_summary" >/dev/null 2>&1; then
    echo "Extracting Fitbit daily sleep ..."
    bq --project_id="${GOOGLE_PROJECT}" query --nouse_legacy_sql --format=csv --max_rows=1000000000 "
      SELECT
        CAST(s.person_id AS STRING) AS person_id,
        s.minute_asleep              AS minute_asleep,
        SAFE_DIVIDE(s.minute_asleep, NULLIF(s.minute_in_bed, 0)) AS sleep_efficiency,
        DATE_DIFF(s.sleep_date, DATE(p.birth_datetime), DAY)/365.25 AS age
      FROM \`${WORKSPACE_CDR}.sleep_daily_summary\` s
      JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
      WHERE s.is_main_sleep = true AND s.minute_asleep IS NOT NULL
    " > "${FITBIT_SLEEP_CSV}" || echo "  WARN: Fitbit sleep extract failed."
  fi
  [[ -s "${FITBIT_ACT_CSV}" ]] || echo "  (no Fitbit activity; set PAN_AOU_SKIP_FITBIT=1 to silence)"
fi

# --- 4. build phenotypes + run GWAS ---------------------------------------- #
PY_ARGS=(
  --bfile "${HM3_BFILE}"
  --keep "${KEEP}"
  --sex "${SEX_COVAR}"
  --pcs "${PROJECTED_PCS}"
  --survey-csv "${SURVEY_CSV}"
  --bhp-csv "${BHP_CSV}"
  --measurements-csv "${MEAS_CSV}"
  --fitbit-activity-csv "${FITBIT_ACT_CSV}"
  --fitbit-sleep-csv "${FITBIT_SLEEP_CSV}"
  --question-manifest "${SCRIPT_DIR}/metadata/survey_question_manifest.tsv"
  --ordinal-manifest "${SCRIPT_DIR}/metadata/ordinal_mapping_manifest.tsv"
  --pfhh-allowlist "${SCRIPT_DIR}/metadata/pfhh_self_allowlist.tsv"
  --composite-manifest "${SCRIPT_DIR}/metadata/composite_items_manifest.tsv"
  --external-scores "${SCRIPT_DIR}/metadata/external_scores.tsv"
  --outdir "${SCRIPT_DIR}"
)
# Where the ea_proxy ETM/proxy score files live (registry paths expand this).
export PAN_AOU_COG_DIR="${PAN_AOU_COG_DIR:-${DX_REGENIE_INPUT_DIR}/ses_ea_proxy/scrap/etm_cog_task_factors}"
[[ "${FORCE}" == 1 ]] && PY_ARGS+=(--force)
[[ "${SETUP_ONLY}" == 1 ]] && PY_ARGS+=(--skip-gwas)
[[ "${SMOKE}" == 1 ]] && PY_ARGS+=(--phenotypes "${SMOKE_PHENOS}")

echo "Running pan_aou_gwas.py ..."
python3 "${SCRIPT_DIR}/scripts/pan_aou_gwas.py" "${PY_ARGS[@]}"
echo "Done."
