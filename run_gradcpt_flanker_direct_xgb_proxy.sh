#!/bin/bash
# run_gradcpt_flanker_direct_xgb_proxy.sh - Scratch survey XGBoost for GradCPT/Flanker.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR
source "${SCRIPT_DIR}/aou_downstream_env.sh"

usage() {
    cat <<'EOF'
Usage: bash run_gradcpt_flanker_direct_xgb_proxy.sh [OPTIONS]

Trains a scratch XGBoost model from survey features to predict the
missing-pattern-aware GradCPT/Flanker target. The final no-teacher calibrated
phenotype is built by run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas.sh.

Options:
  --force                    Recompute and overwrite this output dir.
  --stage-aggregate          Stage aggregate diagnostics to ses_ea_proxy/scrap.
  --no-write-regenie-inputs  Do not write phen/covar/training_iids files.
  -h, --help                 Show this help.

Environment overrides:
  SES_EA_PROXY_GWAS_INPUT_NAME          Default: ses_ea_proxy_v2_kinholdout
  DIRECT_XGB_OUTPUT_NAME                Default: gradcpt_flanker_direct_xgb_proxy_${SES_EA_PROXY_GWAS_INPUT_NAME}
  DIRECT_XGB_ETA                       Default: 0.05
  DIRECT_XGB_MAX_DEPTH                 Default: 6
  DIRECT_XGB_MIN_CHILD_WEIGHT          Default: 20
  DIRECT_XGB_LAMBDA                    Default: 1
  DIRECT_XGB_NUM_BOOST_ROUND           Default: 2000
  DIRECT_XGB_EARLY_STOPPING_ROUNDS     Default: 50
  DIRECT_XGB_CV_FOLDS                  Default: 4
EOF
}

FORCE=0
STAGE_AGGREGATE=0
NO_WRITE_REGENIE_INPUTS=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)
            FORCE=1
            shift
            ;;
        --stage-aggregate)
            STAGE_AGGREGATE=1
            shift
            ;;
        --no-write-regenie-inputs)
            NO_WRITE_REGENIE_INPUTS=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

WORKSPACE_BUCKET_MOUNT="/home/jupyter/workspace/workspace-bucket"
WORKSPACE_BUCKET_URI="gs://$(mount | awk '/ on \/home\/jupyter\/workspace\/workspace-bucket /{print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    exit 1
fi

DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/sbayesrc_genotypes"
DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
SES_EA_PROXY_GWAS_INPUT_NAME="${SES_EA_PROXY_GWAS_INPUT_NAME:-ses_ea_proxy_v2_kinholdout}"
DX_SES_EA_PROXY_REGENIE_INPUT_DIR="${DX_REGENIE_INPUT_DIR}/${SES_EA_PROXY_GWAS_INPUT_NAME}"
LOCAL_REGENIE_DIR="${SCRIPT_DIR}/data/regenie"
LOCAL_SCRAP="${LOCAL_REGENIE_DIR}/ses_ea_proxy_scrap"
LOCAL_ETM_TASK="${LOCAL_SCRAP}/etm_cog_task_factors/etm_cog_task_factors_recommended_wide.tsv"
FINE_TUNED_INPUT_NAME="gradcpt_flanker_finetuned_ea_proxy_${SES_EA_PROXY_GWAS_INPUT_NAME}"
FINE_TUNED_SCORES="${DX_REGENIE_INPUT_DIR}/${FINE_TUNED_INPUT_NAME}/g4_finetuned_ea_proxy_scores_wide.tsv"
METADATA="${SCRIPT_DIR}/data/aou_metadata/aou_ds_survey_question_concepts.tsv"
WORKSPACE_EXTRACT_SCRAP="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap"

DIRECT_XGB_OUTPUT_NAME="${DIRECT_XGB_OUTPUT_NAME:-gradcpt_flanker_direct_xgb_proxy_${SES_EA_PROXY_GWAS_INPUT_NAME}}"
LOCAL_OUT="${LOCAL_SCRAP}/${DIRECT_XGB_OUTPUT_NAME}"
WORKSPACE_OUT="${DX_REGENIE_INPUT_DIR}/${DIRECT_XGB_OUTPUT_NAME}"
WORKSPACE_DIAG="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap/${DIRECT_XGB_OUTPUT_NAME}"

ALL_SCORES="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/all_scores.tsv"

required=(
    "${ALL_SCORES}"
    "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/base_covar.txt"
    "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/covar.txt"
    "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_feature_columns.json"
    "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_model_manifest.tsv"
    "${LOCAL_ETM_TASK}"
    "${FINE_TUNED_SCORES}"
    "${METADATA}"
)
for f in "${required[@]}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing required input ${f}" >&2
        echo "  This command does not query BigQuery. Rerun the SES-EA proxy and GradCPT/Flanker fine-tune setup first." >&2
        exit 1
    fi
done

mkdir -p "${SCRIPT_DIR}/logs" "${LOCAL_OUT}"
if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
    mkdir -p "${WORKSPACE_DIAG}"
fi

LOG_FILE="${SCRIPT_DIR}/logs/run_gradcpt_flanker_direct_xgb_proxy_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

ensure_extract() {
    local name="$1"
    local local_path="${LOCAL_SCRAP}/${name}"
    local workspace_path="${WORKSPACE_EXTRACT_SCRAP}/${name}"
    local stem="${name%.csv}"
    if [[ -s "${local_path}" ]]; then
        return 0
    fi
    if [[ -s "${workspace_path}" ]]; then
        echo "  Localizing ${name} from workspace exact file ..."
        cp "${workspace_path}" "${local_path}"
        return 0
    fi
    shopt -s nullglob
    local parts=("${WORKSPACE_EXTRACT_SCRAP}/${stem}-"*.csv)
    shopt -u nullglob
    if [[ "${#parts[@]}" -gt 0 ]]; then
        echo "  Reconstructing ${name} from ${#parts[@]} workspace shard(s) ..."
        awk 'FNR == 1 && NR != 1 {next} {print}' "${parts[@]}" > "${local_path}"
        return 0
    fi
    echo "ERROR: missing ${name} locally and in ${WORKSPACE_EXTRACT_SCRAP}" >&2
    echo "  Rerun bash run_ses_ea_proxy_gwas.sh --setup-only to regenerate survey extracts." >&2
    exit 1
}

ensure_extract "ea_query.csv"
ensure_extract "main_survey_features.csv"
ensure_extract "bhp_survey_features.csv"
ensure_extract "area_ses.csv"

if ! python3 - <<'PY' >/dev/null 2>&1
import xgboost
PY
then
    echo "  xgboost not found; installing with pip as approved by the SES-EA setup command ..."
    python3 -m pip install --user xgboost
fi

echo "================================================"
echo "AoU direct GradCPT/Flanker XGBoost proxy — $(date)"
echo "================================================"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
echo "  SES-EA input dir     = ${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}"
echo "  ETM task score file  = ${LOCAL_ETM_TASK}"
echo "  Fine-tuned score file= ${FINE_TUNED_SCORES}"
echo "  Local output dir     = ${LOCAL_OUT}"
echo "  Workspace output dir = ${WORKSPACE_OUT}"
echo "  Workspace diag dir   = ${WORKSPACE_DIAG}"
echo "  FORCE                = ${FORCE}"
echo "  STAGE_AGGREGATE      = ${STAGE_AGGREGATE}"
echo "  WRITE_REGENIE_INPUTS = $((1 - NO_WRITE_REGENIE_INPUTS))"
echo "  ETA                  = ${DIRECT_XGB_ETA:-0.05}"
echo "  MAX_DEPTH            = ${DIRECT_XGB_MAX_DEPTH:-6}"
echo "  MIN_CHILD_WEIGHT     = ${DIRECT_XGB_MIN_CHILD_WEIGHT:-20}"
echo "  LAMBDA               = ${DIRECT_XGB_LAMBDA:-1}"
echo "  NUM_BOOST_ROUND      = ${DIRECT_XGB_NUM_BOOST_ROUND:-2000}"
echo "  EARLY_STOPPING       = ${DIRECT_XGB_EARLY_STOPPING_ROUNDS:-50}"
echo "  CV_FOLDS             = ${DIRECT_XGB_CV_FOLDS:-4}"
echo "  THREADS              = ${DIRECT_XGB_THREADS:-$(nproc)}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo ""

cmd=(
    python3 "${SCRIPT_DIR}/train_gradcpt_flanker_direct_xgb_proxy.py"
    --ses-ea-dir "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}"
    --all-scores-file "${ALL_SCORES}"
    --task-score-file "${LOCAL_ETM_TASK}"
    --fine-tuned-score-file "${FINE_TUNED_SCORES}"
    --ea-query "${LOCAL_SCRAP}/ea_query.csv"
    --main-survey "${LOCAL_SCRAP}/main_survey_features.csv"
    --bhp-survey "${LOCAL_SCRAP}/bhp_survey_features.csv"
    --area-ses "${LOCAL_SCRAP}/area_ses.csv"
    --metadata "${METADATA}"
    --output-dir "${LOCAL_OUT}"
    --workspace-output-dir "${WORKSPACE_OUT}"
    --workspace-scrap-dir "${WORKSPACE_DIAG}"
)
if [[ "${FORCE}" -eq 1 ]]; then
    cmd+=(--force)
fi
if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
    cmd+=(--stage-aggregate)
fi
if [[ "${NO_WRITE_REGENIE_INPUTS}" -eq 1 ]]; then
    cmd+=(--no-write-regenie-inputs)
fi

"${cmd[@]}"

echo ""
echo "=== Direct GradCPT/Flanker XGBoost proxy command complete ==="
