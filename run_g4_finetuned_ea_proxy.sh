#!/bin/bash
# run_g4_finetuned_ea_proxy.sh - Fine-tune SES-EA proxy boosters toward ETM-derived targets.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR
source "${SCRIPT_DIR}/aou_downstream_env.sh"

usage() {
    cat <<'EOF'
Usage: bash run_g4_finetuned_ea_proxy.sh [OPTIONS]

Continues training the saved SES-EA proxy XGBoost boosters toward an ETM-
derived cognitive target. The default target is the four-domain ETM-g4 score;
the final selected phenotype path uses --target gradcpt-flanker-mean through
run_gradcpt_flanker_finetuned_ea_proxy.sh. This command does not query BigQuery,
run GWAS, commit, or push.

Options:
  --force                    Recompute and overwrite this fine-tuning output dir.
  --stage-aggregate          Stage aggregate diagnostics to ses_ea_proxy/scrap.
  --reuse-feature-extracts   Require existing local survey extract CSVs.
  --target MODE              strong-task-g4 (default), all-etm-g4, or gradcpt-flanker-mean.
  --no-write-regenie-inputs  Do not write phen/covar/training_iids files.
  -h, --help                 Show this help.
EOF
}

FORCE=0
STAGE_AGGREGATE=0
REUSE_FEATURE_EXTRACTS=0
TARGET_MODE="strong-task-g4"
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
        --reuse-feature-extracts)
            REUSE_FEATURE_EXTRACTS=1
            shift
            ;;
        --target)
            TARGET_MODE="${2:?--target requires MODE}"
            shift 2
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

case "${TARGET_MODE}" in
    strong-task-g4|all-etm-g4|gradcpt-flanker-mean) ;;
    *)
        echo "ERROR: --target must be strong-task-g4, all-etm-g4, or gradcpt-flanker-mean" >&2
        exit 1
        ;;
esac

export WORKSPACE_BUCKET_MOUNT="/home/jupyter/workspace/workspace-bucket"
WORKSPACE_BUCKET_URI="gs://$(mount | awk '/ on \/home\/jupyter\/workspace\/workspace-bucket /{print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    exit 1
fi
export WORKSPACE_BUCKET_URI

export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/sbayesrc_genotypes"
export DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
export SES_EA_PROXY_GWAS_INPUT_NAME="${SES_EA_PROXY_GWAS_INPUT_NAME:-ses_ea_proxy}"
export DX_SES_EA_PROXY_REGENIE_INPUT_DIR="${DX_REGENIE_INPUT_DIR}/${SES_EA_PROXY_GWAS_INPUT_NAME}"
export LOCAL_REGENIE_DIR="${SCRIPT_DIR}/data/regenie"

LOCAL_SCRAP="${LOCAL_REGENIE_DIR}/ses_ea_proxy_scrap"
LOCAL_ETM_TASK="${LOCAL_SCRAP}/etm_cog_task_factors/etm_cog_task_factors_recommended_wide.tsv"
LOCAL_ETM_G="${LOCAL_SCRAP}/etm_cog_task_factors/etm_general_factor/etm_general_factor_scores_wide.tsv"
WORKSPACE_EXTRACT_SCRAP="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap"

case "${TARGET_MODE}" in
    strong-task-g4)
        DEFAULT_OUTPUT_NAME="g4_finetuned_ea_proxy"
        ;;
    all-etm-g4)
        DEFAULT_OUTPUT_NAME="g4_finetuned_ea_proxy_all_etm_g4"
        ;;
    gradcpt-flanker-mean)
        DEFAULT_OUTPUT_NAME="gradcpt_flanker_finetuned_ea_proxy"
        ;;
esac
if [[ "${SES_EA_PROXY_GWAS_INPUT_NAME}" != "ses_ea_proxy" ]]; then
    DEFAULT_OUTPUT_NAME="${DEFAULT_OUTPUT_NAME}_${SES_EA_PROXY_GWAS_INPUT_NAME}"
fi
export G4_FINETUNED_OUTPUT_NAME="${G4_FINETUNED_OUTPUT_NAME:-${DEFAULT_OUTPUT_NAME}}"

LOCAL_OUT="${LOCAL_SCRAP}/${G4_FINETUNED_OUTPUT_NAME}"
WORKSPACE_OUT="${DX_REGENIE_INPUT_DIR}/${G4_FINETUNED_OUTPUT_NAME}"
WORKSPACE_DIAG="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap/${G4_FINETUNED_OUTPUT_NAME}"

ALL_SCORES="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/all_scores.tsv"
BASE_COVAR="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/base_covar.txt"
COVAR="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/covar.txt"
MODEL_MANIFEST="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_model_manifest.tsv"
FEATURE_COLUMNS="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_feature_columns.json"
METADATA="${SCRIPT_DIR}/data/aou_metadata/aou_ds_survey_question_concepts.tsv"

required=(
    "${ALL_SCORES}"
    "${BASE_COVAR}"
    "${COVAR}"
    "${MODEL_MANIFEST}"
    "${FEATURE_COLUMNS}"
    "${LOCAL_ETM_TASK}"
    "${METADATA}"
)
if [[ "${TARGET_MODE}" == "strong-task-g4" || "${TARGET_MODE}" == "all-etm-g4" ]]; then
    required+=("${LOCAL_ETM_G}")
fi
for fold in 0 1 2 3 4; do
    required+=("${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_models/fold_${fold}.json")
done
required+=("${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/xgboost_models/final_model.json")
for f in "${required[@]}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing required input ${f}" >&2
        if [[ "${f}" == "${LOCAL_ETM_G}" ]]; then
            echo "  Run bash run_etm_g_from_task_scores.sh --stage-aggregate first." >&2
        elif [[ "${f}" == "${ALL_SCORES}" ]]; then
            echo "  Run bash run_ses_ea_proxy_gwas.sh --setup-only first." >&2
        fi
        exit 1
    fi
done

mkdir -p "${SCRIPT_DIR}/logs" "${LOCAL_SCRAP}" "${LOCAL_OUT}"
if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
    mkdir -p "${WORKSPACE_DIAG}"
fi

LOG_FILE="${SCRIPT_DIR}/logs/run_g4_finetuned_ea_proxy_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

ensure_extract() {
    local name="$1"
    local local_path="${LOCAL_SCRAP}/${name}"
    local workspace_path="${WORKSPACE_EXTRACT_SCRAP}/${name}"
    local stem="${name%.csv}"
    if [[ -s "${local_path}" ]]; then
        return 0
    fi
    if [[ "${REUSE_FEATURE_EXTRACTS}" -eq 1 ]]; then
        echo "ERROR: missing local feature extract ${local_path} with --reuse-feature-extracts" >&2
        exit 1
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

echo "============================================"
echo "AoU fine-tuned SES-EA proxy — $(date)"
echo "============================================"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
echo "  SES-EA input dir     = ${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}"
echo "  ETM task score file  = ${LOCAL_ETM_TASK}"
if [[ "${TARGET_MODE}" == "strong-task-g4" || "${TARGET_MODE}" == "all-etm-g4" ]]; then
    echo "  ETM-g score file     = ${LOCAL_ETM_G}"
else
    echo "  ETM-g score file     = not required for ${TARGET_MODE}"
fi
echo "  Local output dir     = ${LOCAL_OUT}"
echo "  Workspace output dir = ${WORKSPACE_OUT}"
echo "  Workspace diag dir   = ${WORKSPACE_DIAG}"
echo "  TARGET_MODE          = ${TARGET_MODE}"
echo "  FORCE                = ${FORCE}"
echo "  STAGE_AGGREGATE      = ${STAGE_AGGREGATE}"
echo "  REUSE_EXTRACTS       = ${REUSE_FEATURE_EXTRACTS}"
echo "  WRITE_REGENIE_INPUTS = $((1 - NO_WRITE_REGENIE_INPUTS))"
echo "  ETA                  = ${G4_FINETUNE_ETA:-0.01}"
echo "  MAX_DEPTH            = ${G4_FINETUNE_MAX_DEPTH:-3}"
echo "  MIN_CHILD_WEIGHT     = ${G4_FINETUNE_MIN_CHILD_WEIGHT:-10}"
echo "  LAMBDA               = ${G4_FINETUNE_LAMBDA:-2.0}"
echo "  MAX_ROUNDS           = ${G4_FINETUNE_MAX_ROUNDS:-500}"
echo "  EARLY_STOPPING       = ${G4_FINETUNE_EARLY_STOPPING_ROUNDS:-25}"
echo "  VALID_FRACTION       = ${G4_FINETUNE_VALID_FRACTION:-0.20}"
echo "  MIN_TRAIN_SAMPLES    = ${G4_FINETUNE_MIN_TRAIN_SAMPLES:-1000}"
echo "  SEED                 = ${G4_FINETUNE_SEED:-2026}"
echo "  THREADS              = ${G4_FINETUNE_THREADS:-$(nproc)}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo ""

cmd=(
    python3 "${SCRIPT_DIR}/fine_tune_g4_ea_proxy.py"
    --ses-ea-dir "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}"
    --all-scores-file "${ALL_SCORES}"
    --task-score-file "${LOCAL_ETM_TASK}"
    --etm-g-file "${LOCAL_ETM_G}"
    --ea-query "${LOCAL_SCRAP}/ea_query.csv"
    --main-survey "${LOCAL_SCRAP}/main_survey_features.csv"
    --bhp-survey "${LOCAL_SCRAP}/bhp_survey_features.csv"
    --area-ses "${LOCAL_SCRAP}/area_ses.csv"
    --metadata "${METADATA}"
    --output-dir "${LOCAL_OUT}"
    --workspace-output-dir "${WORKSPACE_OUT}"
    --workspace-scrap-dir "${WORKSPACE_DIAG}"
    --target "${TARGET_MODE}"
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
echo "=== Fine-tuned SES-EA proxy command complete ==="
