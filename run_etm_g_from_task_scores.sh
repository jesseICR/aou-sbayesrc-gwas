#!/bin/bash
# run_etm_g_from_task_scores.sh - Build ETM general factors from task scores.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage: bash run_etm_g_from_task_scores.sh [OPTIONS]

Builds three-domain and, when Emotional Recognition is available, four-domain
ETM general cognitive/performance factors from the already-created ETM
task-score table. This command does not query ETM tables or run GWAS.

Options:
  --task-score-file PATH       Override recommended task-score wide TSV.
  --all-scores-file PATH       Override SES-EA proxy all_scores.tsv.
  --output-dir PATH            Override local ETM-g output directory.
  --stage-aggregate            Stage aggregate diagnostics to workspace bucket scrap.
  --force                      Recompute and overwrite existing ETM-g outputs.
  --flanker-input VALUE        auto, flanker_efficiency_z_age_sex, or flanker_perf_z_age_sex.
  --min-complete-case-n N      Minimum complete-case reference N (default: 500).
  --force-three-domain-g       Write etm_g_z_forced if acceptance criteria fail.
  -h, --help                   Show this help.
EOF
}

TASK_SCORE_FILE=""
ALL_SCORES_FILE=""
OUTPUT_DIR=""
STAGE_AGGREGATE=0
FORCE=0
FLANKER_INPUT="auto"
MIN_COMPLETE_CASE_N=500
FORCE_THREE_DOMAIN_G=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --task-score-file)
            TASK_SCORE_FILE="${2:?--task-score-file requires PATH}"
            shift 2
            ;;
        --all-scores-file)
            ALL_SCORES_FILE="${2:?--all-scores-file requires PATH}"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="${2:?--output-dir requires PATH}"
            shift 2
            ;;
        --stage-aggregate)
            STAGE_AGGREGATE=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --flanker-input)
            FLANKER_INPUT="${2:?--flanker-input requires a value}"
            shift 2
            ;;
        --min-complete-case-n)
            MIN_COMPLETE_CASE_N="${2:?--min-complete-case-n requires N}"
            shift 2
            ;;
        --force-three-domain-g)
            FORCE_THREE_DOMAIN_G=1
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

TASK_SCORE_FILE="${TASK_SCORE_FILE:-${LOCAL_REGENIE_DIR}/ses_ea_proxy_scrap/etm_cog_task_factors/etm_cog_task_factors_recommended_wide.tsv}"
ALL_SCORES_FILE="${ALL_SCORES_FILE:-${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/all_scores.tsv}"
BASE_COVAR_FILE="${ALL_SCORES_FILE%/*}/base_covar.txt"
OUTPUT_DIR="${OUTPUT_DIR:-${LOCAL_REGENIE_DIR}/ses_ea_proxy_scrap/etm_cog_task_factors/etm_general_factor}"
WORKSPACE_SCRAP="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap/etm_cog_task_factors/etm_general_factor"

if [[ ! -s "${TASK_SCORE_FILE}" ]]; then
    echo "ERROR: missing task-score file ${TASK_SCORE_FILE}" >&2
    echo "  Run bash run_etm_cog_task_factors.sh --stage-aggregate first." >&2
    exit 1
fi
if [[ ! -s "${ALL_SCORES_FILE}" ]]; then
    echo "ERROR: missing ${ALL_SCORES_FILE}" >&2
    echo "  Run bash run_ses_ea_proxy_gwas.sh --setup-only first." >&2
    exit 1
fi
if [[ ! -s "${BASE_COVAR_FILE}" ]]; then
    echo "ERROR: missing ${BASE_COVAR_FILE}" >&2
    exit 1
fi

mkdir -p "${SCRIPT_DIR}/logs" "${OUTPUT_DIR}"
if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
    mkdir -p "${WORKSPACE_SCRAP}"
fi

LOG_FILE="${SCRIPT_DIR}/logs/run_etm_g_from_task_scores_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "AoU ETM general factor scoring — $(date)"
echo "============================================"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
echo "  Task score file      = ${TASK_SCORE_FILE}"
echo "  All scores file      = ${ALL_SCORES_FILE}"
echo "  Base covar file      = ${BASE_COVAR_FILE}"
echo "  Local output dir     = ${OUTPUT_DIR}"
echo "  Workspace scrap dir  = ${WORKSPACE_SCRAP}"
echo "  FLANKER_INPUT        = ${FLANKER_INPUT}"
echo "  MIN_COMPLETE_CASE_N  = ${MIN_COMPLETE_CASE_N}"
echo "  FORCE                = ${FORCE}"
echo "  FORCE_THREE_DOMAIN_G = ${FORCE_THREE_DOMAIN_G}"
echo "  STAGE_AGGREGATE      = ${STAGE_AGGREGATE}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo ""

cmd=(
    python3 "${SCRIPT_DIR}/score_etm_general_factor.py"
    --task-score-file "${TASK_SCORE_FILE}"
    --all-scores-file "${ALL_SCORES_FILE}"
    --base-covar-file "${BASE_COVAR_FILE}"
    --output-dir "${OUTPUT_DIR}"
    --workspace-scrap-dir "${WORKSPACE_SCRAP}"
    --flanker-input "${FLANKER_INPUT}"
    --min-complete-case-n "${MIN_COMPLETE_CASE_N}"
)
if [[ "${FORCE}" -eq 1 ]]; then
    cmd+=(--force)
fi
if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
    cmd+=(--stage-aggregate)
fi
if [[ "${FORCE_THREE_DOMAIN_G}" -eq 1 ]]; then
    cmd+=(--force-three-domain-g)
fi

"${cmd[@]}"

echo ""
echo "=== ETM general factor command complete ==="
