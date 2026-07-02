#!/bin/bash
# run_g_ea_proxy_v9_pipeline.sh - Build and run the cdrv9 g-EA proxy GWAS.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR
source "${SCRIPT_DIR}/aou_downstream_env.sh"

usage() {
    cat <<'EOF'
Usage: bash run_g_ea_proxy_v9_pipeline.sh [OPTIONS]

Builds the v9 SES-EA/ETM/GradCPT-Flanker proxy chain and runs the final
continuous-trait REGENIE GWAS for:

  g_ea_proxy_sbayesrc7m

Options:
  --setup-only       Run proxy setup/scoring and write final REGENIE inputs,
                     but do not submit REGENIE.
  --preflight-only   Validate existing final-source inputs and write final
                     REGENIE inputs; do not run upstream setup or REGENIE.
  --skip-setup       Skip upstream setup/scoring and run the final GWAS stage
                     from existing direct-XGB source outputs.
  --smoke            Run upstream setup if needed, then chr22 smoke GWAS.
  --chroms LIST      Chromosomes for the final GWAS, e.g. 22, 1,2,3, or 1-22.
  --force-final      Pass --force to the final calibration/GWAS runner.
  -h, --help         Show this help.
EOF
}

RUN_SETUP=1
RUN_REGENIE=1
PREFLIGHT_ONLY=0
SMOKE=0
FORCE_FINAL=0
FINAL_CHROMS="${REGENIE_CHROMS:-1-22}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --setup-only)
            RUN_SETUP=1
            RUN_REGENIE=0
            PREFLIGHT_ONLY=0
            shift
            ;;
        --preflight-only)
            RUN_SETUP=0
            RUN_REGENIE=0
            PREFLIGHT_ONLY=1
            shift
            ;;
        --skip-setup)
            RUN_SETUP=0
            shift
            ;;
        --smoke)
            SMOKE=1
            FINAL_CHROMS="22"
            shift
            ;;
        --chroms)
            FINAL_CHROMS="${2:?--chroms requires a value}"
            shift 2
            ;;
        --force-final)
            FORCE_FINAL=1
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

export AOU_DATA_VERSION="${AOU_DATA_VERSION:-v9}"
export SES_EA_PROXY_GWAS_INPUT_NAME="${SES_EA_PROXY_GWAS_INPUT_NAME:-ses_ea_proxy_v2_kinholdout}"
export SES_EA_PROXY_GWAS_OUTPUT_NAME="${SES_EA_PROXY_GWAS_OUTPUT_NAME:-${SES_EA_PROXY_GWAS_INPUT_NAME}}"
export FACTOR18_GWAS_SOURCE_NAME="${FACTOR18_GWAS_SOURCE_NAME:-gradcpt_flanker_direct_xgb_proxy_${SES_EA_PROXY_GWAS_INPUT_NAME}}"
export FACTOR18_GWAS_INPUT_NAME="${FACTOR18_GWAS_INPUT_NAME:-g_ea_proxy_sbayesrc7m}"
export FACTOR18_GWAS_OUTPUT_NAME="${FACTOR18_GWAS_OUTPUT_NAME:-g_ea_proxy_sbayesrc7m_gwas}"

mkdir -p "${SCRIPT_DIR}/logs"
LOG_FILE="${SCRIPT_DIR}/logs/run_g_ea_proxy_v9_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "================================================"
echo "AoU cdrv9 g-EA proxy GWAS pipeline - $(date)"
echo "================================================"
echo "  GOOGLE_PROJECT        = ${GOOGLE_PROJECT}"
echo "  WORKSPACE_CDR         = ${WORKSPACE_CDR}"
echo "  WORKSPACE_MHWB_CDR    = ${WORKSPACE_MHWB_CDR}"
echo "  WORKSPACE_ETM_CDR     = ${WORKSPACE_ETM_CDR}"
echo "  WORKSPACE_BUCKET_URI  = ${WORKSPACE_BUCKET_URI}"
echo "  SES-EA input/output   = ${SES_EA_PROXY_GWAS_INPUT_NAME}"
echo "  Final source          = ${FACTOR18_GWAS_SOURCE_NAME}"
echo "  Final input           = ${FACTOR18_GWAS_INPUT_NAME}"
echo "  Final output          = ${FACTOR18_GWAS_OUTPUT_NAME}"
echo "  Run setup             = ${RUN_SETUP}"
echo "  Run REGENIE           = ${RUN_REGENIE}"
echo "  Preflight only        = ${PREFLIGHT_ONLY}"
echo "  Smoke                 = ${SMOKE}"
echo "  Final chromosomes     = ${FINAL_CHROMS}"
echo "  Force final           = ${FORCE_FINAL}"
echo "  LOG_FILE              = ${LOG_FILE}"
echo ""

if [[ "${RUN_SETUP}" -eq 1 ]]; then
    echo "=== 1. SES-EA proxy setup ==="
    bash "${SCRIPT_DIR}/run_ses_ea_proxy_gwas.sh" --setup-only

    echo ""
    echo "=== 2. ETM cognitive task factors ==="
    bash "${SCRIPT_DIR}/run_etm_cog_task_factors.sh" --stage-aggregate

    echo ""
    echo "=== 3. GradCPT/Flanker fine-tuned EA proxy ==="
    bash "${SCRIPT_DIR}/run_gradcpt_flanker_finetuned_ea_proxy.sh" --stage-aggregate

    echo ""
    echo "=== 4. Direct GradCPT/Flanker XGBoost proxy ==="
    bash "${SCRIPT_DIR}/run_gradcpt_flanker_direct_xgb_proxy.sh" --stage-aggregate
else
    echo "=== Upstream proxy setup ==="
    echo "  Skipping upstream setup/scoring."
fi

final_args=()
if [[ "${FORCE_FINAL}" -eq 1 ]]; then
    final_args+=(--force)
fi

echo ""
if [[ "${PREFLIGHT_ONLY}" -eq 1 || "${RUN_REGENIE}" -eq 0 ]]; then
    echo "=== 5. Final g-EA proxy REGENIE input preflight ==="
    bash "${SCRIPT_DIR}/run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas.sh" \
        --preflight-only "${final_args[@]}"
elif [[ "${SMOKE}" -eq 1 ]]; then
    echo "=== 5. Final g-EA proxy chr22 smoke GWAS ==="
    bash "${SCRIPT_DIR}/run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas.sh" \
        --smoke "${final_args[@]}"
else
    echo "=== 5. Final g-EA proxy GWAS ==="
    bash "${SCRIPT_DIR}/run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas.sh" \
        --chroms "${FINAL_CHROMS}" "${final_args[@]}"
fi

echo ""
echo "=== cdrv9 g-EA proxy pipeline command complete ==="
