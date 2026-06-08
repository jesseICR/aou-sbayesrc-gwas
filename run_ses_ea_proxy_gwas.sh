#!/bin/bash
# run_ses_ea_proxy_gwas.sh - Build AoU ses_ea_proxy scores; GWAS is opt-in.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR

usage() {
    cat <<'EOF'
Usage: bash run_ses_ea_proxy_gwas.sh [OPTIONS]

Builds primary ses_ea_proxy XGBoost scores and REGENIE input files. By default
this command stops after setup/scoring and does not submit REGENIE.

Options:
  --run-gwas         After setup, submit the continuous-trait REGENIE GWAS.
  --regenie-only     Skip setup and run REGENIE against existing inputs.
  --setup-only       Explicitly stop after setup/scoring (default).
  --chroms LIST      Chromosomes for REGENIE if --run-gwas/--regenie-only is used.
  --apply-rint       Apply rank-inverse normal transform for REGENIE (default).
  --no-apply-rint    Disable rank-inverse normal transform for REGENIE.
  -h, --help         Show this help.
EOF
}

RUN_SETUP=1
RUN_REGENIE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --setup-only)
            RUN_SETUP=1
            RUN_REGENIE=0
            shift
            ;;
        --run-gwas)
            RUN_REGENIE=1
            shift
            ;;
        --regenie-only)
            RUN_SETUP=0
            RUN_REGENIE=1
            shift
            ;;
        --chroms)
            export REGENIE_CHROMS="${2:?--chroms requires a value}"
            shift 2
            ;;
        --apply-rint)
            export REGENIE_APPLY_RINT=1
            shift
            ;;
        --no-apply-rint)
            export REGENIE_APPLY_RINT=0
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

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set — are you running inside an AoU Verily Jupyter session?}"
: "${WORKSPACE_CDR:?WORKSPACE_CDR not set — are you running inside an AoU Verily Jupyter session?}"

export WORKSPACE_BUCKET_MOUNT="/home/jupyter/workspace/workspace-bucket"
WORKSPACE_BUCKET_URI="gs://$(mount | awk '/ on \/home\/jupyter\/workspace\/workspace-bucket /{print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    exit 1
fi
export WORKSPACE_BUCKET_URI

export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/sbayesrc_genotypes"
export DX_EUROPEANS_DIR="${DX_OUTPUT_DIR}/europeans"
export DX_KINSHIP_DIR="${DX_OUTPUT_DIR}/kinship"
export DX_PCA_EUR_DIR="${DX_OUTPUT_DIR}/pca_eur"
export DX_GENETIC_SEX_DIR="${DX_OUTPUT_DIR}/genetic_sex"
export DX_SAMPLE_QC_DIR="${DX_OUTPUT_DIR}/sample_qc"
export DX_GWAS_GENOTYPES_DIR="${DX_OUTPUT_DIR}/gwas_genotypes"
export DX_GWAS_STEP1_BFILE_DIR="${DX_GWAS_GENOTYPES_DIR}/step1_direct"
export DX_GWAS_STEP2_PFILE_DIR="${DX_GWAS_GENOTYPES_DIR}/step2_wgs_pfiles"
export DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
export DX_REGENIE_OUTPUT_DIR="${DX_OUTPUT_DIR}/regenie_output"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"

export DX_GWAS_GENOTYPES_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/gwas_genotypes"
export DX_GWAS_STEP1_BFILE_URI="${DX_GWAS_GENOTYPES_URI}/step1_direct"
export DX_GWAS_STEP2_PFILE_URI="${DX_GWAS_GENOTYPES_URI}/step2_wgs_pfiles"
export DX_REGENIE_INPUT_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/regenie_input"
export DX_REGENIE_OUTPUT_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/regenie_output"

export SES_EA_PROXY_GWAS_INPUT_NAME="${SES_EA_PROXY_GWAS_INPUT_NAME:-ses_ea_proxy}"
export SES_EA_PROXY_GWAS_OUTPUT_NAME="${SES_EA_PROXY_GWAS_OUTPUT_NAME:-ses_ea_proxy}"
export DX_SES_EA_PROXY_REGENIE_INPUT_DIR="${DX_REGENIE_INPUT_DIR}/${SES_EA_PROXY_GWAS_INPUT_NAME}"
export DX_SES_EA_PROXY_REGENIE_INPUT_URI="${DX_REGENIE_INPUT_URI}/${SES_EA_PROXY_GWAS_INPUT_NAME}"
export LOCAL_REGENIE_DIR="${SCRIPT_DIR}/data/regenie"

export IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE="${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE:-3}"
export SES_EA_PROXY_N_PCS="${SES_EA_PROXY_N_PCS:-10}"
export SES_EA_PROXY_OUTER_FOLDS="${SES_EA_PROXY_OUTER_FOLDS:-5}"
export SES_EA_PROXY_SEED="${SES_EA_PROXY_SEED:-2026}"
export SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD="${SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD:-0.0441941}"
export GWAS_MIN_AGE_AT_SURVEY="${GWAS_MIN_AGE_AT_SURVEY:-26}"
export REGENIE_APPLY_RINT="${REGENIE_APPLY_RINT:-1}"
export REGENIE_STEP1_BLOCK_SIZE="${REGENIE_STEP1_BLOCK_SIZE:-1000}"
export REGENIE_STEP2_BLOCK_SIZE="${REGENIE_STEP2_BLOCK_SIZE:-200}"
export REGENIE_CHROMS="${REGENIE_CHROMS:-1-22}"
export REGENIE_PHENO_COL="${REGENIE_PHENO_COL:-ses_ea_proxy_z}"
if [[ -z "${REGENIE_COVAR_COLS:-}" ]]; then
    REGENIE_COVAR_COLS="yob_c,sex_c,yob_c_sex_c_inter"
    for i in $(seq 1 "${SES_EA_PROXY_N_PCS}"); do
        REGENIE_COVAR_COLS+=",PC${i}_AVG"
    done
    export REGENIE_COVAR_COLS
fi

export DSUB_PROVIDER="${DSUB_PROVIDER:-google-batch}"
export DSUB_REGION="${DSUB_REGION:-us-central1}"
export DSUB_NETWORK="projects/${GOOGLE_PROJECT}/global/networks/network"
export DSUB_SUBNETWORK="projects/${GOOGLE_PROJECT}/regions/${DSUB_REGION}/subnetworks/subnetwork"
export DSUB_IMAGE="${DSUB_IMAGE:-marketplace.gcr.io/google/ubuntu2204}"
DSUB_PET_SA="$(gcloud config get-value account 2>/dev/null || true)"
if [[ -z "${DSUB_PET_SA}" ]]; then
    echo "ERROR: could not determine the pod's pet service account via 'gcloud config get-value account'." >&2
    exit 1
fi
export DSUB_PET_SA
export DSUB_BOOT_DISK_SIZE="${DSUB_BOOT_DISK_SIZE:-50}"
export REGENIE_STEP1_DSUB_MIN_CORES="${REGENIE_STEP1_DSUB_MIN_CORES:-16}"
export REGENIE_STEP1_DSUB_MIN_RAM="${REGENIE_STEP1_DSUB_MIN_RAM:-64}"
export REGENIE_STEP1_DSUB_DISK_SIZE="${REGENIE_STEP1_DSUB_DISK_SIZE:-300}"
export REGENIE_STEP1_DSUB_DISK_TYPE="${REGENIE_STEP1_DSUB_DISK_TYPE:-pd-ssd}"
export REGENIE_STEP2_DSUB_MIN_CORES="${REGENIE_STEP2_DSUB_MIN_CORES:-8}"
export REGENIE_STEP2_DSUB_MIN_RAM="${REGENIE_STEP2_DSUB_MIN_RAM:-32}"
export REGENIE_STEP2_DSUB_DISK_SIZE="${REGENIE_STEP2_DSUB_DISK_SIZE:-300}"
export REGENIE_STEP2_DSUB_DISK_TYPE="${REGENIE_STEP2_DSUB_DISK_TYPE:-pd-ssd}"

export DSUB_BIN_URI="${WORKSPACE_BUCKET_URI}/bin"
export DSUB_REGENIE_BUNDLE_URI="${DSUB_BIN_URI}/regenie_bundle"
export DSUB_LOG_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/logs/dsub"

if [[ "${RUN_REGENIE}" -eq 1 ]]; then
    if [[ -n "${REGENIE:-}" && -x "${REGENIE}" ]]; then
        export REGENIE
    elif command -v regenie >/dev/null 2>&1; then
        export REGENIE="$(command -v regenie)"
    elif [[ -x /opt/workbench-tools/binaries/bin/regenie ]]; then
        export REGENIE="/opt/workbench-tools/binaries/bin/regenie"
    else
        echo "ERROR: could not find regenie; set REGENIE to an executable path." >&2
        exit 1
    fi
else
    export REGENIE="${REGENIE:-}"
fi

mkdir -p "${SCRIPT_DIR}/logs" "${LOCAL_REGENIE_DIR}" "${DX_REGENIE_INPUT_DIR}" \
    "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}" "${DX_REGENIE_OUTPUT_DIR}" "${DX_LOGS_DIR}"

LOG_FILE="${SCRIPT_DIR}/logs/run_ses_ea_proxy_gwas_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "AoU ses_ea_proxy setup — $(date)"
echo "============================================"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
echo "  WORKSPACE_CDR        = ${WORKSPACE_CDR}"
echo "  WORKSPACE_MHWB_CDR   = ${WORKSPACE_MHWB_CDR:-${WORKSPACE_CDR%%.*}.C_V8_R2_offcycle_mhwb}"
echo "  proxy input          = ${SES_EA_PROXY_GWAS_INPUT_NAME}"
echo "  proxy output         = ${SES_EA_PROXY_GWAS_OUTPUT_NAME}"
echo "  Samples              = classified Europeans, confirmed genetic sex, excluding identical components size >= ${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}, EA label age >= ${GWAS_MIN_AGE_AT_SURVEY}"
echo "  Extended design      = ${SES_EA_PROXY_OUTER_FOLDS} OOF folds over fit_pca_iids, sixth model applied to classified-EUR extras"
echo "  Applied model holdout= exclude fit_pca relatives of applied samples at KING >= ${SES_EA_PROXY_FINAL_KINSHIP_THRESHOLD}"
echo "  Covariates for GWAS  = ${REGENIE_COVAR_COLS}"
echo "  RUN_REGENIE          = ${RUN_REGENIE}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo ""

if [[ "${RUN_SETUP}" -eq 1 ]]; then
    echo "=== Setup ses_ea_proxy scores and GWAS inputs ==="
    bash "${SCRIPT_DIR}/setup_ses_ea_proxy_gwas.sh"
else
    echo "=== Setup ses_ea_proxy scores and GWAS inputs ==="
    echo "  Skipping setup because --regenie-only was specified."
fi

if [[ "${RUN_REGENIE}" -eq 1 ]]; then
    echo ""
    echo "=== Run ses_ea_proxy GWAS with REGENIE ==="
    bash "${SCRIPT_DIR}/run_continuous_regenie_gwas.sh" \
        "${SES_EA_PROXY_GWAS_INPUT_NAME}" "${SES_EA_PROXY_GWAS_OUTPUT_NAME}"
else
    echo ""
    echo "=== Run ses_ea_proxy GWAS with REGENIE ==="
    echo "  Skipping REGENIE. Review OOF performance and covariate correlations first."
fi

echo ""
echo "=== ses_ea_proxy command complete ==="
