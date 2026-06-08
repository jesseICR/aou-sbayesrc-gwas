#!/bin/bash
# run_etm_cog_task_factors.sh - Build ETM cognitive task factor scores.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage: bash run_etm_cog_task_factors.sh [OPTIONS]

Builds Flanker, GradCPT, Delay Discounting, and Emotional Recognition cognitive
task scores for the existing ses_ea_proxy phenotype cohort. This command does
not run GWAS.

Options:
  --reuse-extracts              Reuse existing local ETM extract; fail if absent.
  --force                       Re-query ETM and overwrite outputs.
  --stage-aggregate             Stage aggregate diagnostics to workspace bucket scrap.
  --make-etm-g                  Run downstream ETM general-factor scoring afterward.
  --etm-dataset PROJECT.DATASET Override WORKSPACE_ETM_CDR.
  --bq-temp-dataset DATASET     Existing writable BigQuery dataset for temp tables.
  -h, --help                    Show this help.
EOF
}

REUSE_EXTRACTS=0
FORCE=0
STAGE_AGGREGATE=0
MAKE_ETM_G=0
ETM_DATASET_OVERRIDE=""
BQ_TMP_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reuse-extracts)
            REUSE_EXTRACTS=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --stage-aggregate)
            STAGE_AGGREGATE=1
            shift
            ;;
        --make-etm-g)
            MAKE_ETM_G=1
            shift
            ;;
        --etm-dataset)
            ETM_DATASET_OVERRIDE="${2:?--etm-dataset requires PROJECT.DATASET}"
            shift 2
            ;;
        --bq-temp-dataset)
            BQ_TMP_OVERRIDE="${2:?--bq-temp-dataset requires a dataset name}"
            shift 2
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
export DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
export SES_EA_PROXY_GWAS_INPUT_NAME="${SES_EA_PROXY_GWAS_INPUT_NAME:-ses_ea_proxy}"
export DX_SES_EA_PROXY_REGENIE_INPUT_DIR="${DX_REGENIE_INPUT_DIR}/${SES_EA_PROXY_GWAS_INPUT_NAME}"
export LOCAL_REGENIE_DIR="${SCRIPT_DIR}/data/regenie"

if [[ -n "${ETM_DATASET_OVERRIDE}" ]]; then
    WORKSPACE_ETM_CDR="${ETM_DATASET_OVERRIDE}"
else
    WORKSPACE_ETM_CDR="${WORKSPACE_ETM_CDR:-${WORKSPACE_CDR%%.*}.C_V8_R2_offcycle_etm}"
fi
export WORKSPACE_ETM_CDR

choose_bq_tmp_dataset() {
    local requested="${1:-}" candidate
    if [[ -n "${requested}" ]]; then
        printf '%s\n' "${requested}"
        return 0
    fi
    for candidate in "${ETM_COG_BQ_TMP_DATASET:-}" sbayesrc_tmp high_quality_cohort dataset_test2; do
        if [[ -n "${candidate}" ]] && bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${candidate}" >/dev/null 2>&1; then
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
    echo "  Set ETM_COG_BQ_TMP_DATASET or pass --bq-temp-dataset to an existing writable dataset." >&2
    return 1
}

BQ_TMP_DATASET="$(choose_bq_tmp_dataset "${BQ_TMP_OVERRIDE}")"
export BQ_TMP_DATASET

if [[ ! -s "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/all_scores.tsv" ]]; then
    echo "ERROR: missing ${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/all_scores.tsv" >&2
    echo "  Run bash run_ses_ea_proxy_gwas.sh --setup-only first." >&2
    exit 1
fi
if [[ ! -s "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/base_covar.txt" ]]; then
    echo "ERROR: missing ${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/base_covar.txt" >&2
    exit 1
fi

LOCAL_OUT="${LOCAL_REGENIE_DIR}/ses_ea_proxy_scrap/etm_cog_task_factors"
WORKSPACE_SCRAP="${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}/scrap/etm_cog_task_factors"
mkdir -p "${SCRIPT_DIR}/logs" "${LOCAL_OUT}" "${WORKSPACE_SCRAP}"

LOG_FILE="${SCRIPT_DIR}/logs/run_etm_cog_task_factors_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "AoU ETM cognitive task factor scoring — $(date)"
echo "============================================"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
echo "  WORKSPACE_CDR        = ${WORKSPACE_CDR}"
echo "  WORKSPACE_ETM_CDR    = ${WORKSPACE_ETM_CDR}"
echo "  BQ_TMP_DATASET       = ${BQ_TMP_DATASET}"
echo "  SES-EA input dir     = ${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}"
echo "  Local output dir     = ${LOCAL_OUT}"
echo "  Workspace scrap dir  = ${WORKSPACE_SCRAP}"
echo "  REUSE_EXTRACTS       = ${REUSE_EXTRACTS}"
echo "  FORCE                = ${FORCE}"
echo "  STAGE_AGGREGATE      = ${STAGE_AGGREGATE}"
echo "  MAKE_ETM_G           = ${MAKE_ETM_G}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo ""

cmd=(
    python3 "${SCRIPT_DIR}/score_etm_cog_task_factors.py"
    --etm-dataset "${WORKSPACE_ETM_CDR}"
    --workspace-cdr "${WORKSPACE_CDR}"
    --bq-temp-dataset "${BQ_TMP_DATASET}"
    --ses-ea-dir "${DX_SES_EA_PROXY_REGENIE_INPUT_DIR}"
    --work-dir "${LOCAL_OUT}"
    --workspace-scrap-dir "${WORKSPACE_SCRAP}"
)
if [[ "${REUSE_EXTRACTS}" -eq 1 ]]; then
    cmd+=(--reuse-extracts)
fi
if [[ "${FORCE}" -eq 1 ]]; then
    cmd+=(--force)
fi
if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
    cmd+=(--stage-aggregate)
fi

"${cmd[@]}"

if [[ "${MAKE_ETM_G}" -eq 1 ]]; then
    echo ""
    echo "=== Running downstream ETM general factor scoring ==="
    g_cmd=(bash "${SCRIPT_DIR}/run_etm_g_from_task_scores.sh")
    if [[ "${STAGE_AGGREGATE}" -eq 1 ]]; then
        g_cmd+=(--stage-aggregate)
    fi
    if [[ "${FORCE}" -eq 1 ]]; then
        g_cmd+=(--force)
    fi
    "${g_cmd[@]}"
fi

echo ""
echo "=== ETM cognitive task factor command complete ==="
