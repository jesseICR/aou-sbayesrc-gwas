#!/bin/bash
# get_genotypes.sh — AoU SBayesRC genotype-extraction pipeline (early phase).
#
# Steps:
#   1. Generate per-chromosome SBayesRC variant ID + idmap files (local)
#   2. Upload ID files to the workspace bucket on GCS
#   3. Extract matching variants from AoU acaf_threshold plink_bed via dsub,
#      remapping IDs to SBayesRC rsids in the same plink2 invocation
#
# Idempotent: each step checks for existing outputs and skips if already done.
# Submits Google Batch jobs via dsub; never downloads bulk genetics data; the
# workspace bucket and Google project are derived dynamically (no hardcoded
# user-specific values) so any AoU researcher can run this pipeline.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR

# ---------------------------------------------------------------------------
# Workspace / project detection (portable — do not hardcode user-specific IDs)
# ---------------------------------------------------------------------------
# Active workspace must be set beforehand: `wb workspace set --id=<your-id>`.
# Override any of these via env var if auto-detection fails.

if [[ -z "${GOOGLE_PROJECT:-}" ]]; then
    GOOGLE_PROJECT="$(wb workspace describe --format=json 2>/dev/null \
        | jq -r '.googleProjectId // .gcpProjectId // empty' 2>/dev/null || true)"
fi
if [[ -z "${GOOGLE_PROJECT:-}" ]]; then
    echo "ERROR: could not detect GOOGLE_PROJECT. Set it manually:"
    echo "  export GOOGLE_PROJECT=<your-workspace-google-project-id>"
    echo "Or run 'wb workspace set --id=<workspace-id>' first."
    exit 1
fi
export GOOGLE_PROJECT

if [[ -z "${WORKSPACE_BUCKET_URI:-}" ]]; then
    # `wb resource resolve --id=workspace-bucket` returns gs://<bucket-name>
    # for the workspace's primary writable bucket.
    WORKSPACE_BUCKET_URI="$(wb resource resolve --id=workspace-bucket 2>/dev/null || true)"
fi
if [[ -z "${WORKSPACE_BUCKET_URI:-}" ]]; then
    echo "ERROR: could not detect workspace bucket via 'wb resource resolve --id=workspace-bucket'."
    echo "Set it manually:"
    echo "  export WORKSPACE_BUCKET_URI=gs://<your-workspace-bucket-name>"
    exit 1
fi
export WORKSPACE_BUCKET_URI

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
# Output prefix on the workspace bucket — change to relocate pipeline output.
export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_URI%/}/sbayesrc_genotypes"

# AoU shared dataset bucket (requester-pays). Both vwb- and fc- mirrors work.
export AOU_DATASETS_BUCKET="${AOU_DATASETS_BUCKET:-gs://vwb-aou-datasets-controlled}"
export AOU_PLINK_BED_DIR="${AOU_DATASETS_BUCKET}/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed"

# Local paths
export LOCAL_SBAYESRC_ID_DIR="${SCRIPT_DIR}/data/sbayesrc_ids"
export ALIGNMENT_FILE="${SCRIPT_DIR}/data/support/sbayesrc_hg38.csv"

# GCS output paths (all under DX_OUTPUT_DIR)
export DX_SBAYESRC_ID_DIR="${DX_OUTPUT_DIR}/sbayesrc_ids"
export DX_WGS_PFILE_DIR="${DX_OUTPUT_DIR}/wgs_pfiles"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"

# dsub / Google Batch settings
export DSUB_REGION="${DSUB_REGION:-us-central1}"
export DSUB_MACHINE_TYPE="${DSUB_MACHINE_TYPE:-n2-standard-8}"
export PLINK2_IMAGE="${PLINK2_IMAGE:-quay.io/biocontainers/plink2:2.0.0a.6.16--h9948957_0}"

mkdir -p "${SCRIPT_DIR}/data/support" "${LOCAL_SBAYESRC_ID_DIR}" "${SCRIPT_DIR}/logs"

# ---------------------------------------------------------------------------
# Logging — tee terminal + file
# ---------------------------------------------------------------------------
LOG_FILE="${SCRIPT_DIR}/logs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "AoU SBayesRC Pipeline — $(date)"
echo "============================================"
echo "  GOOGLE_PROJECT       = ${GOOGLE_PROJECT}"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
echo "  DX_OUTPUT_DIR        = ${DX_OUTPUT_DIR}"
echo "  AOU_PLINK_BED_DIR  = ${AOU_PLINK_BED_DIR}"
echo "  DSUB_REGION        = ${DSUB_REGION}"
echo "  PLINK2_IMAGE       = ${PLINK2_IMAGE}"
echo "  LOG_FILE           = ${LOG_FILE}"
echo ""

# ---------------------------------------------------------------------------
# Setup: Python deps
# ---------------------------------------------------------------------------
echo "=== Setup: Python dependencies ==="
pip install -r "${SCRIPT_DIR}/requirements.txt" --quiet

# ---------------------------------------------------------------------------
# Setup: SBayesRC alignment file (hg38, ~7.35M SNPs)
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup: SBayesRC alignment file ==="
if [[ -s "${ALIGNMENT_FILE}" ]]; then
    echo "  Already cached at ${ALIGNMENT_FILE} — skipping download"
else
    echo "  Downloading sbayesrc_hg38.csv ..."
    curl -fsSL -o "${ALIGNMENT_FILE}" \
        "https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_hg38.csv"
    echo "  Downloaded ($(wc -l < "${ALIGNMENT_FILE}") lines)"
fi

# ---------------------------------------------------------------------------
# Step 1: Generate SBayesRC variant IDs + idmap (local)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 1: Generate SBayesRC variant IDs ==="
python3 "${SCRIPT_DIR}/store_sbayesrc_ids.py" \
    --input-file "${ALIGNMENT_FILE}" \
    --output-dir "${LOCAL_SBAYESRC_ID_DIR}"

# ---------------------------------------------------------------------------
# Step 2: Upload ID files to workspace bucket
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Upload SBayesRC IDs to GCS ==="
bash "${SCRIPT_DIR}/upload_sbayesrc_ids.sh"

# ---------------------------------------------------------------------------
# Step 3: Extract SBayesRC variants per chromosome (dsub / Google Batch)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 3: Extract SBayesRC variants per chromosome ==="
bash "${SCRIPT_DIR}/wgs_extract_variants.sh"

echo ""
echo "=== Pipeline complete ==="
