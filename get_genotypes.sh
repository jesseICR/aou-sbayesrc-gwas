#!/bin/bash
# get_genotypes.sh — AoU SBayesRC genotype-extraction pipeline.
#
# Runs interactively in an AoU Verily Jupyter terminal. Reads the locally
# FUSE-mounted controlled-tier dataset (read-only) and writes per-chromosome
# PLINK2 pfiles holding the ~7.35M SBayesRC SNPs to the locally FUSE-mounted
# workspace bucket (read-write), then builds the direct-SNP bfile used by
# REGENIE step 1.
#
# Steps:
#   1. Generate per-chromosome SBayesRC variant ID + idmap files (local)
#   2. Extract matching variants from AoU acaf_threshold pgen via plink2,
#      remapping IDs to SBayesRC rsids. Per-chrom work runs as parallel dsub
#      tasks on Google Batch in us-central1 (one Batch worker per chrom).
#   3. Prepare/extract the UKBB direct-SNP set from the extracted pfiles,
#      track absent direct SNPs for reporting, then merge present SNPs across
#      chr1..22 into direct_bfile/chr1_22_merged.{bed,bim,fam}.
#
# Idempotent: each step checks for existing outputs and skips if already done.
#
# Smoke-test override:
#   SBAYESRC_TEST_CHROM=21 bash get_genotypes.sh 2>&1
#     → runs only chr21 in step 2 and skips step 3, which requires all 22
#       extracted pfiles.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR

# ---------------------------------------------------------------------------
# Environment — populated by the AoU Verily Jupyter session
# ---------------------------------------------------------------------------
# $GOOGLE_PROJECT is preset by the Workbench. We do NOT trust $WORKSPACE_BUCKET
# — on at least some AoU pods it points to a non-existent cloned-bucket URI
# (e.g. gs://cloned-mybucket-...). We derive the real writable workspace
# bucket URI from the FUSE mount table instead (matches `wb resource resolve
# --id=workspace-bucket`).

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set — are you running inside an AoU Verily Jupyter session?}"

# ---------------------------------------------------------------------------
# Paths — FUSE mount points are stable across all AoU Verily users
# ---------------------------------------------------------------------------
# Controlled-tier dataset bucket (ro).
export AOU_DATA_MOUNT="/home/jupyter/workspace/data_controlled/vwb-aou-datasets-controlled"
# We use the pgen variant of the AoU acaf_threshold callset (not plink_bed):
# pgen is ~10× smaller than the equivalent bed (chr22 18 GB vs 168 GB), making
# gcsfuse-backed reads tractable. Pgen is multi-allelic with ID=".", so we
# split + assign IDs in pass 1 before extracting in pass 2 (see
# wgs_extract_variants.sh).
export AOU_PGEN_DIR="${AOU_DATA_MOUNT}/v8/wgs/short_read/snpindel/acaf_threshold/pgen"

# Workspace bucket (rw). All durable pipeline output lives under here.
export WORKSPACE_BUCKET_MOUNT="/home/jupyter/workspace/workspace-bucket"
WORKSPACE_BUCKET_URI="gs://$(mount | awk '/ on \/home\/jupyter\/workspace\/workspace-bucket /{print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    echo "  Is the workspace bucket FUSE-mounted? Check 'mount | grep workspace-bucket'." >&2
    exit 1
fi
export WORKSPACE_BUCKET_URI
# Mount-side paths (used for idempotency checks + small text writes).
export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/sbayesrc_genotypes"
export DX_SBAYESRC_ID_DIR="${DX_OUTPUT_DIR}/sbayesrc_ids"
export DX_WGS_PFILE_DIR="${DX_OUTPUT_DIR}/wgs_pfiles"
export DX_DIRECT_PFILE_DIR="${DX_OUTPUT_DIR}/direct_pfiles"
export DX_DIRECT_BFILE_DIR="${DX_OUTPUT_DIR}/direct_bfile"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"
# gs:// path for `gcloud storage cp` (large-pfile uploads — bypasses gcsfuse).
export DX_WGS_PFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/wgs_pfiles"
export DX_DIRECT_PFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/direct_pfiles"
export DX_DIRECT_BFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/direct_bfile"

# Local paths
export LOCAL_SBAYESRC_ID_DIR="${SCRIPT_DIR}/data/sbayesrc_ids"
export LOCAL_WGS_PFILE_DIR="${SCRIPT_DIR}/data/wgs_pfiles"
export LOCAL_DIRECT_SNPS_DIR="${SCRIPT_DIR}/data/support/direct_snps"
export LOCAL_DIRECT_SNPS_FILE="${LOCAL_DIRECT_SNPS_DIR}/ukbb_500k_qc_pass_direct_snps.txt"
export LOCAL_DIRECT_PREP_DIR="${SCRIPT_DIR}/data/direct_snps"
export ALIGNMENT_FILE="${SCRIPT_DIR}/data/support/sbayesrc_hg38.csv"
export DIRECT_SNPS_URL="https://raw.githubusercontent.com/jesseICR/ukbb-sbayesrc-gwas/main/data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt"

# Tools — plink2 is preinstalled on the AoU Verily Jupyter VM.
export PLINK2="${PLINK2:-/opt/workbench-tools/binaries/bin/plink2}"

# Threading — plink2 is multithreaded; default to all cores.
export THREADS="${THREADS:-$(nproc)}"

# ---------------------------------------------------------------------------
# dsub / Google Batch fan-out config
# ---------------------------------------------------------------------------
# Extraction runs as parallel dsub tasks against Google Batch (one task per
# chromosome). See CLAUDE.md "dsub from inside Jupyter does work" section for
# why each of these knobs is required.

export DSUB_PROVIDER="${DSUB_PROVIDER:-google-batch}"
export DSUB_REGION="${DSUB_REGION:-us-central1}"   # same region as workspace bucket
export DSUB_NETWORK="projects/${GOOGLE_PROJECT}/global/networks/network"
export DSUB_SUBNETWORK="projects/${GOOGLE_PROJECT}/regions/${DSUB_REGION}/subnetworks/subnetwork"
export DSUB_IMAGE="${DSUB_IMAGE:-marketplace.gcr.io/google/ubuntu2204}"
# Pet SA of this pod — the only identity dsub workers can actAs.
DSUB_PET_SA="$(gcloud config get-value account 2>/dev/null || true)"
if [[ -z "${DSUB_PET_SA}" ]]; then
    echo "ERROR: could not determine the pod's pet service account via 'gcloud config get-value account'." >&2
    exit 1
fi
export DSUB_PET_SA

# Per-task resource defaults (override via env). chr1 peak disk ~250 GB.
export DSUB_MIN_CORES="${DSUB_MIN_CORES:-4}"
export DSUB_MIN_RAM="${DSUB_MIN_RAM:-32}"
export DSUB_BOOT_DISK_SIZE="${DSUB_BOOT_DISK_SIZE:-50}"
export DSUB_DISK_SIZE="${DSUB_DISK_SIZE:-300}"

# Direct-bfile merge is a single large binary pmerge/write. Use an SSD data
# disk; pd-standard was much slower for this workload.
export DIRECT_BFILE_DSUB_MIN_CORES="${DIRECT_BFILE_DSUB_MIN_CORES:-8}"
export DIRECT_BFILE_DSUB_MIN_RAM="${DIRECT_BFILE_DSUB_MIN_RAM:-32}"
export DIRECT_BFILE_DSUB_DISK_SIZE="${DIRECT_BFILE_DSUB_DISK_SIZE:-300}"
export DIRECT_BFILE_DSUB_DISK_TYPE="${DIRECT_BFILE_DSUB_DISK_TYPE:-pd-ssd}"

# Worker-staging paths on the workspace bucket
export DSUB_BIN_URI="${WORKSPACE_BUCKET_URI}/bin"
export DSUB_PLINK2_GS="${DSUB_BIN_URI}/plink2"
export DSUB_SBAYESRC_ID_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/sbayesrc_ids"
export DSUB_LOG_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/logs/dsub"

# AoU controlled-tier pgen as gs:// (workers can't see the FUSE mount).
export AOU_PGEN_GS_DIR="gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen"

# ---------------------------------------------------------------------------
# Sanity checks (cheap, fail fast)
# ---------------------------------------------------------------------------
if [[ ! -d "${AOU_PGEN_DIR}" ]]; then
    echo "ERROR: ${AOU_PGEN_DIR} is not present."
    echo "  Is the controlled-tier dataset bucket mounted? Check 'mount | grep data_controlled'."
    exit 1
fi
if [[ ! -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
    echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is not writable."
    echo "  Is the workspace bucket FUSE-mounted? Check 'mount | grep workspace-bucket'."
    exit 1
fi
if [[ ! -x "${PLINK2}" ]]; then
    echo "ERROR: plink2 not found at ${PLINK2}."
    exit 1
fi

mkdir -p \
    "${SCRIPT_DIR}/data/support" \
    "${LOCAL_DIRECT_SNPS_DIR}" \
    "${LOCAL_DIRECT_PREP_DIR}" \
    "${LOCAL_SBAYESRC_ID_DIR}" \
    "${LOCAL_WGS_PFILE_DIR}" \
    "${SCRIPT_DIR}/logs/extract" \
    "${DX_OUTPUT_DIR}" \
    "${DX_SBAYESRC_ID_DIR}" \
    "${DX_WGS_PFILE_DIR}" \
    "${DX_DIRECT_PFILE_DIR}" \
    "${DX_DIRECT_BFILE_DIR}" \
    "${DX_LOGS_DIR}"

# ---------------------------------------------------------------------------
# Logging — tee terminal + per-run log file
# ---------------------------------------------------------------------------
LOG_FILE="${SCRIPT_DIR}/logs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "AoU SBayesRC Pipeline — $(date)"
echo "============================================"
echo "  GOOGLE_PROJECT       = ${GOOGLE_PROJECT}"
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}  (from mount table)"
echo "  WORKSPACE_BUCKET_MNT = ${WORKSPACE_BUCKET_MOUNT}"
echo "  AOU_PGEN_DIR         = ${AOU_PGEN_DIR}"
echo "  DX_OUTPUT_DIR        = ${DX_OUTPUT_DIR}"
echo "  LOCAL_WGS_PFILE_DIR  = ${LOCAL_WGS_PFILE_DIR}"
echo "  LOCAL_DIRECT_SNPS    = ${LOCAL_DIRECT_SNPS_FILE}"
echo "  PLINK2               = ${PLINK2} ($("${PLINK2}" --version 2>&1 | head -1))"
echo "  THREADS              = ${THREADS}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo "  DSUB_PROVIDER        = ${DSUB_PROVIDER}"
echo "  DSUB_REGION          = ${DSUB_REGION}"
echo "  DSUB_PET_SA          = ${DSUB_PET_SA}"
echo "  DSUB_IMAGE           = ${DSUB_IMAGE}"
echo "  DSUB worker resources= ${DSUB_MIN_CORES} vCPU, ${DSUB_MIN_RAM} GB RAM, ${DSUB_BOOT_DISK_SIZE}+${DSUB_DISK_SIZE} GB disk per task"
echo "  DIRECT_BFILE worker = ${DIRECT_BFILE_DSUB_MIN_CORES} vCPU, ${DIRECT_BFILE_DSUB_MIN_RAM} GB RAM, ${DIRECT_BFILE_DSUB_DISK_SIZE} GB ${DIRECT_BFILE_DSUB_DISK_TYPE}"
if [[ -n "${SBAYESRC_TEST_CHROM:-}" ]]; then
    echo "  SBAYESRC_TEST_CHROM  = ${SBAYESRC_TEST_CHROM}  (smoke-test mode)"
fi
echo ""

# ---------------------------------------------------------------------------
# Setup: Python deps
# ---------------------------------------------------------------------------
echo "=== Setup: Python dependencies ==="
pip install -r "${SCRIPT_DIR}/requirements.txt" --quiet

# ---------------------------------------------------------------------------
# Setup: SBayesRC alignment file (hg38, ~7.35M SNPs, ~200 MB)
# Hosted on a user GitHub release — public reference data, not AoU data.
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup: SBayesRC alignment file ==="
if [[ -s "${ALIGNMENT_FILE}" ]]; then
    echo "  Already cached at ${ALIGNMENT_FILE} ($(wc -l < "${ALIGNMENT_FILE}") lines) — skipping download"
else
    echo "  Downloading sbayesrc_hg38.csv ..."
    curl -fsSL -o "${ALIGNMENT_FILE}" \
        "https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_hg38.csv"
    echo "  Downloaded ($(wc -l < "${ALIGNMENT_FILE}") lines)"
fi

# ---------------------------------------------------------------------------
# Setup: direct SNP list for REGENIE step 1
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup: direct SNP list ==="
if [[ -s "${LOCAL_DIRECT_SNPS_FILE}" ]]; then
    echo "  Already cached at ${LOCAL_DIRECT_SNPS_FILE} ($(wc -l < "${LOCAL_DIRECT_SNPS_FILE}") lines) — skipping download"
else
    echo "  Downloading direct SNP list ..."
    curl -fsSL -o "${LOCAL_DIRECT_SNPS_FILE}" "${DIRECT_SNPS_URL}"
    echo "  Downloaded ($(wc -l < "${LOCAL_DIRECT_SNPS_FILE}") lines)"
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
# Step 2: Extract SBayesRC variants per chromosome (local plink2)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Extract SBayesRC variants per chromosome ==="
bash "${SCRIPT_DIR}/wgs_extract_variants.sh"

if [[ -n "${SBAYESRC_TEST_CHROM:-}" ]]; then
    echo ""
    echo "=== Step 3: Direct-SNP bfile ==="
    echo "  Skipping in SBAYESRC_TEST_CHROM mode; direct bfile requires all 22 chromosomes."
else
    # -----------------------------------------------------------------------
    # Step 3a: Prepare direct-SNP per-chromosome lists + missing metadata
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 3a: Prepare direct SNP metadata ==="
    python3 "${SCRIPT_DIR}/prepare_direct_snps.py" \
        --direct-snps "${LOCAL_DIRECT_SNPS_FILE}" \
        --alignment "${ALIGNMENT_FILE}" \
        --wgs-pfile-dir "${DX_WGS_PFILE_DIR}" \
        --output-dir "${LOCAL_DIRECT_PREP_DIR}"

    # -----------------------------------------------------------------------
    # Step 3b: Extract present direct SNPs
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 3b: Extract direct SNP pfiles ==="
    bash "${SCRIPT_DIR}/extract_direct_snps.sh"

    # -----------------------------------------------------------------------
    # Step 3c: Merge direct-SNP pfiles into one bfile
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 3c: Merge direct SNP bfile ==="
    bash "${SCRIPT_DIR}/make_direct_bfile.sh"
fi

echo ""
echo "=== Pipeline complete ==="
