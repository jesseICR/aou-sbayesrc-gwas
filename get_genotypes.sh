#!/bin/bash
# get_genotypes.sh — AoU SBayesRC genotype-extraction pipeline.
#
# Runs interactively in an AoU Verily Jupyter terminal. Reads the locally
# FUSE-mounted controlled-tier dataset (read-only) and writes per-chromosome
# PLINK2 pfiles holding the ~7.35M SBayesRC SNPs to the locally FUSE-mounted
# workspace bucket (read-write), then builds direct-SNP bfiles used by
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
#   4. Build a higher-quality direct bfile by filtering direct SNPs on AoU
#      EUR-only ALT-frequency concordance, EUR MAF, and EUR missingness.
#   5. Run ADMIXTURE K=6 projection from the HQ direct bfile, with an
#      ADMIXTURE-specific all-sample missingness filter and allele alignment.
#   6. Compare our ADMIXTURE fractions to AoU-provided ancestry fractions and
#      write the European keep-list used by downstream REGENIE steps.
#   7. Run KING kinship from HQ direct SNPs, compare to AoU's provided
#      relatedness table, and classify close relationships.
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
export DX_HQ_DIRECT_BFILE_DIR="${DX_OUTPUT_DIR}/direct_bfile_hq"
export DX_STATGEN_DIR="${DX_OUTPUT_DIR}/statgen"
export DX_ADMIXTURE_SCRAP_DIR="${DX_OUTPUT_DIR}/statgen/scrap"
export DX_ADMIXTURE_BATCH_DIR="${DX_OUTPUT_DIR}/statgen/scrap/batches"
export DX_ADMIXTURE_Q_DIR="${DX_OUTPUT_DIR}/statgen/scrap/q"
export DX_AOU_VS_OURS_DIR="${DX_OUTPUT_DIR}/statgen/aou_vs_ours"
export DX_EUROPEANS_DIR="${DX_OUTPUT_DIR}/europeans"
export DX_KINSHIP_DIR="${DX_OUTPUT_DIR}/kinship"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"
# gs:// path for `gcloud storage cp` (large-pfile uploads — bypasses gcsfuse).
export DX_WGS_PFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/wgs_pfiles"
export DX_DIRECT_PFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/direct_pfiles"
export DX_DIRECT_BFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/direct_bfile"
export DX_HQ_DIRECT_BFILE_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/direct_bfile_hq"
export DX_STATGEN_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/statgen"
export DX_ADMIXTURE_SCRAP_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/statgen/scrap"
export DX_ADMIXTURE_BATCH_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/statgen/scrap/batches"
export DX_ADMIXTURE_Q_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/statgen/scrap/q"
export DX_KINSHIP_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/kinship"

# Local paths
export LOCAL_SBAYESRC_ID_DIR="${SCRIPT_DIR}/data/sbayesrc_ids"
export LOCAL_WGS_PFILE_DIR="${SCRIPT_DIR}/data/wgs_pfiles"
export LOCAL_DIRECT_SNPS_DIR="${SCRIPT_DIR}/data/support/direct_snps"
export LOCAL_DIRECT_SNPS_FILE="${LOCAL_DIRECT_SNPS_DIR}/ukbb_500k_qc_pass_direct_snps.txt"
export LOCAL_DIRECT_PREP_DIR="${SCRIPT_DIR}/data/direct_snps"
export LOCAL_HQ_DIRECT_DIR="${SCRIPT_DIR}/data/high_quality_direct"
export LOCAL_ADMIXTURE_DIR="${SCRIPT_DIR}/data/admixture"
export LOCAL_ANCESTRY_COMPARE_DIR="${SCRIPT_DIR}/data/ancestry_compare"
export LOCAL_KINSHIP_DIR="${SCRIPT_DIR}/data/kinship"
export LOCAL_SNP_QC_FILE="${SCRIPT_DIR}/data/support/ukb_snp_qc.txt"
export ALIGNMENT_FILE="${SCRIPT_DIR}/data/support/sbayesrc_hg38.csv"
export SBAYESRC_LIFTOVER_FILE="${SCRIPT_DIR}/data/support/sbayesrc_liftover_results.csv"
export DIRECT_SNPS_URL="https://raw.githubusercontent.com/jesseICR/ukbb-sbayesrc-gwas/main/data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt"
export SBAYESRC_LIFTOVER_URL="https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_liftover_results.csv"
export ADMIXTURE_TSV_URL="https://raw.githubusercontent.com/jesseICR/public-statgen/main/outputs/admixture-global-6/admixture_allele_freqs.tsv"
export ADMIXTURE_DOWNLOAD_URL="https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz"
export UKB_SNP_QC_URL="https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt"

# AoU computed ancestry predictions; used to make the EUR keep-list for
# direct-SNP QC metrics. This file stays inside the AoU environment.
export AOU_ANCESTRY_PRED_FILE="${AOU_DATA_MOUNT}/v8/wgs/short_read/snpindel/aux/ancestry/echo_v4_r2.ancestry_preds.tsv"
export AOU_ADMIXTURE_Q_FILE="${AOU_DATA_MOUNT}/v8/wgs/short_read/snpindel/aux/admixture_estimates/aou_admixture_estimates_rye_v8.Q"
export AOU_RELATEDNESS_FILE="${AOU_DATA_MOUNT}/v8/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv"

# High-quality direct-bfile thresholds.
export HQ_AF_DIFF_MAX="${HQ_AF_DIFF_MAX:-0.04}"          # absolute ALT-frequency difference
export HQ_EUR_MAF_MIN="${HQ_EUR_MAF_MIN:-0.007}"        # AoU EUR MAF >= 0.7%
export HQ_EUR_MISSING_MAX="${HQ_EUR_MISSING_MAX:-0.05}" # AoU EUR variant missingness <= 5%

# ADMIXTURE K=6 projection settings.
export ADMIXTURE_K="${ADMIXTURE_K:-6}"
export ADMIXTURE_BATCH_SIZE="${ADMIXTURE_BATCH_SIZE:-20000}"
export ADMIXTURE_GENO_MAX="${ADMIXTURE_GENO_MAX:-0.05}" # all-sample variant missingness <= 5%

# European classifier and AoU-vs-ours ancestry comparison settings.
export OURS_EUR_MIN="${OURS_EUR_MIN:-0.8}"
export OURS_AFR_MAX="${OURS_AFR_MAX:-0.1}"
export OURS_AMR_MAX="${OURS_AMR_MAX:-0.1}"
export OURS_EAS_MAX="${OURS_EAS_MAX:-0.1}"
export OURS_OCE_MAX="${OURS_OCE_MAX:-0.1}"
export AOU_MID_THRESHOLDS="${AOU_MID_THRESHOLDS:-0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2}"

# KING kinship settings. The SNP missingness filter is all-sample, not EUR-only.
export KINSHIP_MISSING_MAX="${KINSHIP_MISSING_MAX:-0.01}"
export KING_TABLE_FILTER="${KING_TABLE_FILTER:-0.035}"
export KINSHIP_CLOSE_LOWER="${KINSHIP_CLOSE_LOWER:-0.1767}"
export KINSHIP_FIRST_DEGREE_UPPER="${KINSHIP_FIRST_DEGREE_UPPER:-0.3535}"
export KINSHIP_IBS0_CUTOFF="${KINSHIP_IBS0_CUTOFF:-0.0012}"
export KINSHIP_PROCEED_AFTER_SNP_REVIEW="${KINSHIP_PROCEED_AFTER_SNP_REVIEW:-0}"

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

# High-quality direct-bfile jobs scan/build a single large bfile.
export HQ_DIRECT_DSUB_MIN_CORES="${HQ_DIRECT_DSUB_MIN_CORES:-8}"
export HQ_DIRECT_DSUB_MIN_RAM="${HQ_DIRECT_DSUB_MIN_RAM:-32}"
export HQ_DIRECT_DSUB_DISK_SIZE="${HQ_DIRECT_DSUB_DISK_SIZE:-300}"
export HQ_DIRECT_DSUB_DISK_TYPE="${HQ_DIRECT_DSUB_DISK_TYPE:-pd-ssd}"

# ADMIXTURE prep/split scan large bfiles; projection runs small independent batches.
export ADMIXTURE_PREP_DSUB_MIN_CORES="${ADMIXTURE_PREP_DSUB_MIN_CORES:-8}"
export ADMIXTURE_PREP_DSUB_MIN_RAM="${ADMIXTURE_PREP_DSUB_MIN_RAM:-32}"
export ADMIXTURE_PREP_DSUB_DISK_SIZE="${ADMIXTURE_PREP_DSUB_DISK_SIZE:-300}"
export ADMIXTURE_PREP_DSUB_DISK_TYPE="${ADMIXTURE_PREP_DSUB_DISK_TYPE:-pd-ssd}"
export ADMIXTURE_SPLIT_DSUB_MIN_CORES="${ADMIXTURE_SPLIT_DSUB_MIN_CORES:-8}"
export ADMIXTURE_SPLIT_DSUB_MIN_RAM="${ADMIXTURE_SPLIT_DSUB_MIN_RAM:-32}"
export ADMIXTURE_SPLIT_DSUB_DISK_SIZE="${ADMIXTURE_SPLIT_DSUB_DISK_SIZE:-300}"
export ADMIXTURE_SPLIT_DSUB_DISK_TYPE="${ADMIXTURE_SPLIT_DSUB_DISK_TYPE:-pd-ssd}"
export ADMIXTURE_PROJECT_DSUB_MIN_CORES="${ADMIXTURE_PROJECT_DSUB_MIN_CORES:-8}"
export ADMIXTURE_PROJECT_DSUB_MIN_RAM="${ADMIXTURE_PROJECT_DSUB_MIN_RAM:-16}"
export ADMIXTURE_PROJECT_DSUB_DISK_SIZE="${ADMIXTURE_PROJECT_DSUB_DISK_SIZE:-100}"
export ADMIXTURE_PROJECT_DSUB_DISK_TYPE="${ADMIXTURE_PROJECT_DSUB_DISK_TYPE:-pd-ssd}"
export ADMIXTURE_CONCAT_DSUB_MIN_CORES="${ADMIXTURE_CONCAT_DSUB_MIN_CORES:-2}"
export ADMIXTURE_CONCAT_DSUB_MIN_RAM="${ADMIXTURE_CONCAT_DSUB_MIN_RAM:-8}"
export ADMIXTURE_CONCAT_DSUB_DISK_SIZE="${ADMIXTURE_CONCAT_DSUB_DISK_SIZE:-50}"
export ADMIXTURE_CONCAT_DSUB_DISK_TYPE="${ADMIXTURE_CONCAT_DSUB_DISK_TYPE:-pd-ssd}"

# Kinship subset/KING jobs scan the HQ direct bfile.
export KINSHIP_SUBSET_DSUB_MIN_CORES="${KINSHIP_SUBSET_DSUB_MIN_CORES:-8}"
export KINSHIP_SUBSET_DSUB_MIN_RAM="${KINSHIP_SUBSET_DSUB_MIN_RAM:-32}"
export KINSHIP_SUBSET_DSUB_DISK_SIZE="${KINSHIP_SUBSET_DSUB_DISK_SIZE:-300}"
export KINSHIP_SUBSET_DSUB_DISK_TYPE="${KINSHIP_SUBSET_DSUB_DISK_TYPE:-pd-ssd}"
export KING_DSUB_MIN_CORES="${KING_DSUB_MIN_CORES:-16}"
export KING_DSUB_MIN_RAM="${KING_DSUB_MIN_RAM:-64}"
export KING_DSUB_DISK_SIZE="${KING_DSUB_DISK_SIZE:-300}"
export KING_DSUB_DISK_TYPE="${KING_DSUB_DISK_TYPE:-pd-ssd}"

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
    "${LOCAL_HQ_DIRECT_DIR}" \
    "${LOCAL_ADMIXTURE_DIR}" \
    "${LOCAL_ANCESTRY_COMPARE_DIR}" \
    "${LOCAL_KINSHIP_DIR}" \
    "${LOCAL_SBAYESRC_ID_DIR}" \
    "${LOCAL_WGS_PFILE_DIR}" \
    "${SCRIPT_DIR}/logs/extract" \
    "${DX_OUTPUT_DIR}" \
    "${DX_SBAYESRC_ID_DIR}" \
    "${DX_WGS_PFILE_DIR}" \
    "${DX_DIRECT_PFILE_DIR}" \
    "${DX_DIRECT_BFILE_DIR}" \
    "${DX_HQ_DIRECT_BFILE_DIR}" \
    "${DX_STATGEN_DIR}" \
    "${DX_ADMIXTURE_SCRAP_DIR}" \
    "${DX_ADMIXTURE_BATCH_DIR}" \
    "${DX_ADMIXTURE_Q_DIR}" \
    "${DX_AOU_VS_OURS_DIR}" \
    "${DX_EUROPEANS_DIR}" \
    "${DX_KINSHIP_DIR}" \
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
echo "  LOCAL_HQ_DIRECT_DIR  = ${LOCAL_HQ_DIRECT_DIR}"
echo "  LOCAL_ADMIXTURE_DIR  = ${LOCAL_ADMIXTURE_DIR}"
echo "  LOCAL_ANCESTRY_CMP   = ${LOCAL_ANCESTRY_COMPARE_DIR}"
echo "  LOCAL_KINSHIP_DIR    = ${LOCAL_KINSHIP_DIR}"
echo "  PLINK2               = ${PLINK2} ($("${PLINK2}" --version 2>&1 | head -1))"
echo "  THREADS              = ${THREADS}"
echo "  LOG_FILE             = ${LOG_FILE}"
echo "  DSUB_PROVIDER        = ${DSUB_PROVIDER}"
echo "  DSUB_REGION          = ${DSUB_REGION}"
echo "  DSUB_PET_SA          = ${DSUB_PET_SA}"
echo "  DSUB_IMAGE           = ${DSUB_IMAGE}"
echo "  DSUB worker resources= ${DSUB_MIN_CORES} vCPU, ${DSUB_MIN_RAM} GB RAM, ${DSUB_BOOT_DISK_SIZE}+${DSUB_DISK_SIZE} GB disk per task"
echo "  DIRECT_BFILE worker = ${DIRECT_BFILE_DSUB_MIN_CORES} vCPU, ${DIRECT_BFILE_DSUB_MIN_RAM} GB RAM, ${DIRECT_BFILE_DSUB_DISK_SIZE} GB ${DIRECT_BFILE_DSUB_DISK_TYPE}"
echo "  HQ_DIRECT thresholds = |AoU_EUR_AF - SBayesRC_AF| <= ${HQ_AF_DIFF_MAX}, EUR MAF >= ${HQ_EUR_MAF_MIN}, EUR missingness <= ${HQ_EUR_MISSING_MAX}"
echo "  HQ_DIRECT worker     = ${HQ_DIRECT_DSUB_MIN_CORES} vCPU, ${HQ_DIRECT_DSUB_MIN_RAM} GB RAM, ${HQ_DIRECT_DSUB_DISK_SIZE} GB ${HQ_DIRECT_DSUB_DISK_TYPE}"
echo "  ADMIXTURE settings   = K=${ADMIXTURE_K}, batch_size=${ADMIXTURE_BATCH_SIZE}, source=direct_bfile_hq, all-sample geno <= ${ADMIXTURE_GENO_MAX}"
echo "  ADMIXTURE prep/split = ${ADMIXTURE_PREP_DSUB_MIN_CORES}/${ADMIXTURE_SPLIT_DSUB_MIN_CORES} vCPU, ${ADMIXTURE_PREP_DSUB_MIN_RAM}/${ADMIXTURE_SPLIT_DSUB_MIN_RAM} GB RAM, ${ADMIXTURE_PREP_DSUB_DISK_SIZE}/${ADMIXTURE_SPLIT_DSUB_DISK_SIZE} GB disk"
echo "  ADMIXTURE project    = ${ADMIXTURE_PROJECT_DSUB_MIN_CORES} vCPU, ${ADMIXTURE_PROJECT_DSUB_MIN_RAM} GB RAM, ${ADMIXTURE_PROJECT_DSUB_DISK_SIZE} GB ${ADMIXTURE_PROJECT_DSUB_DISK_TYPE}"
echo "  European classifier  = European >= ${OURS_EUR_MIN}, African/American/East_Asian/Oceanian <= ${OURS_AFR_MAX}/${OURS_AMR_MAX}/${OURS_EAS_MAX}/${OURS_OCE_MAX}"
echo "  Kinship settings     = source=direct_bfile_hq, UKB in_Relatedness SNPs, all-sample missingness < ${KINSHIP_MISSING_MAX}, KING filter >= ${KING_TABLE_FILTER}"
echo "  Kinship review gate  = KINSHIP_PROCEED_AFTER_SNP_REVIEW=${KINSHIP_PROCEED_AFTER_SNP_REVIEW}"
echo "  Kinship resources    = subset ${KINSHIP_SUBSET_DSUB_MIN_CORES} vCPU/${KINSHIP_SUBSET_DSUB_MIN_RAM} GB; KING ${KING_DSUB_MIN_CORES} vCPU/${KING_DSUB_MIN_RAM} GB"
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
# Setup: SBayesRC liftover/frequency file (public reference data)
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup: SBayesRC liftover frequency file ==="
if [[ -s "${SBAYESRC_LIFTOVER_FILE}" ]]; then
    echo "  Already cached at ${SBAYESRC_LIFTOVER_FILE} ($(wc -l < "${SBAYESRC_LIFTOVER_FILE}") lines) — skipping download"
else
    echo "  Downloading sbayesrc_liftover_results.csv ..."
    curl -fL --retry 3 --retry-delay 5 -o "${SBAYESRC_LIFTOVER_FILE}" "${SBAYESRC_LIFTOVER_URL}"
    echo "  Downloaded ($(wc -l < "${SBAYESRC_LIFTOVER_FILE}") lines)"
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

    # -----------------------------------------------------------------------
    # Step 4: Build high-quality direct-SNP bfile
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 4: Build high-quality direct SNP bfile ==="
    bash "${SCRIPT_DIR}/make_hq_direct_bfile.sh"

    # -----------------------------------------------------------------------
    # Step 5: ADMIXTURE K=6 projection
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 5a: Prepare ADMIXTURE projection inputs ==="
    bash "${SCRIPT_DIR}/admixture_prep.sh"

    echo ""
    echo "=== Step 5b: Split ADMIXTURE batches ==="
    bash "${SCRIPT_DIR}/admixture_split_batches.sh"

    echo ""
    echo "=== Step 5c: Run ADMIXTURE projection ==="
    bash "${SCRIPT_DIR}/admixture_run_projection.sh"

    # -----------------------------------------------------------------------
    # Step 6: AoU-vs-ours ancestry comparison + European classifier
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 6: Compare AoU ancestry fractions to ours ==="
    bash "${SCRIPT_DIR}/compare_aou_ancestry.sh"

    # -----------------------------------------------------------------------
    # Step 7: KING kinship + close relationship classification
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 7a: Build kinship SNP subset ==="
    bash "${SCRIPT_DIR}/subset_kinship_snps.sh"

    if [[ "${KINSHIP_PROCEED_AFTER_SNP_REVIEW}" != "1" ]]; then
        echo ""
        echo "=== Step 7: Pausing before KING kinship ==="
        if [[ -s "${DX_KINSHIP_DIR}/kinship_snp_subset_summary.tsv" ]]; then
            awk -F'\t' 'NR > 1 {printf "  %s = %s\n", $1, $2}' "${DX_KINSHIP_DIR}/kinship_snp_subset_summary.tsv"
        fi
        echo "  Review n_intersection_and_missing_lt_${KINSHIP_MISSING_MAX} before launching the large KING run."
        echo "  To continue with these settings, rerun with KINSHIP_PROCEED_AFTER_SNP_REVIEW=1."
        echo "  To override the SNP missingness threshold, rerun with e.g. KINSHIP_MISSING_MAX=0.02."
        echo ""
        echo "=== Pipeline paused for kinship SNP-count review ==="
        exit 0
    fi

    echo ""
    echo "=== Step 7b: Run KING kinship ==="
    bash "${SCRIPT_DIR}/run_king_kinship.sh"

    echo ""
    echo "=== Step 7c: QC KING kinship against AoU relatedness ==="
    bash "${SCRIPT_DIR}/kinship_qc.sh"

    echo ""
    echo "=== Step 7d: Classify close relationships ==="
    bash "${SCRIPT_DIR}/classify_relations.sh"
fi

echo ""
echo "=== Pipeline complete ==="
