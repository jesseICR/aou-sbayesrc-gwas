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
#   6. Classify European samples from our ADMIXTURE fractions and, when AoU RYE
#      fractions are available, compare our fractions to AoU-provided fractions.
#   7. Run KING kinship from HQ direct SNPs, compare to AoU's provided
#      relatedness table, and classify close relationships.
#   8. Select unrelated European IIDs for fitting PCA.
#   9. QC SNPs for PCA and build pca_ready.{bed,bim,fam}.
#  10. Fit PCA on unrelated Europeans and project PCs onto all samples.
#  11. Build sex covariate + sex-at-birth/WGS-ploidy concordance QC.
#  12. Build sample-QC exclusions for anomalous identical-genotype components.
#  13. Build final GWAS Step 1/Step 2 genotype inputs.
#  14. Build height GWAS phenotype/covariate inputs for all classified
#      European samples with valid program-collected height.
#  15. Optionally run a height GWAS with REGENIE (set RUN_HEIGHT_GWAS=1).
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
# Controlled-tier dataset bucket (ro). In current AoU Verily workspaces, the
# catalog resource is mounted directly under /home/jupyter/workspace with a
# versioned resource name, while the bucket still contains v7/v8/v9 subdirs.
export AOU_DATA_MOUNT="${AOU_DATA_MOUNT:-/home/jupyter/workspace/vwb-aou-datasets-controlled-v9}"

# AoU release selection. v9 is the default current release. The pipeline
# overrides WORKSPACE_CDR to the matching CDR because some Workbench sessions
# still export older CDR env vars even when newer CDRs are visible in the
# catalog and controlled-data mount.
export AOU_DATA_VERSION="${AOU_DATA_VERSION:-v9}"
case "${AOU_DATA_VERSION}" in
    v8|v9) ;;
    *)
        echo "ERROR: AOU_DATA_VERSION must be v8 or v9; got '${AOU_DATA_VERSION}'." >&2
        exit 1
        ;;
esac

infer_cdr_project() {
    if [[ -n "${WORKSPACE_CDR:-}" && "${WORKSPACE_CDR}" == *.* ]]; then
        printf '%s\n' "${WORKSPACE_CDR%%.*}"
    else
        printf '%s\n' "wb-silky-artichoke-2408"
    fi
}

export AOU_CDR_PROJECT="${AOU_CDR_PROJECT:-$(infer_cdr_project)}"
case "${AOU_DATA_VERSION}" in
    v9)
        export AOU_CDR_DATASET="${AOU_CDR_DATASET:-C2025Q4R6}"
        export AOU_WGS_RELEASE_ROOT="${AOU_DATA_MOUNT}/v9/wgs/short_read/snpindel"
        export AOU_WGS_GS_ROOT="gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel"
        export AOU_PGEN_DIR="${AOU_WGS_RELEASE_ROOT}/acaf_threshold/pgen"
        export AOU_PGEN_GS_DIR="${AOU_WGS_GS_ROOT}/acaf_threshold/pgen"
        export AOU_ANCESTRY_PRED_FILE="${AOU_WGS_RELEASE_ROOT}/aux/ancestry/ancestry_preds.tsv"
        export AOU_ADMIXTURE_Q_FILE="${AOU_ADMIXTURE_Q_FILE:-}"
        export AOU_RYE_COMPARISON_MODE="${AOU_RYE_COMPARISON_MODE:-skip}"
        export AOU_RELATEDNESS_FILE="${AOU_WGS_RELEASE_ROOT}/aux/relatedness/samples_relatedness.tsv"
        export AOU_GENOMIC_METRICS_FILE="${AOU_WGS_RELEASE_ROOT}/aux/qc/genomics_metrics_May042026_1724_12_tz0000.tsv"
        export SBAYESRC_OUTPUT_PREFIX="${SBAYESRC_OUTPUT_PREFIX:-sbayesrc_genotypes}"
        ;;
    v8)
        if [[ -z "${AOU_CDR_DATASET:-}" && -n "${WORKSPACE_CDR:-}" && "${WORKSPACE_CDR}" == *.* ]]; then
            export AOU_CDR_DATASET="${WORKSPACE_CDR##*.}"
        fi
        if [[ -z "${AOU_CDR_DATASET:-}" ]]; then
            export AOU_CDR_DATASET="C2024Q3R9"
        fi
        export AOU_WGS_RELEASE_ROOT="${AOU_DATA_MOUNT}/v8/wgs/short_read/snpindel"
        export AOU_WGS_GS_ROOT="gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel"
        export AOU_PGEN_DIR="${AOU_WGS_RELEASE_ROOT}/acaf_threshold/pgen"
        export AOU_PGEN_GS_DIR="${AOU_WGS_GS_ROOT}/acaf_threshold/pgen"
        export AOU_ANCESTRY_PRED_FILE="${AOU_WGS_RELEASE_ROOT}/aux/ancestry/echo_v4_r2.ancestry_preds.tsv"
        export AOU_ADMIXTURE_Q_FILE="${AOU_ADMIXTURE_Q_FILE:-${AOU_WGS_RELEASE_ROOT}/aux/admixture_estimates/aou_admixture_estimates_rye_v8.Q}"
        export AOU_RYE_COMPARISON_MODE="${AOU_RYE_COMPARISON_MODE:-auto}"
        export AOU_RELATEDNESS_FILE="${AOU_WGS_RELEASE_ROOT}/aux/relatedness/samples_relatedness.tsv"
        export AOU_GENOMIC_METRICS_FILE="${AOU_WGS_RELEASE_ROOT}/aux/qc/genomics_metrics_Dec142023_1859_02_tz0000.tsv"
        export SBAYESRC_OUTPUT_PREFIX="${SBAYESRC_OUTPUT_PREFIX:-sbayesrc_genotypes}"
        ;;
esac

export AOU_TARGET_WORKSPACE_CDR="${AOU_CDR_PROJECT}.${AOU_CDR_DATASET}"
if [[ "${AOU_STRICT_WORKSPACE_CDR:-0}" == "1" && -n "${WORKSPACE_CDR:-}" && "${WORKSPACE_CDR}" != "${AOU_TARGET_WORKSPACE_CDR}" ]]; then
    echo "ERROR: WORKSPACE_CDR=${WORKSPACE_CDR}, but AOU_DATA_VERSION=${AOU_DATA_VERSION} expects ${AOU_TARGET_WORKSPACE_CDR}." >&2
    echo "  Unset AOU_STRICT_WORKSPACE_CDR or set AOU_CDR_DATASET/AOU_DATA_VERSION consistently." >&2
    exit 1
fi
export WORKSPACE_CDR_ORIGINAL="${WORKSPACE_CDR:-}"
export WORKSPACE_CDR="${AOU_TARGET_WORKSPACE_CDR}"

# Workspace bucket (rw). All durable pipeline output lives under here.
export WORKSPACE_BUCKET_MOUNT="/home/jupyter/workspace/workspace-bucket"

is_workspace_bucket_fuse_mounted() {
    mount | awk -v target="${WORKSPACE_BUCKET_MOUNT}" '
        $2 == "on" && $3 == target && $4 == "type" && $5 == "fuse.gcsfuse" {
            found = 1
        }
        END { exit found ? 0 : 1 }
    '
}

run_wb() {
    local timeout_seconds="${WB_RESOURCE_TIMEOUT_SECONDS:-300}"
    if command -v timeout >/dev/null 2>&1; then
        timeout "${timeout_seconds}" wb "$@"
    else
        wb "$@"
    fi
}

ensure_workspace_bucket_mount() {
    if is_workspace_bucket_fuse_mounted && [[ -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
        return 0
    fi

    echo "Workspace bucket is not currently mounted at ${WORKSPACE_BUCKET_MOUNT}; attempting Workbench setup ..."
    if ! command -v wb >/dev/null 2>&1; then
        echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is not a writable gcsfuse mount, and 'wb' is unavailable." >&2
        echo "  This pipeline needs the AoU/Verily workspace bucket mounted as durable GCS-backed storage." >&2
        exit 1
    fi

    echo "  Mounting existing Workbench bucket resources ..."
    if run_wb resource mount --allow-other &&
       is_workspace_bucket_fuse_mounted && [[ -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
        return 0
    fi

    if run_wb resource resolve --id=workspace-bucket >/dev/null 2>&1; then
        echo "  Found Workbench resource workspace-bucket."
    else
        echo "  Workbench resource workspace-bucket not found; attempting to create it ..."
        if ! run_wb resource create gcs-bucket --id=workspace-bucket --location=US-CENTRAL1; then
            if run_wb resource resolve --id=workspace-bucket >/dev/null 2>&1; then
                echo "  Workbench resource workspace-bucket is now available after create returned nonzero; continuing."
            else
                echo "ERROR: could not create Workbench GCS resource workspace-bucket." >&2
                echo "  Create or attach a writable workspace bucket resource in AoU Workbench, then restart/remount the session." >&2
                exit 1
            fi
        fi
    fi

    echo "  Mounting Workbench bucket resources ..."
    if ! run_wb resource mount --allow-other; then
        echo "ERROR: 'wb resource mount' failed." >&2
        echo "  The pipeline cannot safely continue without a real gcsfuse workspace bucket mount." >&2
        exit 1
    fi

    if ! is_workspace_bucket_fuse_mounted || [[ ! -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
        echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is still not a writable gcsfuse mount after Workbench setup." >&2
        echo "  Do not create this path manually with mkdir; that would write outputs to local ephemeral disk." >&2
        echo "  Check 'mount | grep workspace-bucket' or recreate/restart the AoU Jupyter environment." >&2
        exit 1
    fi
}

ensure_workspace_bucket_mount
WORKSPACE_BUCKET_URI="gs://$(mount | awk '/ on \/home\/jupyter\/workspace\/workspace-bucket /{print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    echo "  Is the workspace bucket FUSE-mounted? Check 'mount | grep workspace-bucket'." >&2
    exit 1
fi
export WORKSPACE_BUCKET_URI
# Mount-side paths (used for idempotency checks + small text writes).
export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/${SBAYESRC_OUTPUT_PREFIX}"
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
export DX_PCA_EUR_DIR="${DX_OUTPUT_DIR}/pca_eur"
export DX_GENETIC_SEX_DIR="${DX_OUTPUT_DIR}/genetic_sex"
export DX_SAMPLE_QC_DIR="${DX_OUTPUT_DIR}/sample_qc"
export DX_GWAS_GENOTYPES_DIR="${DX_OUTPUT_DIR}/gwas_genotypes"
export DX_GWAS_STEP1_BFILE_DIR="${DX_GWAS_GENOTYPES_DIR}/step1_direct"
export DX_GWAS_STEP2_PFILE_DIR="${DX_GWAS_GENOTYPES_DIR}/step2_wgs_pfiles"
export DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
export DX_HEIGHT_REGENIE_INPUT_DIR="${DX_REGENIE_INPUT_DIR}/height_example"
export DX_REGENIE_OUTPUT_DIR="${DX_OUTPUT_DIR}/regenie_output"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"
# gs:// path for `gcloud storage cp` (large-pfile uploads — bypasses gcsfuse).
export DX_OUTPUT_URI="${WORKSPACE_BUCKET_URI}/${SBAYESRC_OUTPUT_PREFIX}"
export DX_WGS_PFILE_URI="${DX_OUTPUT_URI}/wgs_pfiles"
export DX_DIRECT_PFILE_URI="${DX_OUTPUT_URI}/direct_pfiles"
export DX_DIRECT_BFILE_URI="${DX_OUTPUT_URI}/direct_bfile"
export DX_HQ_DIRECT_BFILE_URI="${DX_OUTPUT_URI}/direct_bfile_hq"
export DX_STATGEN_URI="${DX_OUTPUT_URI}/statgen"
export DX_ADMIXTURE_SCRAP_URI="${DX_OUTPUT_URI}/statgen/scrap"
export DX_ADMIXTURE_BATCH_URI="${DX_OUTPUT_URI}/statgen/scrap/batches"
export DX_ADMIXTURE_Q_URI="${DX_OUTPUT_URI}/statgen/scrap/q"
export DX_KINSHIP_URI="${DX_OUTPUT_URI}/kinship"
export DX_PCA_EUR_URI="${DX_OUTPUT_URI}/pca_eur"
export DX_GENETIC_SEX_URI="${DX_OUTPUT_URI}/genetic_sex"
export DX_SAMPLE_QC_URI="${DX_OUTPUT_URI}/sample_qc"
export DX_EUROPEANS_URI="${DX_OUTPUT_URI}/europeans"
export DX_GWAS_GENOTYPES_URI="${DX_OUTPUT_URI}/gwas_genotypes"
export DX_GWAS_STEP1_BFILE_URI="${DX_GWAS_GENOTYPES_URI}/step1_direct"
export DX_GWAS_STEP2_PFILE_URI="${DX_GWAS_GENOTYPES_URI}/step2_wgs_pfiles"
export DX_REGENIE_INPUT_URI="${DX_OUTPUT_URI}/regenie_input"
export DX_HEIGHT_REGENIE_INPUT_URI="${DX_REGENIE_INPUT_URI}/height_example"
export DX_REGENIE_OUTPUT_URI="${DX_OUTPUT_URI}/regenie_output"

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
export LOCAL_PCA_QC_DIR="${SCRIPT_DIR}/data/pca_qc"
export LOCAL_GWAS_GENOTYPES_DIR="${SCRIPT_DIR}/data/gwas_genotypes"
export LOCAL_REGENIE_DIR="${SCRIPT_DIR}/data/regenie"
export LOCAL_SNP_QC_FILE="${SCRIPT_DIR}/data/support/ukb_snp_qc.txt"
export ALIGNMENT_FILE="${SCRIPT_DIR}/data/support/sbayesrc_hg38.csv"
export SBAYESRC_LIFTOVER_FILE="${SCRIPT_DIR}/data/support/sbayesrc_liftover_results.csv"
export DIRECT_SNPS_URL="https://raw.githubusercontent.com/jesseICR/ukbb-sbayesrc-gwas/main/data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt"
export SBAYESRC_LIFTOVER_URL="https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_liftover_results.csv"
export ADMIXTURE_TSV_URL="https://raw.githubusercontent.com/jesseICR/public-statgen/main/outputs/admixture-global-6/admixture_allele_freqs.tsv"
export ADMIXTURE_DOWNLOAD_URL="https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz"
export UKB_SNP_QC_URL="https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt"
export PCA_HIGH_LD_URL="https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed"

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

# PCA European training sample selection.
export PCA_KINSHIP_THRESHOLD="${PCA_KINSHIP_THRESHOLD:-0.0441941}" # 0.5^(9/2), third-degree lower bound
export PCA_SEED_RELATIONSHIPS="${PCA_SEED_RELATIONSHIPS:-sibling,identical}"

# PCA SNP QC settings.
export PCA_AF_DIFF_MAX="${PCA_AF_DIFF_MAX:-0.03}"
export PCA_MAF_MIN="${PCA_MAF_MIN:-0.01}"
export PCA_GENO_MAX="${PCA_GENO_MAX:-0.01}"
export PCA_MIND_MAX="${PCA_MIND_MAX:-0.01}"
export PCA_LD_WINDOW="${PCA_LD_WINDOW:-1000}"
export PCA_LD_STEP="${PCA_LD_STEP:-80}"
export PCA_LD_R2="${PCA_LD_R2:-0.1}"

# PCA fitting/projection settings.
export PCA_NPCS="${PCA_NPCS:-20}"
export PCA_SEED="${PCA_SEED:-0}"

# Sex covariate/QC settings. sex_covar.txt keeps binary sex-at-birth samples
# with WGS sex ploidy concordance by default.
export GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE="${GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE:-1}"
export SBAYESRC_BQ_TMP_DATASET="${SBAYESRC_BQ_TMP_DATASET:-}"

# Sample QC exclusions. Components of 3+ samples with genetically identical
# profiles are not plausible ordinary twin pairs, so they are excluded from
# downstream GWAS sample sets.
export IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE="${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE:-3}"

# Final GWAS genotype input filters.
# Step 1 direct bfile: geno in our classified Europeans, AF/MAF in fit_pca_iids.
export GWAS_STEP1_GENO_MAX="${GWAS_STEP1_GENO_MAX:-0.01}"
export GWAS_STEP1_AF_DIFF_MAX="${GWAS_STEP1_AF_DIFF_MAX:-0.03}"
export GWAS_STEP1_MAF_MIN="${GWAS_STEP1_MAF_MIN:-0.007}"
# Step 2 WGS pfiles: geno in our classified Europeans, AF/MAF in fit_pca_iids.
export GWAS_STEP2_GENO_MAX="${GWAS_STEP2_GENO_MAX:-0.03}"
export GWAS_STEP2_AF_DIFF_MAX="${GWAS_STEP2_AF_DIFF_MAX:-0.04}"
export GWAS_STEP2_MAF_MIN="${GWAS_STEP2_MAF_MIN:-0.007}"

# Height GWAS example settings. Defaults use program-collected AoU height:
# measurement 3036277 (Body height), source 903133 (Height), type 44818701
# (From physical examination), unit 8582 (centimeter).
export HEIGHT_MEASUREMENT_CONCEPT_ID="${HEIGHT_MEASUREMENT_CONCEPT_ID:-3036277}"
export HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID="${HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID:-903133}"
export HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID="${HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID:-44818701}"
export HEIGHT_UNIT_CONCEPT_ID="${HEIGHT_UNIT_CONCEPT_ID:-8582}"
export HEIGHT_MIN_CM="${HEIGHT_MIN_CM:-140}"
export HEIGHT_N_PCS="${HEIGHT_N_PCS:-10}"

# REGENIE height GWAS settings. The full GWAS is intentionally gated because it
# launches one Step 1 job plus per-chromosome Step 2 jobs.
export RUN_HEIGHT_GWAS="${RUN_HEIGHT_GWAS:-0}"
export HEIGHT_GWAS_INPUT_NAME="${HEIGHT_GWAS_INPUT_NAME:-height_example}"
export HEIGHT_GWAS_OUTPUT_NAME="${HEIGHT_GWAS_OUTPUT_NAME:-height_example}"
export REGENIE_APPLY_RINT="${REGENIE_APPLY_RINT:-1}"
export REGENIE_STEP1_BLOCK_SIZE="${REGENIE_STEP1_BLOCK_SIZE:-1000}"
export REGENIE_STEP2_BLOCK_SIZE="${REGENIE_STEP2_BLOCK_SIZE:-200}"
export REGENIE_CHROMS="${REGENIE_CHROMS:-1-22}"

# Tools — plink2 is preinstalled on the AoU Verily Jupyter VM.
export PLINK2="${PLINK2:-/opt/workbench-tools/binaries/bin/plink2}"
if [[ -z "${REGENIE:-}" ]]; then
    if command -v regenie >/dev/null 2>&1; then
        export REGENIE="$(command -v regenie)"
    elif [[ -x "/opt/workbench-tools/binaries/bin/regenie" ]]; then
        export REGENIE="/opt/workbench-tools/binaries/bin/regenie"
    else
        export REGENIE=""
    fi
fi

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
export KING_DSUB_MIN_CORES="${KING_DSUB_MIN_CORES:-32}"
export KING_DSUB_MIN_RAM="${KING_DSUB_MIN_RAM:-256}"
export KING_DSUB_DISK_SIZE="${KING_DSUB_DISK_SIZE:-300}"
export KING_DSUB_DISK_TYPE="${KING_DSUB_DISK_TYPE:-pd-ssd}"

# PCA SNP QC scans/builds bfiles from the HQ direct bfile.
export PCA_SNP_QC_DSUB_MIN_CORES="${PCA_SNP_QC_DSUB_MIN_CORES:-8}"
export PCA_SNP_QC_DSUB_MIN_RAM="${PCA_SNP_QC_DSUB_MIN_RAM:-32}"
export PCA_SNP_QC_DSUB_DISK_SIZE="${PCA_SNP_QC_DSUB_DISK_SIZE:-300}"
export PCA_SNP_QC_DSUB_DISK_TYPE="${PCA_SNP_QC_DSUB_DISK_TYPE:-pd-ssd}"

# PCA fitting/projection localizes pca_ready and the all-sample HQ direct bfile.
export PCA_PROJECT_DSUB_MIN_CORES="${PCA_PROJECT_DSUB_MIN_CORES:-16}"
export PCA_PROJECT_DSUB_MIN_RAM="${PCA_PROJECT_DSUB_MIN_RAM:-64}"
export PCA_PROJECT_DSUB_DISK_SIZE="${PCA_PROJECT_DSUB_DISK_SIZE:-300}"
export PCA_PROJECT_DSUB_DISK_TYPE="${PCA_PROJECT_DSUB_DISK_TYPE:-pd-ssd}"

# REGENIE Step 1 localizes the final GWAS direct bfile. Step 2 localizes one
# final GWAS WGS pfile per chromosome plus the Step 1 LOCO predictions.
export REGENIE_STEP1_DSUB_MIN_CORES="${REGENIE_STEP1_DSUB_MIN_CORES:-16}"
export REGENIE_STEP1_DSUB_MIN_RAM="${REGENIE_STEP1_DSUB_MIN_RAM:-64}"
export REGENIE_STEP1_DSUB_DISK_SIZE="${REGENIE_STEP1_DSUB_DISK_SIZE:-300}"
export REGENIE_STEP1_DSUB_DISK_TYPE="${REGENIE_STEP1_DSUB_DISK_TYPE:-pd-ssd}"
export REGENIE_STEP2_DSUB_MIN_CORES="${REGENIE_STEP2_DSUB_MIN_CORES:-8}"
export REGENIE_STEP2_DSUB_MIN_RAM="${REGENIE_STEP2_DSUB_MIN_RAM:-32}"
export REGENIE_STEP2_DSUB_DISK_SIZE="${REGENIE_STEP2_DSUB_DISK_SIZE:-300}"
export REGENIE_STEP2_DSUB_DISK_TYPE="${REGENIE_STEP2_DSUB_DISK_TYPE:-pd-ssd}"

# Final GWAS genotype input jobs.
export GWAS_METRICS_DSUB_MIN_CORES="${GWAS_METRICS_DSUB_MIN_CORES:-4}"
export GWAS_METRICS_DSUB_MIN_RAM="${GWAS_METRICS_DSUB_MIN_RAM:-24}"
export GWAS_METRICS_DSUB_DISK_SIZE="${GWAS_METRICS_DSUB_DISK_SIZE:-150}"
export GWAS_METRICS_DSUB_DISK_TYPE="${GWAS_METRICS_DSUB_DISK_TYPE:-pd-ssd}"
export GWAS_DIRECT_DSUB_MIN_CORES="${GWAS_DIRECT_DSUB_MIN_CORES:-8}"
export GWAS_DIRECT_DSUB_MIN_RAM="${GWAS_DIRECT_DSUB_MIN_RAM:-32}"
export GWAS_DIRECT_DSUB_DISK_SIZE="${GWAS_DIRECT_DSUB_DISK_SIZE:-200}"
export GWAS_DIRECT_DSUB_DISK_TYPE="${GWAS_DIRECT_DSUB_DISK_TYPE:-pd-ssd}"
export GWAS_WGS_DSUB_MIN_CORES="${GWAS_WGS_DSUB_MIN_CORES:-4}"
export GWAS_WGS_DSUB_MIN_RAM="${GWAS_WGS_DSUB_MIN_RAM:-24}"
export GWAS_WGS_DSUB_DISK_SIZE="${GWAS_WGS_DSUB_DISK_SIZE:-180}"
export GWAS_WGS_DSUB_DISK_TYPE="${GWAS_WGS_DSUB_DISK_TYPE:-pd-ssd}"

# Worker-staging paths on the workspace bucket
export DSUB_BIN_URI="${WORKSPACE_BUCKET_URI}/bin"
export DSUB_PLINK2_GS="${DSUB_BIN_URI}/plink2"
export DSUB_REGENIE_BUNDLE_URI="${DSUB_BIN_URI}/regenie_bundle"
export DSUB_SBAYESRC_ID_URI="${DX_OUTPUT_URI}/sbayesrc_ids"
export DSUB_LOG_URI="${DX_OUTPUT_URI}/logs/dsub"

# ---------------------------------------------------------------------------
# Sanity checks (cheap, fail fast)
# ---------------------------------------------------------------------------
if [[ ! -d "${AOU_PGEN_DIR}" ]]; then
    echo "ERROR: ${AOU_PGEN_DIR} is not present."
    echo "  Is the controlled-tier dataset bucket mounted? Check 'mount | grep vwb-aou-datasets-controlled'."
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
if [[ -z "${WORKSPACE_CDR:-}" ]]; then
    echo "ERROR: WORKSPACE_CDR is not set — are you running inside an AoU Verily Jupyter session?" >&2
    exit 1
fi
for required_aou_file in "${AOU_ANCESTRY_PRED_FILE}" "${AOU_RELATEDNESS_FILE}" "${AOU_GENOMIC_METRICS_FILE}"; do
    if [[ ! -s "${required_aou_file}" ]]; then
        echo "ERROR: missing AoU ${AOU_DATA_VERSION} auxiliary file ${required_aou_file}" >&2
        exit 1
    fi
done
if [[ "${AOU_RYE_COMPARISON_MODE}" == "required" && ! -s "${AOU_ADMIXTURE_Q_FILE:-}" ]]; then
    echo "ERROR: AOU_RYE_COMPARISON_MODE=required but AOU_ADMIXTURE_Q_FILE is missing or empty." >&2
    echo "  AOU_ADMIXTURE_Q_FILE=${AOU_ADMIXTURE_Q_FILE:-unset}" >&2
    exit 1
fi

regenie_version_text() {
    local help_text
    help_text="$("$1" --help 2>&1 || true)"
    printf '%s\n' "${help_text}" | sed -n '2{s/^ *//;p;q;}'
}

mkdir -p \
    "${SCRIPT_DIR}/data/support" \
    "${LOCAL_DIRECT_SNPS_DIR}" \
    "${LOCAL_DIRECT_PREP_DIR}" \
    "${LOCAL_HQ_DIRECT_DIR}" \
    "${LOCAL_ADMIXTURE_DIR}" \
    "${LOCAL_ANCESTRY_COMPARE_DIR}" \
    "${LOCAL_KINSHIP_DIR}" \
    "${LOCAL_PCA_QC_DIR}" \
    "${LOCAL_GWAS_GENOTYPES_DIR}" \
    "${LOCAL_REGENIE_DIR}" \
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
    "${DX_PCA_EUR_DIR}" \
    "${DX_GENETIC_SEX_DIR}" \
    "${DX_SAMPLE_QC_DIR}" \
    "${DX_GWAS_GENOTYPES_DIR}" \
    "${DX_GWAS_STEP1_BFILE_DIR}" \
    "${DX_GWAS_STEP2_PFILE_DIR}" \
    "${DX_REGENIE_INPUT_DIR}" \
    "${DX_HEIGHT_REGENIE_INPUT_DIR}" \
    "${DX_REGENIE_OUTPUT_DIR}" \
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
echo "  AOU_DATA_VERSION     = ${AOU_DATA_VERSION}"
echo "  AOU_TARGET_CDR       = ${AOU_TARGET_WORKSPACE_CDR}"
if [[ -n "${WORKSPACE_CDR_ORIGINAL}" && "${WORKSPACE_CDR_ORIGINAL}" != "${WORKSPACE_CDR}" ]]; then
    echo "  WORKSPACE_CDR_ORIG   = ${WORKSPACE_CDR_ORIGINAL}  (overridden for ${AOU_DATA_VERSION})"
fi
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}  (from mount table)"
echo "  WORKSPACE_BUCKET_MNT = ${WORKSPACE_BUCKET_MOUNT}"
echo "  OUTPUT_PREFIX        = ${SBAYESRC_OUTPUT_PREFIX}"
echo "  AOU_PGEN_DIR         = ${AOU_PGEN_DIR}"
echo "  AOU_PGEN_GS_DIR      = ${AOU_PGEN_GS_DIR}"
echo "  DX_OUTPUT_DIR        = ${DX_OUTPUT_DIR}"
echo "  WORKSPACE_CDR        = ${WORKSPACE_CDR}"
echo "  LOCAL_WGS_PFILE_DIR  = ${LOCAL_WGS_PFILE_DIR}"
echo "  LOCAL_DIRECT_SNPS    = ${LOCAL_DIRECT_SNPS_FILE}"
echo "  LOCAL_HQ_DIRECT_DIR  = ${LOCAL_HQ_DIRECT_DIR}"
echo "  LOCAL_ADMIXTURE_DIR  = ${LOCAL_ADMIXTURE_DIR}"
echo "  LOCAL_ANCESTRY_CMP   = ${LOCAL_ANCESTRY_COMPARE_DIR}"
echo "  LOCAL_KINSHIP_DIR    = ${LOCAL_KINSHIP_DIR}"
echo "  LOCAL_REGENIE_DIR    = ${LOCAL_REGENIE_DIR}"
echo "  PLINK2               = ${PLINK2} ($("${PLINK2}" --version 2>&1 | head -1))"
if [[ -n "${REGENIE}" && -x "${REGENIE}" ]]; then
    echo "  REGENIE              = ${REGENIE} ($(regenie_version_text "${REGENIE}"))"
else
    echo "  REGENIE              = not found (required only when RUN_HEIGHT_GWAS=1)"
fi
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
echo "  AoU RYE comparison   = ${AOU_RYE_COMPARISON_MODE}; q=${AOU_ADMIXTURE_Q_FILE:-none}"
echo "  Kinship settings     = source=direct_bfile_hq, UKB in_Relatedness SNPs, all-sample missingness < ${KINSHIP_MISSING_MAX}, KING filter >= ${KING_TABLE_FILTER}"
echo "  Kinship resources    = subset ${KINSHIP_SUBSET_DSUB_MIN_CORES} vCPU/${KINSHIP_SUBSET_DSUB_MIN_RAM} GB; KING ${KING_DSUB_MIN_CORES} vCPU/${KING_DSUB_MIN_RAM} GB"
echo "  PCA EUR selection    = seed relationships ${PCA_SEED_RELATIONSHIPS}, kinship threshold ${PCA_KINSHIP_THRESHOLD}"
echo "  PCA SNP QC           = source=direct_bfile_hq, AF diff <= ${PCA_AF_DIFF_MAX}, MAF >= ${PCA_MAF_MIN}, geno <= ${PCA_GENO_MAX}, mind <= ${PCA_MIND_MAX}, LD ${PCA_LD_WINDOW}/${PCA_LD_STEP}/${PCA_LD_R2}"
echo "  PCA fit/project      = ${PCA_NPCS} PCs, seed ${PCA_SEED}, source=direct_bfile_hq, resources ${PCA_PROJECT_DSUB_MIN_CORES} vCPU/${PCA_PROJECT_DSUB_MIN_RAM} GB/${PCA_PROJECT_DSUB_DISK_SIZE} GB ${PCA_PROJECT_DSUB_DISK_TYPE}"
echo "  Genetic sex QC       = sex-at-birth covariate, require WGS ploidy concordance=${GENETIC_SEX_REQUIRE_PLOIDY_CONCORDANCE}"
echo "  Sample QC exclusions = exclude identical-genotype components with size >= ${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}"
echo "  GWAS genotype QC     = Step1 geno<=${GWAS_STEP1_GENO_MAX} in classified EUR, AF diff<=${GWAS_STEP1_AF_DIFF_MAX} and MAF>=${GWAS_STEP1_MAF_MIN} in fit_pca_iids; Step2 geno<=${GWAS_STEP2_GENO_MAX} in classified EUR, AF diff<=${GWAS_STEP2_AF_DIFF_MAX} and MAF>=${GWAS_STEP2_MAF_MIN} in fit_pca_iids"
echo "  Height setup         = concept/source/type/unit ${HEIGHT_MEASUREMENT_CONCEPT_ID}/${HEIGHT_MEASUREMENT_SOURCE_CONCEPT_ID}/${HEIGHT_MEASUREMENT_TYPE_CONCEPT_ID}/${HEIGHT_UNIT_CONCEPT_ID}, min ${HEIGHT_MIN_CM} cm, PCs ${HEIGHT_N_PCS}"
echo "  Height REGENIE       = RUN_HEIGHT_GWAS=${RUN_HEIGHT_GWAS}, input=${HEIGHT_GWAS_INPUT_NAME}, output=${HEIGHT_GWAS_OUTPUT_NAME}, chroms=${REGENIE_CHROMS}, RINT=${REGENIE_APPLY_RINT}, blocks ${REGENIE_STEP1_BLOCK_SIZE}/${REGENIE_STEP2_BLOCK_SIZE}"
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
    echo "=== Step 6: Classify European ancestry and optional AoU comparison ==="
    bash "${SCRIPT_DIR}/compare_aou_ancestry.sh"

    # -----------------------------------------------------------------------
    # Step 7: KING kinship + close relationship classification
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 7a: Build kinship SNP subset ==="
    bash "${SCRIPT_DIR}/subset_kinship_snps.sh"

    if [[ -s "${DX_KINSHIP_DIR}/kinship_snp_subset_summary.tsv" ]]; then
        echo "  Kinship SNP subset summary:"
        awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${DX_KINSHIP_DIR}/kinship_snp_subset_summary.tsv"
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

    # -----------------------------------------------------------------------
    # Step 8: PCA European training sample selection
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 8: Select PCA European IIDs ==="
    bash "${SCRIPT_DIR}/select_pca_europeans.sh"

    # -----------------------------------------------------------------------
    # Step 9: PCA SNP QC
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 9: PCA SNP QC ==="
    bash "${SCRIPT_DIR}/pca_snp_qc.sh"

    # -----------------------------------------------------------------------
    # Step 10: PCA fitting and all-sample projection
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 10: Fit PCA and project all samples ==="
    bash "${SCRIPT_DIR}/fit_project_pca.sh"

    # -----------------------------------------------------------------------
    # Step 11: Sex covariate + sex/ploidy QC
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 11: Build sex covariate and sex/ploidy QC ==="
    bash "${SCRIPT_DIR}/get_genetic_sex.sh"

    # -----------------------------------------------------------------------
    # Step 12: Sample-QC exclusions from anomalous identical-genotype components
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 12: Build sample-QC exclusions ==="
    bash "${SCRIPT_DIR}/build_identical_component_sample_qc.sh"

    # -----------------------------------------------------------------------
    # Step 13: Final GWAS genotype inputs
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 13: Build final GWAS genotype inputs ==="
    bash "${SCRIPT_DIR}/make_gwas_genotype_inputs.sh"

    # -----------------------------------------------------------------------
    # Step 14: Height GWAS input setup
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 14: Set up height GWAS example ==="
    bash "${SCRIPT_DIR}/setup_height_gwas.sh"

    # -----------------------------------------------------------------------
    # Step 15: Optional height GWAS with REGENIE
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Step 15: Run height GWAS example ==="
    if [[ "${RUN_HEIGHT_GWAS}" == "1" ]]; then
        bash "${SCRIPT_DIR}/run_continuous_regenie_gwas.sh" \
            "${HEIGHT_GWAS_INPUT_NAME}" "${HEIGHT_GWAS_OUTPUT_NAME}"
    else
        echo "  Skipping full REGENIE GWAS because RUN_HEIGHT_GWAS=${RUN_HEIGHT_GWAS}."
        echo "  To launch it, rerun with RUN_HEIGHT_GWAS=1."
    fi
fi

echo ""
echo "=== Pipeline complete ==="
