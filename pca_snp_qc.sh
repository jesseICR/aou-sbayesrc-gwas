#!/bin/bash
# pca_snp_qc.sh - QC SNPs for PCA and build the PCA-ready bfile.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_URI:?DX_HQ_DIRECT_BFILE_URI not set}"
: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_PCA_EUR_URI:?DX_PCA_EUR_URI not set}"
: "${LOCAL_PCA_QC_DIR:?LOCAL_PCA_QC_DIR not set}"
: "${PCA_HIGH_LD_URL:?PCA_HIGH_LD_URL not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

PCA_AF_DIFF_MAX="${PCA_AF_DIFF_MAX:-0.03}"
PCA_MAF_MIN="${PCA_MAF_MIN:-0.01}"
PCA_GENO_MAX="${PCA_GENO_MAX:-0.01}"
PCA_MIND_MAX="${PCA_MIND_MAX:-0.01}"
PCA_LD_WINDOW="${PCA_LD_WINDOW:-1000}"
PCA_LD_STEP="${PCA_LD_STEP:-80}"
PCA_LD_R2="${PCA_LD_R2:-0.1}"
PCA_SNP_QC_DSUB_MIN_CORES="${PCA_SNP_QC_DSUB_MIN_CORES:-8}"
PCA_SNP_QC_DSUB_MIN_RAM="${PCA_SNP_QC_DSUB_MIN_RAM:-32}"
PCA_SNP_QC_DSUB_DISK_SIZE="${PCA_SNP_QC_DSUB_DISK_SIZE:-300}"
PCA_SNP_QC_DSUB_DISK_TYPE="${PCA_SNP_QC_DSUB_DISK_TYPE:-pd-ssd}"

mkdir -p "${LOCAL_PCA_QC_DIR}" "${DX_PCA_EUR_DIR}" "${DX_PCA_EUR_DIR}/scrap"

source_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
source_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"
fit_iids="${DX_PCA_EUR_DIR}/fit_pca_iids.txt"
variant_qc="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.variant_qc.tsv"

for ext in bed bim fam; do
    if [[ ! -s "${source_prefix}.${ext}" ]]; then
        echo "ERROR: missing HQ direct bfile input ${source_prefix}.${ext}" >&2
        exit 1
    fi
done
for input in "${fit_iids}" "${variant_qc}"; do
    if [[ ! -s "${input}" ]]; then
        echo "ERROR: missing PCA SNP QC input ${input}" >&2
        exit 1
    fi
done

high_ld_bed="${LOCAL_PCA_QC_DIR}/high-LD-regions-hg38-GRCh38.bed"
if [[ -s "${high_ld_bed}" ]]; then
    echo "  High-LD region BED already cached (${high_ld_bed}, $(wc -l < "${high_ld_bed}") lines)"
else
    echo "  Downloading hg38 high-LD region BED ..."
    curl -fsSL --retry 3 --retry-delay 5 -o "${high_ld_bed}" "${PCA_HIGH_LD_URL}"
    echo "  Downloaded ${high_ld_bed} ($(wc -l < "${high_ld_bed}") lines)"
fi

source_variants=$(wc -l < "${source_prefix}.bim")
source_samples=$(wc -l < "${source_prefix}.fam")
fit_samples=$(wc -l < "${fit_iids}")
variant_qc_size=$(stat -c%s "${variant_qc}")
high_ld_sha256=$(sha256sum "${high_ld_bed}" | awk '{print $1}')

desired_params="${LOCAL_PCA_QC_DIR}/pca_snp_qc.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'source_variants\t%s\n' "${source_variants}"
    printf 'source_samples\t%s\n' "${source_samples}"
    printf 'fit_pca_iids\t%s\n' "${fit_samples}"
    printf 'fit_pca_iids_size\t%s\n' "$(stat -c%s "${fit_iids}")"
    printf 'variant_qc_size\t%s\n' "${variant_qc_size}"
    printf 'pca_af_diff_max\t%s\n' "${PCA_AF_DIFF_MAX}"
    printf 'pca_maf_min\t%s\n' "${PCA_MAF_MIN}"
    printf 'pca_geno_max\t%s\n' "${PCA_GENO_MAX}"
    printf 'pca_mind_max\t%s\n' "${PCA_MIND_MAX}"
    printf 'pca_ld_window\t%s\n' "${PCA_LD_WINDOW}"
    printf 'pca_ld_step\t%s\n' "${PCA_LD_STEP}"
    printf 'pca_ld_r2\t%s\n' "${PCA_LD_R2}"
    printf 'pca_high_ld_url\t%s\n' "${PCA_HIGH_LD_URL}"
    printf 'pca_high_ld_sha256\t%s\n' "${high_ld_sha256}"
    printf 'plink2_version\t%s\n' "$("${PLINK2}" --version 2>&1 | head -1)"
} > "${desired_params}"

params="${DX_PCA_EUR_DIR}/pca_snp_qc.params.tsv"
summary="${DX_PCA_EUR_DIR}/pca_snp_qc.summary.tsv"
filter_steps="${DX_PCA_EUR_DIR}/pca_snp_qc.filter_steps.tsv"
pca_prefix="${DX_PCA_EUR_DIR}/pca_ready"
if [[ -s "${pca_prefix}.bed" && -s "${pca_prefix}.bim" && -s "${pca_prefix}.fam" &&
      -s "${summary}" && -s "${filter_steps}" && -s "${params}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected_variants=$(awk -F'\t' '$1 == "final_variants" {print $2; exit}' "${summary}")
        expected_samples=$(awk -F'\t' '$1 == "final_samples" {print $2; exit}' "${summary}")
        observed_variants=$(wc -l < "${pca_prefix}.bim")
        observed_samples=$(wc -l < "${pca_prefix}.fam")
        if [[ -n "${expected_variants}" && -n "${expected_samples}" &&
              "${observed_variants}" -eq "${expected_variants}" &&
              "${observed_samples}" -eq "${expected_samples}" ]]; then
            echo "  PCA-ready bfile already exists (${observed_variants} variants, ${observed_samples} samples) — skipping"
            exit 0
        fi
    fi
    echo "  PCA SNP QC outputs exist but params/counts do not match — rebuilding"
fi

echo "  Staging PCA SNP QC inputs to workspace bucket ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
export DSUB_PCA_QC_PREP_URI="${DSUB_PCA_QC_PREP_URI:-${DX_PCA_EUR_URI}/_prep}"
gcloud storage cp \
    "${high_ld_bed}" \
    "${desired_params}" \
    "${DSUB_PCA_QC_PREP_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

pca_log="${SCRIPT_DIR}/logs/dsub_pca_snp_qc_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${pca_log}")"

echo "  Submitting dsub job for PCA SNP QC ..."
dsub \
    --provider "${DSUB_PROVIDER}" \
    --project "${GOOGLE_PROJECT}" \
    --regions "${DSUB_REGION}" \
    --service-account "${DSUB_PET_SA}" \
    --use-private-address \
    --network "${DSUB_NETWORK}" \
    --subnetwork "${DSUB_SUBNETWORK}" \
    --user-project "${GOOGLE_PROJECT}" \
    --logging "${DSUB_LOG_URI}" \
    --name "sbayesrc-pca-snp-qc" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_pca_snp_qc_worker.sh" \
    --env PCA_AF_DIFF_MAX="${PCA_AF_DIFF_MAX}" \
    --env PCA_MAF_MIN="${PCA_MAF_MIN}" \
    --env PCA_GENO_MAX="${PCA_GENO_MAX}" \
    --env PCA_MIND_MAX="${PCA_MIND_MAX}" \
    --env PCA_LD_WINDOW="${PCA_LD_WINDOW}" \
    --env PCA_LD_STEP="${PCA_LD_STEP}" \
    --env PCA_LD_R2="${PCA_LD_R2}" \
    --env EXPECTED_SOURCE_VARIANTS="${source_variants}" \
    --env EXPECTED_SOURCE_SAMPLES="${source_samples}" \
    --env EXPECTED_FIT_SAMPLES="${fit_samples}" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input BED="${source_uri_prefix}.bed" \
    --input BIM="${source_uri_prefix}.bim" \
    --input FAM="${source_uri_prefix}.fam" \
    --input KEEP="${DX_PCA_EUR_URI}/fit_pca_iids.txt" \
    --input VARIANT_QC="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq.variant_qc.tsv" \
    --input HIGH_LD_BED="${DSUB_PCA_QC_PREP_URI}/high-LD-regions-hg38-GRCh38.bed" \
    --input PARAMS="${DSUB_PCA_QC_PREP_URI}/pca_snp_qc.desired_params.tsv" \
    --output-recursive OUTDIR="${DX_PCA_EUR_URI}/" \
    --min-cores "${PCA_SNP_QC_DSUB_MIN_CORES}" \
    --min-ram "${PCA_SNP_QC_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${PCA_SNP_QC_DSUB_DISK_SIZE}" \
    --disk-type "${PCA_SNP_QC_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${pca_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub PCA SNP QC job returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

echo "  Verifying PCA SNP QC outputs ..."
for attempt in $(seq 1 60); do
    if [[ -s "${pca_prefix}.bed" && -s "${pca_prefix}.bim" && -s "${pca_prefix}.fam" &&
          -s "${summary}" && -s "${filter_steps}" && -s "${params}" ]]; then
        expected_variants=$(awk -F'\t' '$1 == "final_variants" {print $2; exit}' "${summary}")
        expected_samples=$(awk -F'\t' '$1 == "final_samples" {print $2; exit}' "${summary}")
        observed_variants=$(wc -l < "${pca_prefix}.bim")
        observed_samples=$(wc -l < "${pca_prefix}.fam")
        if [[ "${observed_variants}" -eq "${expected_variants}" &&
              "${observed_samples}" -eq "${expected_samples}" ]]; then
            echo "  Done: ${pca_prefix}.{bed,bim,fam} (${observed_variants} variants, ${observed_samples} samples)"
            echo "  PCA SNP QC filter steps:"
            sed 's/^/    /' "${filter_steps}"
            exit 0
        fi
        echo "  Outputs visible but counts do not match yet; waiting ..."
    else
        echo "  Waiting for PCA SNP QC outputs to appear (${attempt}/60) ..."
    fi
    sleep 10
done

echo "ERROR: PCA SNP QC outputs were not visible with expected counts after dsub completed" >&2
exit 1
