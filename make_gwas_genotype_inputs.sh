#!/bin/bash
# make_gwas_genotype_inputs.sh - Build final genotype inputs for REGENIE GWAS.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${SBAYESRC_LIFTOVER_FILE:?SBAYESRC_LIFTOVER_FILE not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_URI:?DX_HQ_DIRECT_BFILE_URI not set}"
: "${DX_WGS_PFILE_DIR:?DX_WGS_PFILE_DIR not set}"
: "${DX_WGS_PFILE_URI:?DX_WGS_PFILE_URI not set}"
: "${DX_EUROPEANS_DIR:?DX_EUROPEANS_DIR not set}"
: "${DX_EUROPEANS_URI:?DX_EUROPEANS_URI not set}"
: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_PCA_EUR_URI:?DX_PCA_EUR_URI not set}"
: "${DX_GWAS_GENOTYPES_DIR:?DX_GWAS_GENOTYPES_DIR not set}"
: "${DX_GWAS_GENOTYPES_URI:?DX_GWAS_GENOTYPES_URI not set}"
: "${DX_GWAS_STEP1_BFILE_DIR:?DX_GWAS_STEP1_BFILE_DIR not set}"
: "${DX_GWAS_STEP1_BFILE_URI:?DX_GWAS_STEP1_BFILE_URI not set}"
: "${DX_GWAS_STEP2_PFILE_DIR:?DX_GWAS_STEP2_PFILE_DIR not set}"
: "${DX_GWAS_STEP2_PFILE_URI:?DX_GWAS_STEP2_PFILE_URI not set}"
: "${LOCAL_GWAS_GENOTYPES_DIR:?LOCAL_GWAS_GENOTYPES_DIR not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

GWAS_STEP1_GENO_MAX="${GWAS_STEP1_GENO_MAX:-0.01}"
GWAS_STEP1_AF_DIFF_MAX="${GWAS_STEP1_AF_DIFF_MAX:-0.03}"
GWAS_STEP1_MAF_MIN="${GWAS_STEP1_MAF_MIN:-0.007}"
GWAS_STEP2_GENO_MAX="${GWAS_STEP2_GENO_MAX:-0.03}"
GWAS_STEP2_AF_DIFF_MAX="${GWAS_STEP2_AF_DIFF_MAX:-0.04}"
GWAS_STEP2_MAF_MIN="${GWAS_STEP2_MAF_MIN:-0.007}"

GWAS_METRICS_DSUB_MIN_CORES="${GWAS_METRICS_DSUB_MIN_CORES:-4}"
GWAS_METRICS_DSUB_MIN_RAM="${GWAS_METRICS_DSUB_MIN_RAM:-24}"
GWAS_METRICS_DSUB_DISK_SIZE="${GWAS_METRICS_DSUB_DISK_SIZE:-150}"
GWAS_METRICS_DSUB_DISK_TYPE="${GWAS_METRICS_DSUB_DISK_TYPE:-pd-ssd}"
GWAS_DIRECT_DSUB_MIN_CORES="${GWAS_DIRECT_DSUB_MIN_CORES:-8}"
GWAS_DIRECT_DSUB_MIN_RAM="${GWAS_DIRECT_DSUB_MIN_RAM:-32}"
GWAS_DIRECT_DSUB_DISK_SIZE="${GWAS_DIRECT_DSUB_DISK_SIZE:-200}"
GWAS_DIRECT_DSUB_DISK_TYPE="${GWAS_DIRECT_DSUB_DISK_TYPE:-pd-ssd}"
GWAS_WGS_DSUB_MIN_CORES="${GWAS_WGS_DSUB_MIN_CORES:-4}"
GWAS_WGS_DSUB_MIN_RAM="${GWAS_WGS_DSUB_MIN_RAM:-24}"
GWAS_WGS_DSUB_DISK_SIZE="${GWAS_WGS_DSUB_DISK_SIZE:-180}"
GWAS_WGS_DSUB_DISK_TYPE="${GWAS_WGS_DSUB_DISK_TYPE:-pd-ssd}"

mkdir -p \
    "${LOCAL_GWAS_GENOTYPES_DIR}" \
    "${DX_GWAS_GENOTYPES_DIR}" \
    "${DX_GWAS_GENOTYPES_DIR}/metrics/direct" \
    "${DX_GWAS_GENOTYPES_DIR}/metrics/wgs" \
    "${DX_GWAS_STEP1_BFILE_DIR}" \
    "${DX_GWAS_STEP2_PFILE_DIR}"

is_uint() {
    [[ "${1:-}" =~ ^[0-9]+$ ]]
}

hq_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
hq_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"
fit_pca_iids="${DX_PCA_EUR_DIR}/fit_pca_iids.txt"
eur_iids="${DX_EUROPEANS_DIR}/classified_european_iids.txt"

for ext in bed bim fam; do
    [[ -s "${hq_prefix}.${ext}" ]] || { echo "ERROR: missing ${hq_prefix}.${ext}" >&2; exit 1; }
done
for input in "${fit_pca_iids}" "${eur_iids}" "${SBAYESRC_LIFTOVER_FILE}"; do
    [[ -s "${input}" ]] || { echo "ERROR: missing ${input}" >&2; exit 1; }
done
for c in $(seq 1 22); do
    for ext in pgen pvar psam summary.tsv; do
        [[ -s "${DX_WGS_PFILE_DIR}/chr${c}.${ext}" ]] || { echo "ERROR: missing ${DX_WGS_PFILE_DIR}/chr${c}.${ext}" >&2; exit 1; }
    done
done

direct_source_variants=$(wc -l < "${hq_prefix}.bim")
direct_source_samples=$(wc -l < "${hq_prefix}.fam")
fit_pca_samples=$(wc -l < "${fit_pca_iids}")
eur_samples=$(wc -l < "${eur_iids}")

desired_params="${LOCAL_GWAS_GENOTYPES_DIR}/gwas_genotype_qc.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'step1_source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'step1_geno_max_classified_eur\t%s\n' "${GWAS_STEP1_GENO_MAX}"
    printf 'step1_af_diff_max_fit_pca_vs_snpinfo\t%s\n' "${GWAS_STEP1_AF_DIFF_MAX}"
    printf 'step1_maf_min_fit_pca\t%s\n' "${GWAS_STEP1_MAF_MIN}"
    printf 'step2_source_pfiles\t%s\n' "wgs_pfiles/chr1-22"
    printf 'step2_geno_max_classified_eur\t%s\n' "${GWAS_STEP2_GENO_MAX}"
    printf 'step2_af_diff_max_fit_pca_vs_snpinfo\t%s\n' "${GWAS_STEP2_AF_DIFF_MAX}"
    printf 'step2_maf_min_fit_pca\t%s\n' "${GWAS_STEP2_MAF_MIN}"
    printf 'classified_european_samples\t%s\n' "${eur_samples}"
    printf 'classified_european_iids_size\t%s\n' "$(stat -c%s "${eur_iids}")"
    printf 'fit_pca_samples\t%s\n' "${fit_pca_samples}"
    printf 'fit_pca_iids_size\t%s\n' "$(stat -c%s "${fit_pca_iids}")"
    printf 'direct_source_variants\t%s\n' "${direct_source_variants}"
    printf 'direct_source_samples\t%s\n' "${direct_source_samples}"
    printf 'liftover_size\t%s\n' "$(stat -c%s "${SBAYESRC_LIFTOVER_FILE}")"
    printf 'plink2_version\t%s\n' "$("${PLINK2}" --version 2>&1 | head -1)"
    for c in $(seq 1 22); do
        variants=$(awk -F'\t' 'NR == 2 {print $6; exit}' "${DX_WGS_PFILE_DIR}/chr${c}.summary.tsv")
        printf 'chr%s_source_variants\t%s\n' "${c}" "${variants}"
        printf 'chr%s_pgen_size\t%s\n' "${c}" "$(stat -c%s "${DX_WGS_PFILE_DIR}/chr${c}.pgen")"
        printf 'chr%s_pvar_size\t%s\n' "${c}" "$(stat -c%s "${DX_WGS_PFILE_DIR}/chr${c}.pvar")"
        printf 'chr%s_psam_size\t%s\n' "${c}" "$(stat -c%s "${DX_WGS_PFILE_DIR}/chr${c}.psam")"
    done
} > "${desired_params}"

params="${DX_GWAS_GENOTYPES_DIR}/gwas_genotype_qc.params.tsv"
summary="${DX_GWAS_GENOTYPES_DIR}/gwas_genotype_qc.summary.tsv"
step1_prefix="${DX_GWAS_STEP1_BFILE_DIR}/chr1_22_merged_gwas_step1"

outputs_complete=1
if [[ ! -s "${params}" || ! -s "${summary}" || ! -s "${step1_prefix}.bed" || ! -s "${step1_prefix}.bim" || ! -s "${step1_prefix}.fam" ]]; then
    outputs_complete=0
elif ! diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
    outputs_complete=0
else
    expected_step1=$(awk -F'\t' '$1 == "step1_direct" {print $7; exit}' "${summary}")
    observed_step1=$(wc -l < "${step1_prefix}.bim")
    [[ -n "${expected_step1}" && "${observed_step1}" -eq "${expected_step1}" ]] || outputs_complete=0
    for c in $(seq 1 22); do
        p="${DX_GWAS_STEP2_PFILE_DIR}/chr${c}"
        if [[ ! -s "${p}.pgen" || ! -s "${p}.pvar" || ! -s "${p}.psam" || ! -s "${p}.summary.tsv" ]]; then
            outputs_complete=0
            break
        fi
        expected=$(awk -F'\t' -v label="chr${c}" '$1 == label {print $7; exit}' "${summary}")
        observed=$(grep -vc '^#' "${p}.pvar")
        [[ -n "${expected}" && "${observed}" -eq "${expected}" ]] || { outputs_complete=0; break; }
    done
fi
if [[ "${outputs_complete}" -eq 1 ]]; then
    echo "  GWAS genotype inputs already exist and match current parameters — skipping"
    exit 0
fi

echo "  Staging plink2 and desired params ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
gcloud storage cp "${desired_params}" "${DX_GWAS_GENOTYPES_URI}/gwas_genotype_qc.desired_params.tsv" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

# ---------------------------------------------------------------------------
# Metric scans
# ---------------------------------------------------------------------------
direct_metrics_ok=0
direct_metrics_summary="${DX_GWAS_GENOTYPES_DIR}/metrics/direct/direct_hq.metrics_summary.tsv"
if [[ -s "${DX_GWAS_GENOTYPES_DIR}/metrics/direct/direct_hq.fit_pca.acount" &&
      -s "${DX_GWAS_GENOTYPES_DIR}/metrics/direct/direct_hq.our_eur.vmiss" &&
      -s "${direct_metrics_summary}" ]]; then
    fit_lines=$(awk -F'\t' '$1 == "fit_pca_acount_lines" {print $2; exit}' "${direct_metrics_summary}")
    miss_lines=$(awk -F'\t' '$1 == "our_eur_vmiss_lines" {print $2; exit}' "${direct_metrics_summary}")
    metric_fit_samples=$(awk -F'\t' '$1 == "fit_pca_keep_samples" {print $2; exit}' "${direct_metrics_summary}")
    metric_eur_samples=$(awk -F'\t' '$1 == "eur_keep_samples" {print $2; exit}' "${direct_metrics_summary}")
    if is_uint "${fit_lines}" && is_uint "${miss_lines}" &&
       is_uint "${metric_fit_samples}" && is_uint "${metric_eur_samples}" &&
       [[ "${fit_lines}" -eq $((direct_source_variants + 1)) &&
          "${miss_lines}" -eq $((direct_source_variants + 1)) &&
          "${metric_fit_samples}" -eq "${fit_pca_samples}" &&
          "${metric_eur_samples}" -eq "${eur_samples}" ]]; then
        direct_metrics_ok=1
    fi
fi
if [[ "${direct_metrics_ok}" -eq 0 ]]; then
    echo "  Submitting GWAS Step 1 direct metrics job ..."
    direct_metrics_log="${SCRIPT_DIR}/logs/dsub_gwas_direct_metrics_$(date +%Y%m%d_%H%M%S).dsub.out"
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
        --name "sbayesrc-gwas-direct-metrics" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_gwas_direct_metrics_worker.sh" \
        --input PLINK2="${DSUB_PLINK2_GS}" \
        --input BED="${hq_uri_prefix}.bed" \
        --input BIM="${hq_uri_prefix}.bim" \
        --input FAM="${hq_uri_prefix}.fam" \
        --input FIT_PCA_KEEP="${DX_PCA_EUR_URI}/fit_pca_iids.txt" \
        --input EUR_KEEP="${DX_EUROPEANS_URI}/classified_european_iids.txt" \
        --output-recursive OUTDIR="${DX_GWAS_GENOTYPES_URI}/metrics/direct/" \
        --min-cores "${GWAS_DIRECT_DSUB_MIN_CORES}" \
        --min-ram "${GWAS_DIRECT_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${GWAS_DIRECT_DSUB_DISK_SIZE}" \
        --disk-type "${GWAS_DIRECT_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${direct_metrics_log}"
    dsub_rc=${PIPESTATUS[0]}
    [[ "${dsub_rc}" -eq 0 ]] || exit "${dsub_rc}"
else
    echo "  Direct GWAS metrics already exist — skipping metrics job"
fi

wgs_tasks="${SCRIPT_DIR}/logs/dsub_gwas_wgs_metrics_$(date +%Y%m%d_%H%M%S).tsv"
{
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '--env CHROM' '--input PLINK2' '--input PGEN' '--input PVAR' '--input PSAM' \
        '--input FIT_PCA_KEEP' '--input EUR_KEEP' '--output-recursive OUTDIR'
    for c in $(seq 1 22); do
        chrom="chr${c}"
        source_variants=$(awk -F'\t' -v key="chr${c}_source_variants" '$1 == key {print $2; exit}' "${desired_params}")
        metric_summary="${DX_GWAS_GENOTYPES_DIR}/metrics/wgs/${chrom}.metrics_summary.tsv"
        submit=1
        if [[ -s "${DX_GWAS_GENOTYPES_DIR}/metrics/wgs/${chrom}.fit_pca.acount" &&
              -s "${DX_GWAS_GENOTYPES_DIR}/metrics/wgs/${chrom}.our_eur.vmiss" &&
              -s "${metric_summary}" ]]; then
            fit_lines=$(awk -F'\t' '$2 == "fit_pca_acount_lines" {print $3; exit}' "${metric_summary}")
            miss_lines=$(awk -F'\t' '$2 == "our_eur_vmiss_lines" {print $3; exit}' "${metric_summary}")
            metric_fit_samples=$(awk -F'\t' '$2 == "fit_pca_keep_samples" {print $3; exit}' "${metric_summary}")
            metric_eur_samples=$(awk -F'\t' '$2 == "eur_keep_samples" {print $3; exit}' "${metric_summary}")
            if is_uint "${source_variants}" && is_uint "${fit_lines}" && is_uint "${miss_lines}" &&
               is_uint "${metric_fit_samples}" && is_uint "${metric_eur_samples}" &&
               [[ "${fit_lines}" -eq $((source_variants + 1)) &&
                  "${miss_lines}" -eq $((source_variants + 1)) &&
                  "${metric_fit_samples}" -eq "${fit_pca_samples}" &&
                  "${metric_eur_samples}" -eq "${eur_samples}" ]]; then
                submit=0
            fi
        fi
        if [[ "${submit}" -eq 1 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${chrom}" \
                "${DSUB_PLINK2_GS}" \
                "${DX_WGS_PFILE_URI}/${chrom}.pgen" \
                "${DX_WGS_PFILE_URI}/${chrom}.pvar" \
                "${DX_WGS_PFILE_URI}/${chrom}.psam" \
                "${DX_PCA_EUR_URI}/fit_pca_iids.txt" \
                "${DX_EUROPEANS_URI}/classified_european_iids.txt" \
                "${DX_GWAS_GENOTYPES_URI}/metrics/wgs/"
        fi
    done
} > "${wgs_tasks}"
wgs_metric_tasks=$(( $(wc -l < "${wgs_tasks}") - 1 ))
if [[ "${wgs_metric_tasks}" -gt 0 ]]; then
    echo "  Submitting ${wgs_metric_tasks} WGS metric task(s) ..."
    wgs_metrics_log="${wgs_tasks%.tsv}.dsub.out"
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
        --name "sbayesrc-gwas-wgs-metrics" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_gwas_wgs_metrics_worker.sh" \
        --tasks "${wgs_tasks}" \
        --min-cores "${GWAS_METRICS_DSUB_MIN_CORES}" \
        --min-ram "${GWAS_METRICS_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${GWAS_METRICS_DSUB_DISK_SIZE}" \
        --disk-type "${GWAS_METRICS_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${wgs_metrics_log}"
    dsub_rc=${PIPESTATUS[0]}
    [[ "${dsub_rc}" -eq 0 ]] || exit "${dsub_rc}"
else
    echo "  WGS GWAS metrics already exist — skipping metrics jobs"
fi

# ---------------------------------------------------------------------------
# Build extract lists and frequency/QC summaries locally.
# ---------------------------------------------------------------------------
echo "  Building GWAS genotype filter extract lists and fit_pca allele-frequency files ..."
filter_output_dir="${LOCAL_GWAS_GENOTYPES_DIR}/filter_outputs"
rm -rf "${filter_output_dir}"
mkdir -p "${filter_output_dir}"
python3 "${SCRIPT_DIR}/build_gwas_genotype_filters.py" \
    --liftover "${SBAYESRC_LIFTOVER_FILE}" \
    --metrics-dir "${DX_GWAS_GENOTYPES_DIR}/metrics" \
    --output-dir "${filter_output_dir}" \
    --step1-af-diff-max "${GWAS_STEP1_AF_DIFF_MAX}" \
    --step1-maf-min "${GWAS_STEP1_MAF_MIN}" \
    --step1-geno-max "${GWAS_STEP1_GENO_MAX}" \
    --step2-af-diff-max "${GWAS_STEP2_AF_DIFF_MAX}" \
    --step2-maf-min "${GWAS_STEP2_MAF_MIN}" \
    --step2-geno-max "${GWAS_STEP2_GENO_MAX}"
cp "${desired_params}" "${filter_output_dir}/gwas_genotype_qc.params.tsv"
echo "  Uploading completed filter outputs to ${DX_GWAS_GENOTYPES_URI}/ ..."
gcloud storage cp "${filter_output_dir}/gwas_genotype_qc.summary.tsv" \
    "${DX_GWAS_GENOTYPES_URI}/gwas_genotype_qc.summary.tsv" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
gcloud storage cp "${filter_output_dir}/gwas_genotype_qc.params.tsv" \
    "${DX_GWAS_GENOTYPES_URI}/gwas_genotype_qc.params.tsv" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
gcloud storage rsync -r "${filter_output_dir}/step1_direct" \
    "${DX_GWAS_GENOTYPES_URI}/step1_direct" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
gcloud storage rsync -r "${filter_output_dir}/step2_wgs" \
    "${DX_GWAS_GENOTYPES_URI}/step2_wgs" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

expected_step1=$(awk -F'\t' '$1 == "step1_direct" {print $7; exit}' "${summary}")
if ! is_uint "${expected_step1}" || [[ "${expected_step1}" -le 0 ]]; then
    echo "ERROR: could not determine expected Step 1 variant count from ${summary}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build final Step 1 bfile.
# ---------------------------------------------------------------------------
step1_run=1
if [[ -s "${step1_prefix}.bed" && -s "${step1_prefix}.bim" && -s "${step1_prefix}.fam" ]]; then
    observed=$(wc -l < "${step1_prefix}.bim")
    [[ "${observed}" -eq "${expected_step1}" ]] && step1_run=0
fi
if [[ "${step1_run}" -eq 1 ]]; then
    echo "  Submitting final Step 1 direct bfile job (${expected_step1} variants) ..."
    step1_log="${SCRIPT_DIR}/logs/dsub_gwas_step1_direct_$(date +%Y%m%d_%H%M%S).dsub.out"
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
        --name "sbayesrc-gwas-step1-direct" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_gwas_step1_direct_worker.sh" \
        --env EXPECTED_VARIANTS="${expected_step1}" \
        --input PLINK2="${DSUB_PLINK2_GS}" \
        --input BED="${hq_uri_prefix}.bed" \
        --input BIM="${hq_uri_prefix}.bim" \
        --input FAM="${hq_uri_prefix}.fam" \
        --input EXTRACT="${DX_GWAS_GENOTYPES_URI}/step1_direct/chr1_22_merged_gwas_step1.extract.txt" \
        --input PARAMS="${DX_GWAS_GENOTYPES_URI}/gwas_genotype_qc.params.tsv" \
        --output-recursive OUTDIR="${DX_GWAS_STEP1_BFILE_URI}/" \
        --min-cores "${GWAS_DIRECT_DSUB_MIN_CORES}" \
        --min-ram "${GWAS_DIRECT_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${GWAS_DIRECT_DSUB_DISK_SIZE}" \
        --disk-type "${GWAS_DIRECT_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${step1_log}"
    dsub_rc=${PIPESTATUS[0]}
    [[ "${dsub_rc}" -eq 0 ]] || exit "${dsub_rc}"
else
    echo "  Final Step 1 direct bfile already exists (${expected_step1} variants) — skipping"
fi

# ---------------------------------------------------------------------------
# Build final Step 2 per-chromosome WGS pfiles.
# ---------------------------------------------------------------------------
wgs_build_tasks="${SCRIPT_DIR}/logs/dsub_gwas_step2_wgs_$(date +%Y%m%d_%H%M%S).tsv"
{
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '--env CHROM' '--env EXPECTED_VARIANTS' '--input PLINK2' '--input PGEN' '--input PVAR' '--input PSAM' '--input EXTRACT' '--output-recursive OUTDIR'
    for c in $(seq 1 22); do
        chrom="chr${c}"
        expected=$(awk -F'\t' -v label="${chrom}" '$1 == label {print $7; exit}' "${summary}")
        if ! is_uint "${expected}"; then
            echo "ERROR: could not determine expected ${chrom} variant count from ${summary}" >&2
            exit 1
        fi
        p="${DX_GWAS_STEP2_PFILE_DIR}/${chrom}"
        submit=1
        if [[ -s "${p}.pgen" && -s "${p}.pvar" && -s "${p}.psam" && -s "${p}.summary.tsv" ]]; then
            observed=$(grep -vc '^#' "${p}.pvar")
            [[ "${observed}" -eq "${expected}" ]] && submit=0
        fi
        if [[ "${submit}" -eq 1 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${chrom}" "${expected}" "${DSUB_PLINK2_GS}" \
                "${DX_WGS_PFILE_URI}/${chrom}.pgen" \
                "${DX_WGS_PFILE_URI}/${chrom}.pvar" \
                "${DX_WGS_PFILE_URI}/${chrom}.psam" \
                "${DX_GWAS_GENOTYPES_URI}/step2_wgs/extracts/${chrom}.extract.txt" \
                "${DX_GWAS_STEP2_PFILE_URI}/"
        fi
    done
} > "${wgs_build_tasks}"
wgs_build_count=$(( $(wc -l < "${wgs_build_tasks}") - 1 ))
if [[ "${wgs_build_count}" -gt 0 ]]; then
    echo "  Submitting ${wgs_build_count} final Step 2 WGS pfile task(s) ..."
    wgs_build_log="${wgs_build_tasks%.tsv}.dsub.out"
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
        --name "sbayesrc-gwas-step2-wgs" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_gwas_step2_wgs_worker.sh" \
        --tasks "${wgs_build_tasks}" \
        --input PARAMS="${DX_GWAS_GENOTYPES_URI}/gwas_genotype_qc.params.tsv" \
        --min-cores "${GWAS_WGS_DSUB_MIN_CORES}" \
        --min-ram "${GWAS_WGS_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${GWAS_WGS_DSUB_DISK_SIZE}" \
        --disk-type "${GWAS_WGS_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${wgs_build_log}"
    dsub_rc=${PIPESTATUS[0]}
    [[ "${dsub_rc}" -eq 0 ]] || exit "${dsub_rc}"
else
    echo "  Final Step 2 WGS pfiles already exist — skipping"
fi

echo "  Verifying GWAS genotype outputs ..."
observed_step1=$(wc -l < "${step1_prefix}.bim")
[[ "${observed_step1}" -eq "${expected_step1}" ]] || { echo "ERROR: Step 1 observed ${observed_step1}, expected ${expected_step1}" >&2; exit 1; }
for c in $(seq 1 22); do
    chrom="chr${c}"
    expected=$(awk -F'\t' -v label="${chrom}" '$1 == label {print $7; exit}' "${summary}")
    observed=$(grep -vc '^#' "${DX_GWAS_STEP2_PFILE_DIR}/${chrom}.pvar")
    [[ "${observed}" -eq "${expected}" ]] || { echo "ERROR: ${chrom} observed ${observed}, expected ${expected}" >&2; exit 1; }
done

echo "  Done. GWAS genotype summary:"
sed 's/^/    /' "${summary}"
