#!/bin/bash
# make_hapmap3_bfile_hq.sh - Optional HapMap3 HQ bfile from WGS pfiles.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set — are you running inside an AoU Verily Jupyter session?}"

WORKSPACE_BUCKET_MOUNT="${WORKSPACE_BUCKET_MOUNT:-/home/jupyter/workspace/workspace-bucket}"

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

    echo "Workspace bucket is not mounted at ${WORKSPACE_BUCKET_MOUNT}; attempting Workbench mount ..."
    if ! command -v wb >/dev/null 2>&1; then
        echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is not a writable gcsfuse mount, and 'wb' is unavailable." >&2
        exit 1
    fi
    if ! run_wb resource mount --allow-other; then
        echo "ERROR: 'wb resource mount' failed." >&2
        exit 1
    fi
    if ! is_workspace_bucket_fuse_mounted || [[ ! -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
        echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is still not a writable gcsfuse mount." >&2
        echo "  Do not mkdir this path manually; it must be durable GCS-backed storage." >&2
        exit 1
    fi
}

is_uint() {
    [[ "${1:-}" =~ ^[0-9]+$ ]]
}

require_file() {
    local path="$1"
    if [[ ! -s "${path}" ]]; then
        echo "ERROR: missing required file ${path}" >&2
        exit 1
    fi
}

metric_value() {
    local file="$1"
    local key="$2"
    awk -F'\t' -v key="${key}" '$2 == key {print $3; exit}' "${file}"
}

ensure_workspace_bucket_mount
WORKSPACE_BUCKET_URI="gs://$(mount | awk -v target="${WORKSPACE_BUCKET_MOUNT}" '$2 == "on" && $3 == target {print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    exit 1
fi

export PLINK2="${PLINK2:-/opt/workbench-tools/binaries/bin/plink2}"
export DSUB_PROVIDER="${DSUB_PROVIDER:-google-batch}"
export DSUB_REGION="${DSUB_REGION:-us-central1}"
export DSUB_NETWORK="${DSUB_NETWORK:-projects/${GOOGLE_PROJECT}/global/networks/network}"
export DSUB_SUBNETWORK="${DSUB_SUBNETWORK:-projects/${GOOGLE_PROJECT}/regions/${DSUB_REGION}/subnetworks/subnetwork}"
export DSUB_IMAGE="${DSUB_IMAGE:-marketplace.gcr.io/google/ubuntu2204}"
DSUB_PET_SA="${DSUB_PET_SA:-$(gcloud config get-value account 2>/dev/null || true)}"
if [[ -z "${DSUB_PET_SA}" ]]; then
    echo "ERROR: could not determine the pod pet service account via gcloud config." >&2
    exit 1
fi
export DSUB_PET_SA

DSUB_BOOT_DISK_SIZE="${DSUB_BOOT_DISK_SIZE:-50}"
AOU_DATA_VERSION="${AOU_DATA_VERSION:-v9}"
SBAYESRC_OUTPUT_PREFIX="${SBAYESRC_OUTPUT_PREFIX:-sbayesrc_genotypes}"
DSUB_PLINK2_GS="${DSUB_PLINK2_GS:-${WORKSPACE_BUCKET_URI}/bin/plink2}"
DSUB_LOG_URI="${DSUB_LOG_URI:-${WORKSPACE_BUCKET_URI}/${SBAYESRC_OUTPUT_PREFIX}/logs/dsub}"

HAPMAP3_AF_DIFF_MAX="${HAPMAP3_AF_DIFF_MAX:-0.03}"
HAPMAP3_MAF_MIN="${HAPMAP3_MAF_MIN:-0.007}"
HAPMAP3_GENO_MAX="${HAPMAP3_GENO_MAX:-0.01}"
HAPMAP3_EXTRACT_DSUB_MIN_CORES="${HAPMAP3_EXTRACT_DSUB_MIN_CORES:-4}"
HAPMAP3_EXTRACT_DSUB_MIN_RAM="${HAPMAP3_EXTRACT_DSUB_MIN_RAM:-24}"
HAPMAP3_EXTRACT_DSUB_DISK_SIZE="${HAPMAP3_EXTRACT_DSUB_DISK_SIZE:-180}"
HAPMAP3_EXTRACT_DSUB_DISK_TYPE="${HAPMAP3_EXTRACT_DSUB_DISK_TYPE:-pd-ssd}"
HAPMAP3_MERGE_DSUB_MIN_CORES="${HAPMAP3_MERGE_DSUB_MIN_CORES:-8}"
HAPMAP3_MERGE_DSUB_MIN_RAM="${HAPMAP3_MERGE_DSUB_MIN_RAM:-32}"
HAPMAP3_MERGE_DSUB_DISK_SIZE="${HAPMAP3_MERGE_DSUB_DISK_SIZE:-500}"
HAPMAP3_MERGE_DSUB_DISK_TYPE="${HAPMAP3_MERGE_DSUB_DISK_TYPE:-pd-ssd}"

DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/${SBAYESRC_OUTPUT_PREFIX}"
DX_OUTPUT_URI="${WORKSPACE_BUCKET_URI}/${SBAYESRC_OUTPUT_PREFIX}"
DX_WGS_PFILE_DIR="${DX_OUTPUT_DIR}/wgs_pfiles"
DX_WGS_PFILE_URI="${DX_OUTPUT_URI}/wgs_pfiles"
DX_GWAS_METRICS_DIR="${DX_OUTPUT_DIR}/gwas_genotypes/metrics"
DX_GWAS_METRICS_URI="${DX_OUTPUT_URI}/gwas_genotypes/metrics"
DX_EUROPEANS_DIR="${DX_OUTPUT_DIR}/europeans"
DX_EUROPEANS_URI="${DX_OUTPUT_URI}/europeans"
DX_PCA_EUR_DIR="${DX_OUTPUT_DIR}/pca_eur"
DX_HAPMAP3_DIR="${DX_OUTPUT_DIR}/hapmap3_bfile_hq"
DX_HAPMAP3_URI="${DX_OUTPUT_URI}/hapmap3_bfile_hq"

LOCAL_HAPMAP3_DIR="${SCRIPT_DIR}/data/hapmap3_bfile_hq"
HAPMAP3_RSIDS_FILE="${HAPMAP3_RSIDS_FILE:-${SCRIPT_DIR}/hapmap3_rsids.txt}"
HAPMAP3_EXPECTED_ROWS="${HAPMAP3_EXPECTED_ROWS:-1154522}"
HAPMAP3_EXPECTED_SHA256="${HAPMAP3_EXPECTED_SHA256:-508e1b1739484b52af24b51f58aea833cf01f186b17696f0683ecd2abc687087}"
SBAYESRC_LIFTOVER_FILE="${SBAYESRC_LIFTOVER_FILE:-${SCRIPT_DIR}/data/support/sbayesrc_liftover_results.csv}"
EUR_KEEP="${DX_EUROPEANS_DIR}/classified_european_iids.txt"
FIT_PCA_KEEP="${DX_PCA_EUR_DIR}/fit_pca_iids.txt"

mkdir -p "${SCRIPT_DIR}/logs" "${LOCAL_HAPMAP3_DIR}" "${DX_HAPMAP3_DIR}" "${DX_HAPMAP3_DIR}/pfiles"

for f in "${PLINK2}" "${HAPMAP3_RSIDS_FILE}" "${SBAYESRC_LIFTOVER_FILE}" "${EUR_KEEP}" "${FIT_PCA_KEEP}"; do
    require_file "${f}"
done
for c in $(seq 1 22); do
    for ext in pgen pvar psam; do
        require_file "${DX_WGS_PFILE_DIR}/chr${c}.${ext}"
    done
    require_file "${DX_GWAS_METRICS_DIR}/wgs/chr${c}.fit_pca.acount"
    require_file "${DX_GWAS_METRICS_DIR}/wgs/chr${c}.our_eur.vmiss"
    require_file "${DX_GWAS_METRICS_DIR}/wgs/chr${c}.metrics_summary.tsv"
done

hapmap_rows=$(wc -l < "${HAPMAP3_RSIDS_FILE}")
hapmap_sha=$(sha256sum "${HAPMAP3_RSIDS_FILE}" | awk '{print $1}')
hapmap_duplicates=$(awk 'NF {seen[$1]++} END {d=0; for (id in seen) if (seen[id] > 1) d += seen[id] - 1; print d + 0}' "${HAPMAP3_RSIDS_FILE}")
if [[ "${hapmap_rows}" -ne "${HAPMAP3_EXPECTED_ROWS}" ]]; then
    echo "ERROR: ${HAPMAP3_RSIDS_FILE} has ${hapmap_rows} rows, expected ${HAPMAP3_EXPECTED_ROWS}" >&2
    exit 1
fi
if [[ "${hapmap_sha}" != "${HAPMAP3_EXPECTED_SHA256}" ]]; then
    echo "ERROR: ${HAPMAP3_RSIDS_FILE} SHA256 ${hapmap_sha}, expected ${HAPMAP3_EXPECTED_SHA256}" >&2
    exit 1
fi
if [[ "${hapmap_duplicates}" -ne 0 ]]; then
    echo "ERROR: ${HAPMAP3_RSIDS_FILE} has ${hapmap_duplicates} duplicate rsids" >&2
    exit 1
fi

fit_pca_samples=$(wc -l < "${FIT_PCA_KEEP}")
eur_samples=$(wc -l < "${EUR_KEEP}")
liftover_sha=$(sha256sum "${SBAYESRC_LIFTOVER_FILE}" | awk '{print $1}')
fit_pca_sha=$(sha256sum "${FIT_PCA_KEEP}" | awk '{print $1}')
eur_sha=$(sha256sum "${EUR_KEEP}" | awk '{print $1}')
filter_script_sha=$(sha256sum "${SCRIPT_DIR}/filter_hapmap3_wgs_hq_snps.py" | awk '{print $1}')

echo "AoU optional HapMap3 HQ bfile — $(date -u)"
echo "  HapMap3 rsids       = ${hapmap_rows}"
echo "  classified EUR      = ${eur_samples}"
echo "  fit-PCA samples     = ${fit_pca_samples}"
echo "  output              = ${DX_HAPMAP3_DIR}/hapmap3_bfile_hq.{bed,bim,fam}"
echo "  QC                  = AF diff<=${HAPMAP3_AF_DIFF_MAX}, MAF>=${HAPMAP3_MAF_MIN}, EUR missing<=${HAPMAP3_GENO_MAX}"

echo "  Validating existing WGS metric files ..."
for c in $(seq 1 22); do
    chrom="chr${c}"
    metric_summary="${DX_GWAS_METRICS_DIR}/wgs/${chrom}.metrics_summary.tsv"
    source_variants=$(metric_value "${metric_summary}" "source_variants")
    fit_lines=$(metric_value "${metric_summary}" "fit_pca_acount_lines")
    miss_lines=$(metric_value "${metric_summary}" "our_eur_vmiss_lines")
    metric_fit_samples=$(metric_value "${metric_summary}" "fit_pca_keep_samples")
    metric_eur_samples=$(metric_value "${metric_summary}" "eur_keep_samples")
    if ! is_uint "${source_variants}" || ! is_uint "${fit_lines}" || ! is_uint "${miss_lines}" ||
       ! is_uint "${metric_fit_samples}" || ! is_uint "${metric_eur_samples}"; then
        echo "ERROR: invalid metrics summary values in ${metric_summary}" >&2
        exit 1
    fi
    if [[ "${fit_lines}" -ne $((source_variants + 1)) || "${miss_lines}" -ne $((source_variants + 1)) ]]; then
        echo "ERROR: metric line counts do not match source variants for ${chrom}" >&2
        exit 1
    fi
    if [[ "${metric_fit_samples}" -ne "${fit_pca_samples}" || "${metric_eur_samples}" -ne "${eur_samples}" ]]; then
        echo "ERROR: metric sample counts do not match current keep lists for ${chrom}" >&2
        exit 1
    fi
done

desired_params="${LOCAL_HAPMAP3_DIR}/hapmap3_bfile_hq.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'hapmap3_rsids_file\t%s\n' "hapmap3_rsids.txt"
    printf 'hapmap3_rsids_rows\t%s\n' "${hapmap_rows}"
    printf 'hapmap3_rsids_sha256\t%s\n' "${hapmap_sha}"
    printf 'source_pfiles\t%s\n' "wgs_pfiles/chr1-22"
    printf 'metrics_source\t%s\n' "gwas_genotypes/metrics/wgs"
    printf 'af_diff_max_fit_pca_vs_snpinfo\t%s\n' "${HAPMAP3_AF_DIFF_MAX}"
    printf 'maf_min_fit_pca\t%s\n' "${HAPMAP3_MAF_MIN}"
    printf 'geno_max_classified_eur\t%s\n' "${HAPMAP3_GENO_MAX}"
    printf 'classified_european_samples\t%s\n' "${eur_samples}"
    printf 'classified_european_iids_sha256\t%s\n' "${eur_sha}"
    printf 'fit_pca_samples\t%s\n' "${fit_pca_samples}"
    printf 'fit_pca_iids_sha256\t%s\n' "${fit_pca_sha}"
    printf 'liftover_file_size\t%s\n' "$(stat -c%s "${SBAYESRC_LIFTOVER_FILE}")"
    printf 'liftover_sha256\t%s\n' "${liftover_sha}"
    printf 'filter_script_sha256\t%s\n' "${filter_script_sha}"
    printf 'plink2_version\t%s\n' "$("${PLINK2}" --version 2>&1 | head -1)"
    for c in $(seq 1 22); do
        chrom="chr${c}"
        metric_summary="${DX_GWAS_METRICS_DIR}/wgs/${chrom}.metrics_summary.tsv"
        printf '%s_source_variants\t%s\n' "${chrom}" "$(metric_value "${metric_summary}" "source_variants")"
        printf '%s_metrics_summary_sha256\t%s\n' "${chrom}" "$(sha256sum "${metric_summary}" | awk '{print $1}')"
    done
} > "${desired_params}"

final_prefix="${DX_HAPMAP3_DIR}/hapmap3_bfile_hq"
params="${DX_HAPMAP3_DIR}/hapmap3_bfile_hq.params.tsv"
summary="${DX_HAPMAP3_DIR}/hapmap3_bfile_hq.filter_summary.tsv"
final_ok=0
if [[ -s "${final_prefix}.bed" && -s "${final_prefix}.bim" && -s "${final_prefix}.fam" &&
      -s "${params}" && -s "${summary}" ]] &&
   diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
    expected_final=$(awk -F'\t' '$1 == "final_hapmap3_hq_snps" {print $2; exit}' "${summary}")
    observed_variants=$(wc -l < "${final_prefix}.bim")
    observed_samples=$(wc -l < "${final_prefix}.fam")
    if is_uint "${expected_final}" &&
       [[ "${observed_variants}" -eq "${expected_final}" && "${observed_samples}" -eq "${eur_samples}" ]]; then
        final_ok=1
    fi
fi
if [[ "${final_ok}" -eq 1 ]]; then
    echo "  HapMap3 HQ bfile already exists and matches current parameters — skipping"
    exit 0
fi

echo "  Staging plink2 and HapMap3 rsids ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
gcloud storage cp "${HAPMAP3_RSIDS_FILE}" "${DX_HAPMAP3_URI}/hapmap3_rsids.txt" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

echo "  Building HapMap3 HQ extract lists from existing WGS metrics ..."
filter_output_dir="${LOCAL_HAPMAP3_DIR}/filter_outputs"
rm -rf "${filter_output_dir}"
mkdir -p "${filter_output_dir}"
python3 "${SCRIPT_DIR}/filter_hapmap3_wgs_hq_snps.py" \
    --hapmap-rsids "${HAPMAP3_RSIDS_FILE}" \
    --metrics-dir "${DX_GWAS_METRICS_DIR}" \
    --liftover "${SBAYESRC_LIFTOVER_FILE}" \
    --output-dir "${filter_output_dir}" \
    --af-diff-max "${HAPMAP3_AF_DIFF_MAX}" \
    --maf-min "${HAPMAP3_MAF_MIN}" \
    --geno-max "${HAPMAP3_GENO_MAX}"
cp "${desired_params}" "${filter_output_dir}/hapmap3_bfile_hq.params.tsv"

filter_summary="${filter_output_dir}/hapmap3_bfile_hq.filter_summary.tsv"
per_chrom_summary="${filter_output_dir}/hapmap3_bfile_hq.per_chrom_summary.tsv"
expected_final=$(awk -F'\t' '$1 == "final_hapmap3_hq_snps" {print $2; exit}' "${filter_summary}")
if ! is_uint "${expected_final}" || [[ "${expected_final}" -le 0 ]]; then
    echo "ERROR: invalid final HapMap3 HQ variant count in ${filter_summary}" >&2
    exit 1
fi
echo "  HapMap3 HQ passing variants = ${expected_final}"

echo "  Uploading HapMap3 filter outputs ..."
gcloud storage cp \
    "${filter_output_dir}/hapmap3_bfile_hq.filter_summary.tsv" \
    "${filter_output_dir}/hapmap3_bfile_hq.per_chrom_summary.tsv" \
    "${filter_output_dir}/hapmap3_bfile_hq.params.tsv" \
    "${filter_output_dir}/hapmap3_bfile_hq.absent_from_wgs.txt" \
    "${filter_output_dir}/hapmap3_bfile_hq.variant_qc.tsv.gz" \
    "${DX_HAPMAP3_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
gcloud storage rsync -r "${filter_output_dir}/extracts" "${DX_HAPMAP3_URI}/extracts" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

extract_tasks="${SCRIPT_DIR}/logs/dsub_hapmap3_wgs_extract_$(date +%Y%m%d_%H%M%S).tsv"
{
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '--env CHROM' '--env EXPECTED_VARIANTS' '--env EXPECTED_SAMPLES' \
        '--input PLINK2' '--input PGEN' '--input PVAR' '--input PSAM' \
        '--input EXTRACT' '--input KEEP' '--output-recursive OUTDIR'
    for c in $(seq 1 22); do
        chrom="chr${c}"
        expected=$(awk -F'\t' -v chrom="${chrom}" 'NR > 1 && $1 == chrom {print $7; exit}' "${per_chrom_summary}")
        if ! is_uint "${expected}" || [[ "${expected}" -le 0 ]]; then
            echo "ERROR: invalid ${chrom} passing count in ${per_chrom_summary}: ${expected}" >&2
            exit 1
        fi
        p="${DX_HAPMAP3_DIR}/pfiles/${chrom}"
        submit=1
        if [[ -s "${p}.pgen" && -s "${p}.pvar" && -s "${p}.psam" ]]; then
            observed=$(grep -vc '^#' "${p}.pvar")
            observed_samples=$(grep -vc '^#' "${p}.psam")
            if [[ "${observed}" -eq "${expected}" && "${observed_samples}" -eq "${eur_samples}" ]]; then
                submit=0
            fi
        fi
        if [[ "${submit}" -eq 1 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${chrom}" "${expected}" "${eur_samples}" "${DSUB_PLINK2_GS}" \
                "${DX_WGS_PFILE_URI}/${chrom}.pgen" \
                "${DX_WGS_PFILE_URI}/${chrom}.pvar" \
                "${DX_WGS_PFILE_URI}/${chrom}.psam" \
                "${DX_HAPMAP3_URI}/extracts/${chrom}.extract.txt" \
                "${DX_EUROPEANS_URI}/classified_european_iids.txt" \
                "${DX_HAPMAP3_URI}/pfiles/"
        fi
    done
} > "${extract_tasks}"
extract_count=$(( $(wc -l < "${extract_tasks}") - 1 ))
if [[ "${extract_count}" -gt 0 ]]; then
    echo "  Submitting ${extract_count} HapMap3 chromosome extraction task(s) ..."
    extract_log="${extract_tasks%.tsv}.dsub.out"
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
        --name "sbayesrc-hapmap3-extract" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_hapmap3_wgs_extract_worker.sh" \
        --tasks "${extract_tasks}" \
        --min-cores "${HAPMAP3_EXTRACT_DSUB_MIN_CORES}" \
        --min-ram "${HAPMAP3_EXTRACT_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE}" \
        --disk-size "${HAPMAP3_EXTRACT_DSUB_DISK_SIZE}" \
        --disk-type "${HAPMAP3_EXTRACT_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${extract_log}"
    dsub_rc=${PIPESTATUS[0]}
    [[ "${dsub_rc}" -eq 0 ]] || exit "${dsub_rc}"
else
    echo "  HapMap3 chromosome pfiles already exist — skipping extraction"
fi

for c in $(seq 1 22); do
    chrom="chr${c}"
    expected=$(awk -F'\t' -v chrom="${chrom}" 'NR > 1 && $1 == chrom {print $7; exit}' "${per_chrom_summary}")
    p="${DX_HAPMAP3_DIR}/pfiles/${chrom}"
    observed=$(grep -vc '^#' "${p}.pvar")
    observed_samples=$(grep -vc '^#' "${p}.psam")
    if [[ "${observed}" -ne "${expected}" || "${observed_samples}" -ne "${eur_samples}" ]]; then
        echo "ERROR: ${chrom} HapMap3 pfile has ${observed}/${observed_samples}, expected ${expected}/${eur_samples}" >&2
        exit 1
    fi
done

merge_run=1
if [[ -s "${final_prefix}.bed" && -s "${final_prefix}.bim" && -s "${final_prefix}.fam" &&
      -s "${params}" ]] &&
   diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
    observed=$(wc -l < "${final_prefix}.bim")
    observed_samples=$(wc -l < "${final_prefix}.fam")
    if [[ "${observed}" -eq "${expected_final}" && "${observed_samples}" -eq "${eur_samples}" ]]; then
        merge_run=0
    fi
fi
if [[ "${merge_run}" -eq 1 ]]; then
    echo "  Submitting HapMap3 pfile merge / bfile conversion job ..."
    merge_log="${SCRIPT_DIR}/logs/dsub_hapmap3_bfile_merge_$(date +%Y%m%d_%H%M%S).dsub.out"
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
        --name "sbayesrc-hapmap3-merge" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_hapmap3_bfile_merge_worker.sh" \
        --env EXPECTED_VARIANTS="${expected_final}" \
        --env EXPECTED_SAMPLES="${eur_samples}" \
        --input PLINK2="${DSUB_PLINK2_GS}" \
        --input-recursive PFILES="${DX_HAPMAP3_URI}/pfiles/" \
        --input PARAMS="${DX_HAPMAP3_URI}/hapmap3_bfile_hq.params.tsv" \
        --input FILTER_SUMMARY="${DX_HAPMAP3_URI}/hapmap3_bfile_hq.filter_summary.tsv" \
        --output-recursive OUTDIR="${DX_HAPMAP3_URI}/" \
        --min-cores "${HAPMAP3_MERGE_DSUB_MIN_CORES}" \
        --min-ram "${HAPMAP3_MERGE_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE}" \
        --disk-size "${HAPMAP3_MERGE_DSUB_DISK_SIZE}" \
        --disk-type "${HAPMAP3_MERGE_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${merge_log}"
    dsub_rc=${PIPESTATUS[0]}
    [[ "${dsub_rc}" -eq 0 ]] || exit "${dsub_rc}"
else
    echo "  Final HapMap3 HQ bfile already exists — skipping merge"
fi

observed_final=$(wc -l < "${final_prefix}.bim")
observed_samples=$(wc -l < "${final_prefix}.fam")
if [[ "${observed_final}" -ne "${expected_final}" || "${observed_samples}" -ne "${eur_samples}" ]]; then
    echo "ERROR: final HapMap3 HQ bfile has ${observed_final}/${observed_samples}, expected ${expected_final}/${eur_samples}" >&2
    exit 1
fi

echo "  Done: ${final_prefix}.{bed,bim,fam} (${observed_final} variants, ${observed_samples} samples)"
echo "  Filter summary: ${summary}"
