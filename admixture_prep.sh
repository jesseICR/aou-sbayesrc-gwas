#!/bin/bash
# admixture_prep.sh - Prepare ADMIXTURE K=6 projection inputs on AoU.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_URI:?DX_HQ_DIRECT_BFILE_URI not set}"
: "${DX_ADMIXTURE_SCRAP_DIR:?DX_ADMIXTURE_SCRAP_DIR not set}"
: "${DX_ADMIXTURE_SCRAP_URI:?DX_ADMIXTURE_SCRAP_URI not set}"
: "${LOCAL_ADMIXTURE_DIR:?LOCAL_ADMIXTURE_DIR not set}"
: "${ADMIXTURE_TSV_URL:?ADMIXTURE_TSV_URL not set}"
: "${ADMIXTURE_DOWNLOAD_URL:?ADMIXTURE_DOWNLOAD_URL not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"
: "${ADMIXTURE_GENO_MAX:?ADMIXTURE_GENO_MAX not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

mkdir -p "${LOCAL_ADMIXTURE_DIR}" "${SCRIPT_DIR}/tools" "${DX_ADMIXTURE_SCRAP_DIR}"

ADMIXTURE_PREP_DSUB_MIN_CORES="${ADMIXTURE_PREP_DSUB_MIN_CORES:-8}"
ADMIXTURE_PREP_DSUB_MIN_RAM="${ADMIXTURE_PREP_DSUB_MIN_RAM:-32}"
ADMIXTURE_PREP_DSUB_DISK_SIZE="${ADMIXTURE_PREP_DSUB_DISK_SIZE:-300}"
ADMIXTURE_PREP_DSUB_DISK_TYPE="${ADMIXTURE_PREP_DSUB_DISK_TYPE:-pd-ssd}"

source_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
source_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"
for ext in bed bim fam; do
    if [[ ! -s "${source_prefix}.${ext}" ]]; then
        echo "ERROR: missing HQ direct bfile input ${source_prefix}.${ext}" >&2
        exit 1
    fi
done

ref_tsv="${LOCAL_ADMIXTURE_DIR}/admixture_allele_freqs.tsv"
if [[ -s "${ref_tsv}" ]]; then
    echo "  ADMIXTURE reference TSV already cached (${ref_tsv}, $(wc -l < "${ref_tsv}") lines)"
else
    echo "  Downloading ADMIXTURE K=${ADMIXTURE_K} reference TSV ..."
    curl -fsSL --retry 3 --retry-delay 5 -o "${ref_tsv}" "${ADMIXTURE_TSV_URL}"
    echo "  Downloaded ${ref_tsv} ($(wc -l < "${ref_tsv}") lines)"
fi

admixture_bin="${SCRIPT_DIR}/tools/admixture"
if [[ -x "${admixture_bin}" ]]; then
    echo "  ADMIXTURE binary already cached at ${admixture_bin}"
else
    echo "  Downloading ADMIXTURE binary ..."
    tarball="${SCRIPT_DIR}/tools/admixture_linux.tar.gz"
    curl -fsSL --retry 3 --retry-delay 5 -o "${tarball}" "${ADMIXTURE_DOWNLOAD_URL}"
    tar xzf "${tarball}" -C "${SCRIPT_DIR}/tools"
    found_bin=$(find "${SCRIPT_DIR}/tools" -path '*/admixture' -type f | head -1)
    if [[ -z "${found_bin}" ]]; then
        echo "ERROR: ADMIXTURE binary was not found after unpacking ${tarball}" >&2
        exit 1
    fi
    mv "${found_bin}" "${admixture_bin}"
    chmod +x "${admixture_bin}"
    rm -rf "${SCRIPT_DIR}/tools/dist" "${tarball}"
    echo "  Downloaded ADMIXTURE to ${admixture_bin}"
fi

source_variants=$(wc -l < "${source_prefix}.bim")
source_samples=$(wc -l < "${source_prefix}.fam")
ref_variants=$(awk 'NR > 1 {n++} END {print n + 0}' "${ref_tsv}")
ref_sha256=$(sha256sum "${ref_tsv}" | awk '{print $1}')
admixture_sha256=$(sha256sum "${admixture_bin}" | awk '{print $1}')

desired_params="${LOCAL_ADMIXTURE_DIR}/admixture_prep.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'source_bfile\t%s\n' 'direct_bfile_hq/chr1_22_merged_hq'
    printf 'source_variants\t%s\n' "${source_variants}"
    printf 'source_samples\t%s\n' "${source_samples}"
    printf 'admixture_k\t%s\n' "${ADMIXTURE_K}"
    printf 'admixture_geno_max\t%s\n' "${ADMIXTURE_GENO_MAX}"
    printf 'admixture_tsv_url\t%s\n' "${ADMIXTURE_TSV_URL}"
    printf 'admixture_tsv_sha256\t%s\n' "${ref_sha256}"
    printf 'admixture_reference_variants\t%s\n' "${ref_variants}"
    printf 'admixture_binary_sha256\t%s\n' "${admixture_sha256}"
} > "${desired_params}"

summary="${DX_ADMIXTURE_SCRAP_DIR}/admixture_prep_summary.tsv"
params="${DX_ADMIXTURE_SCRAP_DIR}/admixture_prep.params.tsv"
aligned_prefix="${DX_ADMIXTURE_SCRAP_DIR}/aou_admixture_aligned"
p_file="${DX_ADMIXTURE_SCRAP_DIR}/ref_aligned.P"
if [[ -s "${aligned_prefix}.bed" && -s "${aligned_prefix}.bim" && -s "${aligned_prefix}.fam" && -s "${p_file}" && -s "${summary}" && -s "${params}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        final_variants=$(awk -F'\t' '$1 == "final_aligned_variants" {print $2; exit}' "${summary}")
        observed_variants=$(wc -l < "${aligned_prefix}.bim")
        p_rows=$(wc -l < "${p_file}")
        observed_samples=$(wc -l < "${aligned_prefix}.fam")
        if [[ -n "${final_variants}" && "${observed_variants}" -eq "${final_variants}" && "${p_rows}" -eq "${final_variants}" && "${observed_samples}" -eq "${source_samples}" ]]; then
            echo "  ADMIXTURE prep already complete (${observed_variants} variants, ${observed_samples} samples) — skipping"
            exit 0
        fi
    fi
    echo "  ADMIXTURE prep outputs exist but params/counts do not match — rebuilding"
fi

echo "  Staging ADMIXTURE prep inputs to workspace bucket ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
export DSUB_ADMIXTURE_PREP_URI="${DSUB_ADMIXTURE_PREP_URI:-${DX_ADMIXTURE_SCRAP_URI}/_prep}"
gcloud storage cp \
    "${ref_tsv}" \
    "${admixture_bin}" \
    "${desired_params}" \
    "${DSUB_ADMIXTURE_PREP_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

prep_log="${SCRIPT_DIR}/logs/dsub_admixture_prep_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${prep_log}")"

echo "  Submitting dsub job for ADMIXTURE prep/alignment ..."
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
    --name "sbayesrc-admixture-prep" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_admixture_prep_worker.sh" \
    --env ADMIXTURE_GENO_MAX="${ADMIXTURE_GENO_MAX}" \
    --env EXPECTED_SOURCE_VARIANTS="${source_variants}" \
    --env EXPECTED_SOURCE_SAMPLES="${source_samples}" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input BED="${source_uri_prefix}.bed" \
    --input BIM="${source_uri_prefix}.bim" \
    --input FAM="${source_uri_prefix}.fam" \
    --input REF_TSV="${DSUB_ADMIXTURE_PREP_URI}/admixture_allele_freqs.tsv" \
    --input PARAMS="${DSUB_ADMIXTURE_PREP_URI}/admixture_prep.desired_params.tsv" \
    --output-recursive OUTDIR="${DX_ADMIXTURE_SCRAP_URI}/" \
    --min-cores "${ADMIXTURE_PREP_DSUB_MIN_CORES}" \
    --min-ram "${ADMIXTURE_PREP_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${ADMIXTURE_PREP_DSUB_DISK_SIZE}" \
    --disk-type "${ADMIXTURE_PREP_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${prep_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub ADMIXTURE prep job returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

echo "  Verifying ADMIXTURE prep outputs ..."
for attempt in $(seq 1 60); do
    if [[ -s "${aligned_prefix}.bed" && -s "${aligned_prefix}.bim" && -s "${aligned_prefix}.fam" && -s "${p_file}" && -s "${summary}" ]]; then
        observed_variants=$(wc -l < "${aligned_prefix}.bim")
        observed_samples=$(wc -l < "${aligned_prefix}.fam")
        p_rows=$(wc -l < "${p_file}")
        if [[ "${observed_variants}" -eq "${p_rows}" && "${observed_samples}" -eq "${source_samples}" ]]; then
            echo "  Done: ${aligned_prefix}.{bed,bim,fam} (${observed_variants} variants, ${observed_samples} samples)"
            echo "  Summary: ${summary}"
            exit 0
        fi
        echo "  Outputs visible but counts do not match yet; waiting ..."
    else
        echo "  Waiting for ADMIXTURE prep outputs to appear (${attempt}/60) ..."
    fi
    sleep 10
done

echo "ERROR: ADMIXTURE prep outputs were not visible with expected counts after dsub completed" >&2
exit 1
