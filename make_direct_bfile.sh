#!/bin/bash
# make_direct_bfile.sh — Merge per-chromosome present direct-SNP pfiles into one bfile.
#
# Output:
#   ${DX_DIRECT_BFILE_DIR}/chr1_22_merged.{bed,bim,fam}

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_DIRECT_PFILE_DIR:?DX_DIRECT_PFILE_DIR not set}"
: "${DX_DIRECT_PFILE_URI:?DX_DIRECT_PFILE_URI not set}"
: "${DX_DIRECT_BFILE_DIR:?DX_DIRECT_BFILE_DIR not set}"
: "${DX_DIRECT_BFILE_URI:?DX_DIRECT_BFILE_URI not set}"
: "${LOCAL_DIRECT_PREP_DIR:?LOCAL_DIRECT_PREP_DIR not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

mkdir -p "${DX_DIRECT_BFILE_DIR}"

DIRECT_BFILE_DSUB_MIN_CORES="${DIRECT_BFILE_DSUB_MIN_CORES:-8}"
DIRECT_BFILE_DSUB_MIN_RAM="${DIRECT_BFILE_DSUB_MIN_RAM:-32}"
DIRECT_BFILE_DSUB_DISK_SIZE="${DIRECT_BFILE_DSUB_DISK_SIZE:-300}"
DIRECT_BFILE_DSUB_DISK_TYPE="${DIRECT_BFILE_DSUB_DISK_TYPE:-pd-ssd}"

summary="${LOCAL_DIRECT_PREP_DIR}/summary.tsv"
if [[ ! -s "${summary}" ]]; then
    echo "ERROR: missing direct SNP summary: ${summary}" >&2
    echo "Run prepare_direct_snps.py before make_direct_bfile.sh." >&2
    exit 1
fi

requested_direct=$(awk -F'\t' 'NR > 1 { s += $2 } END { print s + 0 }' "${summary}")
expected_present=$(awk -F'\t' 'NR > 1 { s += $3 } END { print s + 0 }' "${summary}")
expected_missing=$(awk -F'\t' 'NR > 1 { s += $4 } END { print s + 0 }' "${summary}")
if [[ "${requested_direct}" -le 0 || "${expected_present}" -le 0 ]]; then
    echo "ERROR: could not determine expected direct SNP count from ${summary}" >&2
    exit 1
fi

out_prefix="${DX_DIRECT_BFILE_DIR}/chr1_22_merged"
if [[ -s "${out_prefix}.bed" && -s "${out_prefix}.bim" && -s "${out_prefix}.fam" ]]; then
    existing=$(wc -l < "${out_prefix}.bim")
    if [[ "${existing}" -eq "${expected_present}" ]]; then
        echo "  Direct bfile already exists (${existing}/${expected_present}) — skipping"
        exit 0
    fi
    echo "  Direct bfile exists but has ${existing}/${expected_present} variants — rebuilding"
fi

for chrom in $(seq 1 22); do
    prefix="${DX_DIRECT_PFILE_DIR}/chr${chrom}"
    for ext in pgen pvar psam; do
        if [[ ! -s "${prefix}.${ext}" ]]; then
            echo "ERROR: missing direct pfile input ${prefix}.${ext}" >&2
            exit 1
        fi
    done
done

echo "  Staging plink2 for dsub direct-bfile merge ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

export DSUB_DIRECT_PREP_URI="${DSUB_DIRECT_PREP_URI:-${DX_DIRECT_PFILE_URI}/_prep}"
tasks_log="${SCRIPT_DIR}/logs/dsub_direct_bfile_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${tasks_log}")"

echo "  Submitting dsub direct-bfile merge in ${DSUB_REGION} ..."
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
    --name "sbayesrc-direct-bfile" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_direct_bfile_worker.sh" \
    --env OUTPUT_PREFIX="${DX_DIRECT_BFILE_URI}/chr1_22_merged" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input SUMMARY="${DSUB_DIRECT_PREP_URI}/summary.tsv" \
    --input-recursive DIRECT_PFILES="${DX_DIRECT_PFILE_URI}/" \
    --output-recursive OUTDIR="${DX_DIRECT_BFILE_URI}/" \
    --min-cores "${DIRECT_BFILE_DSUB_MIN_CORES}" \
    --min-ram "${DIRECT_BFILE_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${DIRECT_BFILE_DSUB_DISK_SIZE}" \
    --disk-type "${DIRECT_BFILE_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${tasks_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub direct-bfile merge returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

echo "  Verifying direct bfile on bucket ..."
for attempt in $(seq 1 60); do
    if [[ -s "${out_prefix}.bed" && -s "${out_prefix}.bim" && -s "${out_prefix}.fam" ]]; then
        observed=$(wc -l < "${out_prefix}.bim")
        if [[ "${observed}" -eq "${expected_present}" ]]; then
            samples=$(wc -l < "${out_prefix}.fam")
            echo "  Done: ${out_prefix}.{bed,bim,fam} (${observed} variants, ${samples} samples)"
            echo "  Requested direct SNPs: ${requested_direct}; absent from WGS pfiles: ${expected_missing}"
            exit 0
        fi
        echo "  Output visible but has ${observed}/${expected_present} variants; waiting for consistency ..."
    else
        echo "  Waiting for direct bfile outputs to appear (${attempt}/60) ..."
    fi
    sleep 10
done

echo "ERROR: direct bfile outputs were not visible with expected counts after dsub completed" >&2
exit 1
