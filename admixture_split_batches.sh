#!/bin/bash
# admixture_split_batches.sh - Split ADMIXTURE aligned bfile into 20k batches.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_ADMIXTURE_SCRAP_DIR:?DX_ADMIXTURE_SCRAP_DIR not set}"
: "${DX_ADMIXTURE_SCRAP_URI:?DX_ADMIXTURE_SCRAP_URI not set}"
: "${DX_ADMIXTURE_BATCH_DIR:?DX_ADMIXTURE_BATCH_DIR not set}"
: "${DX_ADMIXTURE_BATCH_URI:?DX_ADMIXTURE_BATCH_URI not set}"
: "${LOCAL_ADMIXTURE_DIR:?LOCAL_ADMIXTURE_DIR not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"
: "${ADMIXTURE_BATCH_SIZE:?ADMIXTURE_BATCH_SIZE not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

mkdir -p "${LOCAL_ADMIXTURE_DIR}" "${DX_ADMIXTURE_BATCH_DIR}"

ADMIXTURE_SPLIT_DSUB_MIN_CORES="${ADMIXTURE_SPLIT_DSUB_MIN_CORES:-8}"
ADMIXTURE_SPLIT_DSUB_MIN_RAM="${ADMIXTURE_SPLIT_DSUB_MIN_RAM:-32}"
ADMIXTURE_SPLIT_DSUB_DISK_SIZE="${ADMIXTURE_SPLIT_DSUB_DISK_SIZE:-300}"
ADMIXTURE_SPLIT_DSUB_DISK_TYPE="${ADMIXTURE_SPLIT_DSUB_DISK_TYPE:-pd-ssd}"

aligned_prefix="${DX_ADMIXTURE_SCRAP_DIR}/aou_admixture_aligned"
aligned_uri_prefix="${DX_ADMIXTURE_SCRAP_URI}/aou_admixture_aligned"
p_file="${DX_ADMIXTURE_SCRAP_DIR}/ref_aligned.P"
p_uri="${DX_ADMIXTURE_SCRAP_URI}/ref_aligned.P"
prep_params="${DX_ADMIXTURE_SCRAP_DIR}/admixture_prep.params.tsv"
prep_params_uri="${DX_ADMIXTURE_SCRAP_URI}/admixture_prep.params.tsv"

for ext in bed bim fam; do
    if [[ ! -s "${aligned_prefix}.${ext}" ]]; then
        echo "ERROR: missing ADMIXTURE aligned bfile ${aligned_prefix}.${ext}" >&2
        exit 1
    fi
done
if [[ ! -s "${p_file}" || ! -s "${prep_params}" ]]; then
    echo "ERROR: missing ADMIXTURE prep outputs under ${DX_ADMIXTURE_SCRAP_DIR}" >&2
    exit 1
fi

aligned_variants=$(wc -l < "${aligned_prefix}.bim")
aligned_samples=$(wc -l < "${aligned_prefix}.fam")
p_rows=$(wc -l < "${p_file}")
if [[ "${aligned_variants}" -ne "${p_rows}" ]]; then
    echo "ERROR: aligned BIM rows (${aligned_variants}) != ref_aligned.P rows (${p_rows})" >&2
    exit 1
fi
n_batches=$(( (aligned_samples + ADMIXTURE_BATCH_SIZE - 1) / ADMIXTURE_BATCH_SIZE ))
prep_sha256=$(sha256sum "${prep_params}" | awk '{print $1}')

desired_params="${LOCAL_ADMIXTURE_DIR}/admixture_split.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'admixture_k\t%s\n' "${ADMIXTURE_K}"
    printf 'batch_size\t%s\n' "${ADMIXTURE_BATCH_SIZE}"
    printf 'aligned_variants\t%s\n' "${aligned_variants}"
    printf 'aligned_samples\t%s\n' "${aligned_samples}"
    printf 'ref_aligned_p_rows\t%s\n' "${p_rows}"
    printf 'prep_params_sha256\t%s\n' "${prep_sha256}"
} > "${desired_params}"

manifest="${DX_ADMIXTURE_BATCH_DIR}/batch_manifest.tsv"
params="${DX_ADMIXTURE_BATCH_DIR}/admixture_split.params.tsv"
summary="${DX_ADMIXTURE_BATCH_DIR}/admixture_split_summary.tsv"

verify_batches() {
    [[ -s "${manifest}" && -s "${params}" && -s "${summary}" ]] || return 1
    diff -q "${desired_params}" "${params}" >/dev/null 2>&1 || return 1
    local rows
    rows=$(awk 'NR > 1 {n++} END {print n + 0}' "${manifest}")
    [[ "${rows}" -eq "${n_batches}" ]] || return 1
    while IFS=$'\t' read -r batch start end samples variants rows_p; do
        [[ "${batch}" == "batch" ]] && continue
        [[ -s "${DX_ADMIXTURE_BATCH_DIR}/${batch}.bed" ]] || return 1
        [[ -s "${DX_ADMIXTURE_BATCH_DIR}/${batch}.bim" ]] || return 1
        [[ -s "${DX_ADMIXTURE_BATCH_DIR}/${batch}.fam" ]] || return 1
        [[ -s "${DX_ADMIXTURE_BATCH_DIR}/${batch}.${ADMIXTURE_K}.P.in" ]] || return 1
        [[ "$(wc -l < "${DX_ADMIXTURE_BATCH_DIR}/${batch}.fam")" -eq "${samples}" ]] || return 1
        [[ "$(wc -l < "${DX_ADMIXTURE_BATCH_DIR}/${batch}.bim")" -eq "${variants}" ]] || return 1
        [[ "$(wc -l < "${DX_ADMIXTURE_BATCH_DIR}/${batch}.${ADMIXTURE_K}.P.in")" -eq "${rows_p}" ]] || return 1
    done < "${manifest}"
}

if verify_batches; then
    echo "  ADMIXTURE batches already exist (${n_batches} batches, ${aligned_samples} samples) — skipping"
    exit 0
fi

echo "  Staging plink2 + split params to workspace bucket ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
export DSUB_ADMIXTURE_PREP_URI="${DSUB_ADMIXTURE_PREP_URI:-${DX_ADMIXTURE_SCRAP_URI}/_prep}"
gcloud storage cp "${desired_params}" "${DSUB_ADMIXTURE_PREP_URI}/admixture_split.desired_params.tsv" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

split_log="${SCRIPT_DIR}/logs/dsub_admixture_split_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${split_log}")"

echo "  Submitting dsub job to split ADMIXTURE bfile into ${n_batches} batches ..."
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
    --name "sbayesrc-admixture-split" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_admixture_split_worker.sh" \
    --env ADMIXTURE_K="${ADMIXTURE_K}" \
    --env ADMIXTURE_BATCH_SIZE="${ADMIXTURE_BATCH_SIZE}" \
    --env EXPECTED_VARIANTS="${aligned_variants}" \
    --env EXPECTED_SAMPLES="${aligned_samples}" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input BED="${aligned_uri_prefix}.bed" \
    --input BIM="${aligned_uri_prefix}.bim" \
    --input FAM="${aligned_uri_prefix}.fam" \
    --input PFILE="${p_uri}" \
    --input PARAMS="${DSUB_ADMIXTURE_PREP_URI}/admixture_split.desired_params.tsv" \
    --output-recursive OUTDIR="${DX_ADMIXTURE_BATCH_URI}/" \
    --min-cores "${ADMIXTURE_SPLIT_DSUB_MIN_CORES}" \
    --min-ram "${ADMIXTURE_SPLIT_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${ADMIXTURE_SPLIT_DSUB_DISK_SIZE}" \
    --disk-type "${ADMIXTURE_SPLIT_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${split_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub ADMIXTURE split job returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

echo "  Verifying ADMIXTURE batches ..."
for attempt in $(seq 1 60); do
    if verify_batches; then
        echo "  Done: ${n_batches} ADMIXTURE batches in ${DX_ADMIXTURE_BATCH_DIR}"
        exit 0
    fi
    echo "  Waiting for ADMIXTURE batch outputs to appear (${attempt}/60) ..."
    sleep 10
done

echo "ERROR: ADMIXTURE batch outputs were not visible with expected counts after dsub completed" >&2
exit 1
