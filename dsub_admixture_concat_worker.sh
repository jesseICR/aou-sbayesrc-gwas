#!/bin/bash
# dsub_admixture_concat_worker.sh - Concatenate ADMIXTURE batch .Q files.

set -euo pipefail

: "${FAM:?FAM not set}"
: "${MANIFEST:?MANIFEST not set}"
: "${QDIR:?QDIR not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"
: "${EXPECTED_SAMPLES:?EXPECTED_SAMPLES not set}"
: "${EXPECTED_BATCHES:?EXPECTED_BATCHES not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[admixture-concat $(ts)] $*"; }

mkdir -p "${OUTDIR}"

log "=== starting concat on $(hostname) ==="
fam_samples=$(wc -l < "${FAM}")
if [[ "${fam_samples}" -ne "${EXPECTED_SAMPLES}" ]]; then
    log "ERROR: FAM samples ${fam_samples} != expected ${EXPECTED_SAMPLES}"
    exit 1
fi

q_all="${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.Q"
: > "${q_all}"
observed_batches=0
while IFS=$'\t' read -r batch start end samples variants p_rows; do
    [[ "${batch}" == "batch" ]] && continue
    q_file="${QDIR}/${batch}.${ADMIXTURE_K}.Q"
    if [[ ! -s "${q_file}" ]]; then
        log "ERROR: missing ${q_file}"
        exit 1
    fi
    q_lines=$(wc -l < "${q_file}")
    if [[ "${q_lines}" -ne "${samples}" ]]; then
        log "ERROR: ${q_file} has ${q_lines} lines, expected ${samples}"
        exit 1
    fi
    cat "${q_file}" >> "${q_all}"
    observed_batches=$((observed_batches + 1))
done < "${MANIFEST}"

if [[ "${observed_batches}" -ne "${EXPECTED_BATCHES}" ]]; then
    log "ERROR: observed ${observed_batches} batches, expected ${EXPECTED_BATCHES}"
    exit 1
fi

q_total=$(wc -l < "${q_all}")
if [[ "${q_total}" -ne "${fam_samples}" ]]; then
    log "ERROR: concatenated Q rows ${q_total} != FAM samples ${fam_samples}"
    exit 1
fi

awk '{print $1 "\t" $2}' "${FAM}" > "${OUTDIR}/aou_admixture_iids.tsv"
tr ' ' '\t' < "${q_all}" > "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.Q.tsv"
{
    printf 'FID\tIID\tEuropean\tEast_Asian\tAmerican\tAfrican\tSouth_Asian\tOceanian\n'
    paste "${OUTDIR}/aou_admixture_iids.tsv" "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.Q.tsv"
} > "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.tsv"

cp "${PARAMS}" "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'samples\t%s\n' "${fam_samples}"
    printf 'batches\t%s\n' "${observed_batches}"
    printf 'admixture_k\t%s\n' "${ADMIXTURE_K}"
    printf 'q_rows\t%s\n' "${q_total}"
    printf 'columns\tFID,IID,European,East_Asian,American,African,South_Asian,Oceanian\n'
} > "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.summary.tsv"

rm -f "${OUTDIR}/aou_admixture_iids.tsv" "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.Q.tsv"

log "=== done: ${q_total} samples in aou_admixture_k${ADMIXTURE_K}.tsv ==="
cat "${OUTDIR}/aou_admixture_k${ADMIXTURE_K}.summary.tsv" | sed 's/^/  /'
