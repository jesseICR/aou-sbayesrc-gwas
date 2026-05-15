#!/bin/bash
# dsub_admixture_project_worker.sh - Run ADMIXTURE projection for one batch.

set -euo pipefail

: "${ADMIXTURE:?ADMIXTURE not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${PFILE:?PFILE not set}"
: "${Q:?Q not set}"
: "${BATCH:?BATCH not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"
: "${EXPECTED_SAMPLES:?EXPECTED_SAMPLES not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[admixture-project $(ts)] $*"; }

chmod +x "${ADMIXTURE}"

SCRATCH=/mnt/data/scratch/admixture_project
mkdir -p "${SCRATCH}"
ln -sf "${BED}" "${SCRATCH}/${BATCH}.bed"
ln -sf "${BIM}" "${SCRATCH}/${BATCH}.bim"
ln -sf "${FAM}" "${SCRATCH}/${BATCH}.fam"
ln -sf "${PFILE}" "${SCRATCH}/${BATCH}.${ADMIXTURE_K}.P.in"

log "=== starting ${BATCH} on $(hostname) ==="
log "ADMIXTURE_K=${ADMIXTURE_K}"
admixture_threads="${ADMIXTURE_THREADS:-$(nproc)}"
log "ADMIXTURE_THREADS=${admixture_threads}"
df -h /mnt/data | sed 's/^/  /'

samples=$(wc -l < "${SCRATCH}/${BATCH}.fam")
variants=$(wc -l < "${SCRATCH}/${BATCH}.bim")
if [[ "${samples}" -ne "${EXPECTED_SAMPLES}" || "${variants}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: input counts for ${BATCH}: samples=${samples}/${EXPECTED_SAMPLES}, variants=${variants}/${EXPECTED_VARIANTS}"
    exit 1
fi

(
    cd "${SCRATCH}"
    "${ADMIXTURE}" -j"${admixture_threads}" -P "${BATCH}.bed" "${ADMIXTURE_K}" | tee "${BATCH}.${ADMIXTURE_K}.admixture.log"
)

q_file="${SCRATCH}/${BATCH}.${ADMIXTURE_K}.Q"
if [[ ! -s "${q_file}" ]]; then
    log "ERROR: missing ADMIXTURE Q output ${q_file}"
    exit 1
fi
q_lines=$(wc -l < "${q_file}")
if [[ "${q_lines}" -ne "${samples}" ]]; then
    log "ERROR: ${BATCH} Q lines ${q_lines} != samples ${samples}"
    exit 1
fi

cp "${q_file}" "${Q}"
log "=== done ${BATCH}: ${q_lines} Q rows ==="
