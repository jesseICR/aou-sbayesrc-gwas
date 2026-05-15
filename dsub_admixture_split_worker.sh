#!/bin/bash
# dsub_admixture_split_worker.sh - Split aligned ADMIXTURE bfile into batches.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${PFILE:?PFILE not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"
: "${ADMIXTURE_BATCH_SIZE:?ADMIXTURE_BATCH_SIZE not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"
: "${EXPECTED_SAMPLES:?EXPECTED_SAMPLES not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[admixture-split $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/admixture_split
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/aligned.bed"
ln -sf "${BIM}" "${SCRATCH}/aligned.bim"
ln -sf "${FAM}" "${SCRATCH}/aligned.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "ADMIXTURE_K=${ADMIXTURE_K} batch_size=${ADMIXTURE_BATCH_SIZE}"
df -h /mnt/data | sed 's/^/  /'

variants=$(wc -l < "${SCRATCH}/aligned.bim")
samples=$(wc -l < "${SCRATCH}/aligned.fam")
p_rows=$(wc -l < "${PFILE}")
if [[ "${variants}" -ne "${EXPECTED_VARIANTS}" || "${p_rows}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: variant/P count mismatch: variants=${variants}, p_rows=${p_rows}, expected=${EXPECTED_VARIANTS}"
    exit 1
fi
if [[ "${samples}" -ne "${EXPECTED_SAMPLES}" ]]; then
    log "ERROR: sample count mismatch: samples=${samples}, expected=${EXPECTED_SAMPLES}"
    exit 1
fi

n_batches=$(( (samples + ADMIXTURE_BATCH_SIZE - 1) / ADMIXTURE_BATCH_SIZE ))
manifest="${OUTDIR}/batch_manifest.tsv"
printf 'batch\tstart_line\tend_line\tsamples\tvariants\tp_rows\n' > "${manifest}"

for i in $(seq 1 "${n_batches}"); do
    batch=$(printf 'batch_%03d' "${i}")
    start=$(( (i - 1) * ADMIXTURE_BATCH_SIZE + 1 ))
    end=$(( i * ADMIXTURE_BATCH_SIZE ))
    if [[ "${end}" -gt "${samples}" ]]; then
        end="${samples}"
    fi

    awk -v s="${start}" -v e="${end}" 'NR >= s && NR <= e {print $1, $2}' \
        "${SCRATCH}/aligned.fam" > "${SCRATCH}/${batch}.keep"

    log "building ${batch} (${start}-${end})"
    "${PLINK2}" \
        --bfile "${SCRATCH}/aligned" \
        --keep "${SCRATCH}/${batch}.keep" \
        --make-bed \
        --threads "$(nproc)" \
        --out "${OUTDIR}/${batch}"

    batch_samples=$(wc -l < "${OUTDIR}/${batch}.fam")
    batch_variants=$(wc -l < "${OUTDIR}/${batch}.bim")
    if [[ "${batch_samples}" -ne $((end - start + 1)) || "${batch_variants}" -ne "${variants}" ]]; then
        log "ERROR: ${batch} count mismatch: samples=${batch_samples}, variants=${batch_variants}"
        exit 1
    fi

    cp "${PFILE}" "${OUTDIR}/${batch}.${ADMIXTURE_K}.P.in"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${batch}" "${start}" "${end}" "${batch_samples}" "${batch_variants}" "${p_rows}" >> "${manifest}"
    rm -f "${SCRATCH}/${batch}.keep"
done

cp "${PARAMS}" "${OUTDIR}/admixture_split.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'samples\t%s\n' "${samples}"
    printf 'variants\t%s\n' "${variants}"
    printf 'batch_size\t%s\n' "${ADMIXTURE_BATCH_SIZE}"
    printf 'batches\t%s\n' "${n_batches}"
    printf 'admixture_k\t%s\n' "${ADMIXTURE_K}"
} > "${OUTDIR}/admixture_split_summary.tsv"

log "=== done: ${n_batches} batches, ${samples} samples, ${variants} variants ==="
cat "${OUTDIR}/admixture_split_summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
