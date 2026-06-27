#!/bin/bash
# dsub_hapmap3_wgs_extract_worker.sh - Extract one chromosome of HapMap3 HQ WGS pfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${PGEN:?PGEN not set}"
: "${PVAR:?PVAR not set}"
: "${PSAM:?PSAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${KEEP:?KEEP not set}"
: "${CHROM:?CHROM not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"
: "${EXPECTED_SAMPLES:?EXPECTED_SAMPLES not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[hapmap3-extract ${CHROM} $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH="/mnt/data/scratch/hapmap3_extract_${CHROM}"
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${PGEN}" "${SCRATCH}/${CHROM}.pgen"
ln -sf "${PVAR}" "${SCRATCH}/${CHROM}.pvar"
ln -sf "${PSAM}" "${SCRATCH}/${CHROM}.psam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "extract variants = $(wc -l < "${EXTRACT}")"
log "keep samples = $(wc -l < "${KEEP}")"
df -h /mnt/data | sed 's/^/  /'

"${PLINK2}" \
    --pfile "${SCRATCH}/${CHROM}" \
    --extract "${EXTRACT}" \
    --keep "${KEEP}" \
    --no-pheno \
    --output-chr chrM \
    --make-pgen \
    --sort-vars \
    --indiv-sort none \
    --threads "$(nproc)" \
    --out "${OUTDIR}/${CHROM}"

observed_variants=$(grep -vc '^#' "${OUTDIR}/${CHROM}.pvar")
observed_samples=$(grep -vc '^#' "${OUTDIR}/${CHROM}.psam")
if [[ "${observed_variants}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: ${CHROM} pfile has ${observed_variants} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi
if [[ "${observed_samples}" -ne "${EXPECTED_SAMPLES}" ]]; then
    log "ERROR: ${CHROM} pfile has ${observed_samples} samples, expected ${EXPECTED_SAMPLES}"
    exit 1
fi

{
    printf 'metric\tvalue\n'
    printf 'chrom\t%s\n' "${CHROM}"
    printf 'final_pfile\t%s\n' "hapmap3_bfile_hq/pfiles/${CHROM}"
    printf 'final_variants\t%s\n' "${observed_variants}"
    printf 'final_samples\t%s\n' "${observed_samples}"
} > "${OUTDIR}/${CHROM}.summary.tsv"

log "=== done: ${observed_variants} variants, ${observed_samples} samples ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
