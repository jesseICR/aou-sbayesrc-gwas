#!/bin/bash
# dsub_gwas_step2_wgs_worker.sh - Build final REGENIE Step 2 WGS pfile for one chromosome.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${PGEN:?PGEN not set}"
: "${PVAR:?PVAR not set}"
: "${PSAM:?PSAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${PARAMS:?PARAMS not set}"
: "${CHROM:?CHROM not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[gwas-step2 ${CHROM} $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH="/mnt/data/scratch/gwas_step2_${CHROM}"
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${PGEN}" "${SCRATCH}/${CHROM}.pgen"
ln -sf "${PVAR}" "${SCRATCH}/${CHROM}.pvar"
ln -sf "${PSAM}" "${SCRATCH}/${CHROM}.psam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "extract variants = $(wc -l < "${EXTRACT}")"
df -h /mnt/data | sed 's/^/  /'

"${PLINK2}" \
    --pfile "${SCRATCH}/${CHROM}" \
    --extract "${EXTRACT}" \
    --no-pheno \
    --output-chr chrM \
    --make-pgen \
    --threads "$(nproc)" \
    --out "${OUTDIR}/${CHROM}"

observed=$(grep -vc '^#' "${OUTDIR}/${CHROM}.pvar")
samples=$(grep -vc '^#' "${OUTDIR}/${CHROM}.psam")
if [[ "${observed}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: final ${CHROM} pfile has ${observed} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi

cp "${PARAMS}" "${OUTDIR}/${CHROM}.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'chrom\t%s\n' "${CHROM}"
    printf 'source_pfile\t%s\n' "wgs_pfiles/${CHROM}"
    printf 'final_pfile\t%s\n' "gwas_genotypes/step2_wgs_pfiles/${CHROM}"
    printf 'final_variants\t%s\n' "${observed}"
    printf 'final_samples\t%s\n' "${samples}"
} > "${OUTDIR}/${CHROM}.summary.tsv"

log "=== done: ${observed} variants, ${samples} samples ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
