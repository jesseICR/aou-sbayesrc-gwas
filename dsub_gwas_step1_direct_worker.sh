#!/bin/bash
# dsub_gwas_step1_direct_worker.sh - Build final REGENIE Step 1 direct bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[gwas-step1-direct $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/gwas_step1_direct
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/source.bed"
ln -sf "${BIM}" "${SCRATCH}/source.bim"
ln -sf "${FAM}" "${SCRATCH}/source.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "extract variants = $(wc -l < "${EXTRACT}")"
df -h /mnt/data | sed 's/^/  /'

"${PLINK2}" \
    --bfile "${SCRATCH}/source" \
    --extract "${EXTRACT}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${OUTDIR}/chr1_22_merged_gwas_step1"

observed=$(wc -l < "${OUTDIR}/chr1_22_merged_gwas_step1.bim")
samples=$(wc -l < "${OUTDIR}/chr1_22_merged_gwas_step1.fam")
if [[ "${observed}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: final Step 1 bfile has ${observed} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi

cp "${PARAMS}" "${OUTDIR}/gwas_step1_direct.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'final_bfile\t%s\n' "gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1"
    printf 'final_variants\t%s\n' "${observed}"
    printf 'final_samples\t%s\n' "${samples}"
} > "${OUTDIR}/gwas_step1_direct.summary.tsv"

log "=== done: ${observed} variants, ${samples} samples ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
