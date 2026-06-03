#!/bin/bash
# dsub_gwas_direct_metrics_worker.sh - QC-panel metrics for GWAS Step 1 direct bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${FIT_PCA_KEEP:?FIT_PCA_KEEP not set}"
: "${EUR_KEEP:?EUR_KEEP not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[gwas-direct-metrics $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/gwas_direct_metrics
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/direct_hq.bed"
ln -sf "${BIM}" "${SCRATCH}/direct_hq.bim"
ln -sf "${FAM}" "${SCRATCH}/direct_hq.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "fit_pca samples = $(wc -l < "${FIT_PCA_KEEP}")"
log "classified EUR samples = $(wc -l < "${EUR_KEEP}")"
df -h /mnt/data | sed 's/^/  /'

log "computing fit_pca allele counts"
"${PLINK2}" \
    --bfile "${SCRATCH}/direct_hq" \
    --keep "${FIT_PCA_KEEP}" \
    --freq counts \
    --threads "$(nproc)" \
    --out "${OUTDIR}/direct_hq.fit_pca"

log "computing classified-EUR variant missingness"
"${PLINK2}" \
    --bfile "${SCRATCH}/direct_hq" \
    --keep "${EUR_KEEP}" \
    --missing variant-only \
    --threads "$(nproc)" \
    --out "${OUTDIR}/direct_hq.our_eur"

{
    printf 'metric\tvalue\n'
    printf 'source_variants\t%s\n' "$(wc -l < "${SCRATCH}/direct_hq.bim")"
    printf 'source_samples\t%s\n' "$(wc -l < "${SCRATCH}/direct_hq.fam")"
    printf 'fit_pca_keep_samples\t%s\n' "$(wc -l < "${FIT_PCA_KEEP}")"
    printf 'eur_keep_samples\t%s\n' "$(wc -l < "${EUR_KEEP}")"
    printf 'fit_pca_acount_lines\t%s\n' "$(wc -l < "${OUTDIR}/direct_hq.fit_pca.acount")"
    printf 'our_eur_vmiss_lines\t%s\n' "$(wc -l < "${OUTDIR}/direct_hq.our_eur.vmiss")"
} > "${OUTDIR}/direct_hq.metrics_summary.tsv"

log "=== done ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
