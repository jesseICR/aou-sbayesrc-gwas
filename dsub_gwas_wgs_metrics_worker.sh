#!/bin/bash
# dsub_gwas_wgs_metrics_worker.sh - QC-panel metrics for one GWAS Step 2 WGS pfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${PGEN:?PGEN not set}"
: "${PVAR:?PVAR not set}"
: "${PSAM:?PSAM not set}"
: "${FIT_PCA_KEEP:?FIT_PCA_KEEP not set}"
: "${EUR_KEEP:?EUR_KEEP not set}"
: "${CHROM:?CHROM not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[gwas-wgs-metrics ${CHROM} $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH="/mnt/data/scratch/gwas_wgs_metrics_${CHROM}"
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${PGEN}" "${SCRATCH}/${CHROM}.pgen"
ln -sf "${PVAR}" "${SCRATCH}/${CHROM}.pvar"
ln -sf "${PSAM}" "${SCRATCH}/${CHROM}.psam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "fit_pca samples = $(wc -l < "${FIT_PCA_KEEP}")"
log "classified EUR samples = $(wc -l < "${EUR_KEEP}")"
log "source pvar header = $(awk '/^#CHROM/ {print; exit}' "${SCRATCH}/${CHROM}.pvar")"
df -h /mnt/data | sed 's/^/  /'

log "computing fit_pca allele counts"
"${PLINK2}" \
    --pfile "${SCRATCH}/${CHROM}" \
    --keep "${FIT_PCA_KEEP}" \
    --freq counts \
    --threads "$(nproc)" \
    --out "${OUTDIR}/${CHROM}.fit_pca"

log "computing classified-EUR variant missingness"
"${PLINK2}" \
    --pfile "${SCRATCH}/${CHROM}" \
    --keep "${EUR_KEEP}" \
    --missing variant-only \
    --threads "$(nproc)" \
    --out "${OUTDIR}/${CHROM}.our_eur"

{
    printf 'chrom\tmetric\tvalue\n'
    printf '%s\tsource_variants\t%s\n' "${CHROM}" "$(grep -vc '^#' "${SCRATCH}/${CHROM}.pvar")"
    printf '%s\tsource_samples\t%s\n' "${CHROM}" "$(grep -vc '^#' "${SCRATCH}/${CHROM}.psam")"
    printf '%s\tfit_pca_keep_samples\t%s\n' "${CHROM}" "$(wc -l < "${FIT_PCA_KEEP}")"
    printf '%s\teur_keep_samples\t%s\n' "${CHROM}" "$(wc -l < "${EUR_KEEP}")"
    printf '%s\tfit_pca_acount_lines\t%s\n' "${CHROM}" "$(wc -l < "${OUTDIR}/${CHROM}.fit_pca.acount")"
    printf '%s\tour_eur_vmiss_lines\t%s\n' "${CHROM}" "$(wc -l < "${OUTDIR}/${CHROM}.our_eur.vmiss")"
} > "${OUTDIR}/${CHROM}.metrics_summary.tsv"

log "=== done ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
