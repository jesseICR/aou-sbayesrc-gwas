#!/bin/bash
# dsub_hq_direct_metrics_worker.sh - EUR frequency/missingness metrics for the direct bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${KEEP:?KEEP not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[hq-direct-metrics $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/hq_direct_metrics
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/direct.bed"
ln -sf "${BIM}" "${SCRATCH}/direct.bim"
ln -sf "${FAM}" "${SCRATCH}/direct.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "EUR keep samples = $(wc -l < "${KEEP}")"
df -h /mnt/data | sed 's/^/  /'

out_prefix="${OUTDIR}/eur_direct_qc"
log "computing EUR allele frequencies and variant missingness"
"${PLINK2}" \
    --bfile "${SCRATCH}/direct" \
    --keep "${KEEP}" \
    --freq \
    --missing variant-only \
    --threads "$(nproc)" \
    --out "${out_prefix}"

afreq_lines=$(wc -l < "${out_prefix}.afreq")
vmiss_lines=$(wc -l < "${out_prefix}.vmiss")
{
    printf 'file\tlines\n'
    printf 'eur_direct_qc.afreq\t%s\n' "${afreq_lines}"
    printf 'eur_direct_qc.vmiss\t%s\n' "${vmiss_lines}"
} > "${OUTDIR}/eur_direct_qc.line_counts.tsv"

log "=== done: afreq_lines=${afreq_lines}, vmiss_lines=${vmiss_lines} ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
