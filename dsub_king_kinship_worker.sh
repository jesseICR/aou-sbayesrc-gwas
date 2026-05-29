#!/bin/bash
# dsub_king_kinship_worker.sh - Run plink2 KING kinship.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${KING_TABLE_FILTER:?KING_TABLE_FILTER not set}"
: "${EXPECTED_SNPS:?EXPECTED_SNPS not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[king-kinship $(ts)] $*"; }

param_value() {
    awk -F'\t' -v key="$1" '$1 == key {print $2; exit}' "${PARAMS}"
}

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/king_kinship
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/hq.bed"
ln -sf "${BIM}" "${SCRATCH}/hq.bim"
ln -sf "${FAM}" "${SCRATCH}/hq.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "KING_TABLE_FILTER=${KING_TABLE_FILTER}"
log "EXPECTED_SNPS=${EXPECTED_SNPS}"
df -h /mnt/data | sed 's/^/  /'

out_prefix="${OUTDIR}/aou_hq_direct_rel"
log "running plink2 --make-king-table"
"${PLINK2}" \
    --bfile "${SCRATCH}/hq" \
    --extract "${EXTRACT}" \
    --make-king-table \
    --king-table-filter "${KING_TABLE_FILTER}" \
    --threads "$(nproc)" \
    --out "${out_prefix}"

if [[ ! -s "${out_prefix}.kin0" ]]; then
    log "ERROR: missing KING output ${out_prefix}.kin0"
    exit 1
fi

king_pairs=$(( $(wc -l < "${out_prefix}.kin0") - 1 ))
source_variants=$(param_value "source_bfile_variants")
source_samples=$(param_value "source_bfile_samples")
kinship_snps=$(param_value "kinship_snps")
if [[ "${kinship_snps}" -ne "${EXPECTED_SNPS}" ]]; then
    log "ERROR: expected ${EXPECTED_SNPS} SNPs, params specify ${kinship_snps}"
    exit 1
fi

cp "${PARAMS}" "${OUTDIR}/aou_hq_direct_rel.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'source_bfile_variants\t%s\n' "${source_variants}"
    printf 'source_bfile_samples\t%s\n' "${source_samples}"
    printf 'kinship_snps\t%s\n' "${kinship_snps}"
    printf 'king_table_filter\t%s\n' "${KING_TABLE_FILTER}"
    printf 'king_pairs\t%s\n' "${king_pairs}"
} > "${OUTDIR}/aou_hq_direct_rel.summary.tsv"

log "=== done: ${king_pairs} pairs ==="
cat "${OUTDIR}/aou_hq_direct_rel.summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
