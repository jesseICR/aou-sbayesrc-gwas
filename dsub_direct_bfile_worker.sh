#!/bin/bash
# dsub_direct_bfile_worker.sh — Merge direct-SNP pfiles into one REGENIE step-1 bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${SUMMARY:?SUMMARY not set}"
: "${DIRECT_PFILES:?DIRECT_PFILES not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[direct-bfile $(ts)] $*"; }

chmod +x "${PLINK2}"

DIRECT_PFILES="${DIRECT_PFILES%/}"
SCRATCH=/mnt/data/scratch/direct_bfile
mkdir -p "${SCRATCH}" "${OUTDIR}"

requested_direct=$(awk -F'\t' 'NR > 1 { s += $2 } END { print s + 0 }' "${SUMMARY}")
expected_present=$(awk -F'\t' 'NR > 1 { s += $3 } END { print s + 0 }' "${SUMMARY}")
expected_missing=$(awk -F'\t' 'NR > 1 { s += $4 } END { print s + 0 }' "${SUMMARY}")
if [[ "${requested_direct}" -le 0 || "${expected_present}" -le 0 ]]; then
    log "ERROR: could not determine expected direct SNP count from ${SUMMARY}"
    exit 1
fi

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "requested_direct=${requested_direct} present=${expected_present} missing=${expected_missing}"
df -h /mnt/data | sed 's/^/  /'

for chrom in $(seq 1 22); do
    prefix="${DIRECT_PFILES}/chr${chrom}"
    for ext in pgen pvar psam; do
        if [[ ! -f "${prefix}.${ext}" ]]; then
            log "ERROR: missing ${prefix}.${ext}"
            exit 1
        fi
    done
done

merge_list="${SCRATCH}/direct_pmerge_list.txt"
for chrom in $(seq 2 22); do
    printf '%s/chr%s\n' "${DIRECT_PFILES}" "${chrom}"
done > "${merge_list}"

merged_pgen="${SCRATCH}/chr1_22_merged_pgen"
merged_bed="${SCRATCH}/chr1_22_merged"

log "merging chr1..22 direct pfiles into sorted intermediate pgen"
"${PLINK2}" \
    --pfile "${DIRECT_PFILES}/chr1" \
    --pmerge-list "${merge_list}" pfile \
    --make-pgen \
    --sort-vars \
    --indiv-sort none \
    --threads "$(nproc)" \
    --out "${merged_pgen}"

observed_pgen=$(grep -vc '^#' "${merged_pgen}.pvar")
if [[ "${observed_pgen}" -ne "${expected_present}" ]]; then
    log "ERROR: merged pgen has ${observed_pgen} variants, expected ${expected_present}"
    exit 1
fi

log "converting merged pgen to bed/bim/fam"
"${PLINK2}" \
    --pfile "${merged_pgen}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${merged_bed}"

observed_bed=$(wc -l < "${merged_bed}.bim")
if [[ "${observed_bed}" -ne "${expected_present}" ]]; then
    log "ERROR: merged bfile has ${observed_bed} variants, expected ${expected_present}"
    exit 1
fi
samples=$(wc -l < "${merged_bed}.fam")

{
    printf 'prefix\trequested_direct\tpresent_from_wgs_pfiles\tmissing_from_wgs_pfiles\tvariants\tsamples\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${OUTPUT_PREFIX:-${OUTDIR}/chr1_22_merged}" "${requested_direct}" "${expected_present}" \
        "${expected_missing}" "${observed_bed}" "${samples}"
} > "${SCRATCH}/chr1_22_merged.summary.tsv"

cp "${merged_bed}.bed" "${OUTDIR}/chr1_22_merged.bed"
cp "${merged_bed}.bim" "${OUTDIR}/chr1_22_merged.bim"
cp "${merged_bed}.fam" "${OUTDIR}/chr1_22_merged.fam"
cp "${SCRATCH}/chr1_22_merged.summary.tsv" "${OUTDIR}/chr1_22_merged.summary.tsv"
[[ -f "${merged_bed}.log" ]] && cp "${merged_bed}.log" "${OUTDIR}/chr1_22_merged.log"

log "=== done: ${observed_bed} variants, ${samples} samples ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
