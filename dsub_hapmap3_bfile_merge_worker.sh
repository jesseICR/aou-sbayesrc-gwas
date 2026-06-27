#!/bin/bash
# dsub_hapmap3_bfile_merge_worker.sh - Merge HapMap3 HQ chromosome pfiles to one bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${PFILES:?PFILES not set}"
: "${PARAMS:?PARAMS not set}"
: "${FILTER_SUMMARY:?FILTER_SUMMARY not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"
: "${EXPECTED_SAMPLES:?EXPECTED_SAMPLES not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[hapmap3-merge $(ts)] $*"; }

chmod +x "${PLINK2}"

PFILES="${PFILES%/}"
SCRATCH="/mnt/data/scratch/hapmap3_merge"
mkdir -p "${SCRATCH}" "${OUTDIR}"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "expected_variants=${EXPECTED_VARIANTS} expected_samples=${EXPECTED_SAMPLES}"
df -h /mnt/data | sed 's/^/  /'

sum_variants=0
for chrom in $(seq 1 22); do
    prefix="${PFILES}/chr${chrom}"
    for ext in pgen pvar psam; do
        if [[ ! -s "${prefix}.${ext}" ]]; then
            log "ERROR: missing ${prefix}.${ext}"
            exit 1
        fi
    done
    chrom_variants=$(grep -vc '^#' "${prefix}.pvar")
    chrom_samples=$(grep -vc '^#' "${prefix}.psam")
    if [[ "${chrom_samples}" -ne "${EXPECTED_SAMPLES}" ]]; then
        log "ERROR: ${prefix}.psam has ${chrom_samples} samples, expected ${EXPECTED_SAMPLES}"
        exit 1
    fi
    sum_variants=$((sum_variants + chrom_variants))
done
if [[ "${sum_variants}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: chromosome pfiles sum to ${sum_variants} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi

merge_list="${SCRATCH}/hapmap3_pmerge_list.txt"
for chrom in $(seq 2 22); do
    printf '%s/chr%s\n' "${PFILES}" "${chrom}"
done > "${merge_list}"

merged_pgen="${SCRATCH}/hapmap3_bfile_hq_pgen"
merged_bed="${SCRATCH}/hapmap3_bfile_hq"

log "merging HapMap3 chromosome pfiles"
"${PLINK2}" \
    --pfile "${PFILES}/chr1" \
    --pmerge-list "${merge_list}" pfile \
    --make-pgen \
    --sort-vars \
    --indiv-sort none \
    --threads "$(nproc)" \
    --out "${merged_pgen}"

observed_pgen=$(grep -vc '^#' "${merged_pgen}.pvar")
if [[ "${observed_pgen}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: merged pgen has ${observed_pgen} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi

log "converting merged pgen to bed/bim/fam"
"${PLINK2}" \
    --pfile "${merged_pgen}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${merged_bed}"

observed_bim=$(wc -l < "${merged_bed}.bim")
observed_fam=$(wc -l < "${merged_bed}.fam")
if [[ "${observed_bim}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: final bfile has ${observed_bim} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi
if [[ "${observed_fam}" -ne "${EXPECTED_SAMPLES}" ]]; then
    log "ERROR: final bfile has ${observed_fam} samples, expected ${EXPECTED_SAMPLES}"
    exit 1
fi

cp "${merged_bed}.bed" "${OUTDIR}/hapmap3_bfile_hq.bed"
cp "${merged_bed}.bim" "${OUTDIR}/hapmap3_bfile_hq.bim"
cp "${merged_bed}.fam" "${OUTDIR}/hapmap3_bfile_hq.fam"
[[ -f "${merged_bed}.log" ]] && cp "${merged_bed}.log" "${OUTDIR}/hapmap3_bfile_hq.log"
cp "${PARAMS}" "${OUTDIR}/hapmap3_bfile_hq.params.tsv"
cp "${FILTER_SUMMARY}" "${OUTDIR}/hapmap3_bfile_hq.filter_summary.tsv"

log "computing EUR sample missingness on final HapMap3 HQ bfile"
"${PLINK2}" \
    --bfile "${OUTDIR}/hapmap3_bfile_hq" \
    --missing sample-only \
    --threads "$(nproc)" \
    --out "${OUTDIR}/hapmap3_bfile_hq.sample_missingness_eur"

{
    printf 'metric\tvalue\n'
    printf 'final_bfile\t%s\n' "hapmap3_bfile_hq/hapmap3_bfile_hq"
    printf 'final_variants\t%s\n' "${observed_bim}"
    printf 'final_samples\t%s\n' "${observed_fam}"
    printf 'source_pfiles\t%s\n' "wgs_pfiles/chr1-22"
} > "${OUTDIR}/hapmap3_bfile_hq.summary.tsv"

log "=== done: ${observed_bim} variants, ${observed_fam} samples ==="
ls -lh "${OUTDIR}" | sed 's/^/  /'
