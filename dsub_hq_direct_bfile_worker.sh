#!/bin/bash
# dsub_hq_direct_bfile_worker.sh - Build the filtered high-quality direct bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${KEEP:?KEEP not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[hq-direct-bfile $(ts)] $*"; }

summarize_smiss() {
    local label="$1"
    local smiss="$2"
    local out="$3"
    awk -v label="${label}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "F_MISS") {
                    fmiss_col = i
                }
            }
            next
        }
        fmiss_col > 0 {
            value = $fmiss_col + 0
            n++
            sum += value
            if (n == 1 || value < min) min = value
            if (n == 1 || value > max) max = value
        }
        END {
            if (n == 0) {
                printf "%s\t0\tNA\tNA\tNA\n", label
            } else {
                printf "%s\t%d\t%.10g\t%.10g\t%.10g\n", label, n, min, sum / n, max
            }
        }
    ' "${smiss}" >> "${out}"
}

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/hq_direct_bfile
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/direct.bed"
ln -sf "${BIM}" "${SCRATCH}/direct.bim"
ln -sf "${FAM}" "${SCRATCH}/direct.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "expected_variants=${EXPECTED_VARIANTS}"
log "EUR keep samples = $(wc -l < "${KEEP}")"
df -h /mnt/data | sed 's/^/  /'

out_prefix="${OUTDIR}/chr1_22_merged_hq"

log "building filtered high-quality direct bfile"
"${PLINK2}" \
    --bfile "${SCRATCH}/direct" \
    --extract "${EXTRACT}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${out_prefix}"

observed_variants=$(wc -l < "${out_prefix}.bim")
if [[ "${observed_variants}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: filtered bfile has ${observed_variants} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi
observed_samples=$(wc -l < "${out_prefix}.fam")

log "computing sample missingness over final variant set for all samples"
"${PLINK2}" \
    --bfile "${out_prefix}" \
    --missing sample-only \
    --threads "$(nproc)" \
    --out "${OUTDIR}/chr1_22_merged_hq.sample_missingness_all"

log "computing sample missingness over final variant set for EUR samples"
"${PLINK2}" \
    --bfile "${out_prefix}" \
    --keep "${KEEP}" \
    --missing sample-only \
    --threads "$(nproc)" \
    --out "${OUTDIR}/chr1_22_merged_hq.sample_missingness_eur"

sample_summary="${OUTDIR}/chr1_22_merged_hq.sample_missingness_summary.tsv"
printf 'sample_set\tsamples\tmin_f_miss\tmean_f_miss\tmax_f_miss\n' > "${sample_summary}"
summarize_smiss "all_samples" "${OUTDIR}/chr1_22_merged_hq.sample_missingness_all.smiss" "${sample_summary}"
summarize_smiss "eur_samples" "${OUTDIR}/chr1_22_merged_hq.sample_missingness_eur.smiss" "${sample_summary}"

{
    printf 'prefix\tvariants\tsamples\teur_samples\n'
    printf '%s\t%s\t%s\t%s\n' \
        "${OUTDIR}/chr1_22_merged_hq" \
        "${observed_variants}" \
        "${observed_samples}" \
        "$(wc -l < "${KEEP}")"
} > "${OUTDIR}/chr1_22_merged_hq.summary.tsv"

log "=== done: ${observed_variants} variants, ${observed_samples} samples ==="
cat "${sample_summary}" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
