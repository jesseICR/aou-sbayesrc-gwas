#!/bin/bash
# dsub_kinship_subset_worker.sh - Compute all-sample missingness for kinship SNPs.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${KINSHIP_MISSING_MAX:?KINSHIP_MISSING_MAX not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[kinship-subset $(ts)] $*"; }

param_value() {
    awk -F'\t' -v key="$1" '$1 == key {print $2; exit}' "${PARAMS}"
}

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/kinship_subset
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/hq.bed"
ln -sf "${BIM}" "${SCRATCH}/hq.bim"
ln -sf "${FAM}" "${SCRATCH}/hq.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "KINSHIP_MISSING_MAX=${KINSHIP_MISSING_MAX}"
log "extract SNPs before missingness = $(wc -l < "${EXTRACT}")"
df -h /mnt/data | sed 's/^/  /'

out_prefix="${OUTDIR}/kinship_snp_subset_all_sample_missingness"
log "computing all-sample variant missingness"
"${PLINK2}" \
    --bfile "${SCRATCH}/hq" \
    --extract "${EXTRACT}" \
    --missing variant-only \
    --threads "$(nproc)" \
    --out "${out_prefix}"

final_extract="${OUTDIR}/ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt"
threshold_counts="${OUTDIR}/kinship_snp_missingness_threshold_counts.tsv"
awk -v max_missing="${KINSHIP_MISSING_MAX}" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "ID") id_col = i
            if ($i == "F_MISS") fmiss_col = i
        }
        next
    }
    id_col > 0 && fmiss_col > 0 {
        total++
        fmiss = $fmiss_col + 0
        if (fmiss < max_missing) {
            print $id_col
            pass++
        }
    }
    END {
        if (id_col == 0 || fmiss_col == 0) {
            exit 2
        }
    }
' "${out_prefix}.vmiss" > "${final_extract}"

awk '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "F_MISS") fmiss_col = i
        }
        next
    }
    fmiss_col > 0 {
        total++
        fmiss = $fmiss_col + 0
        if (fmiss < 0.05) pass_005++
        if (fmiss < 0.04) pass_004++
        if (fmiss < 0.03) pass_003++
        if (fmiss < 0.02) pass_002++
        if (fmiss < 0.01) pass_001++
    }
    END {
        if (fmiss_col == 0) {
            exit 2
        }
        print "threshold\tpassing_variants\tfailing_variants"
        printf "missingness_lt_0.05\t%d\t%d\n", pass_005, total - pass_005
        printf "missingness_lt_0.04\t%d\t%d\n", pass_004, total - pass_004
        printf "missingness_lt_0.03\t%d\t%d\n", pass_003, total - pass_003
        printf "missingness_lt_0.02\t%d\t%d\n", pass_002, total - pass_002
        printf "missingness_lt_0.01\t%d\t%d\n", pass_001, total - pass_001
        printf "total_measured\t%d\t0\n", total
    }
' "${out_prefix}.vmiss" > "${threshold_counts}"

if [[ ! -s "${final_extract}" ]]; then
    log "ERROR: final kinship SNP extract is empty"
    exit 1
fi

n_final=$(wc -l < "${final_extract}")
n_vmiss=$(awk 'NR > 1 {n++} END {print n + 0}' "${out_prefix}.vmiss")
source_variants=$(param_value "source_bfile_variants")
source_samples=$(param_value "source_bfile_samples")
ukb_relatedness_snps=$(param_value "ukb_relatedness_snps")
n_intersection=$(param_value "n_intersection_hq_direct")

cp "${PARAMS}" "${OUTDIR}/kinship_snp_subset.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'source_bfile_variants\t%s\n' "${source_variants}"
    printf 'source_bfile_samples\t%s\n' "${source_samples}"
    printf 'ukb_relatedness_snps\t%s\n' "${ukb_relatedness_snps}"
    printf 'n_intersection_hq_direct\t%s\n' "${n_intersection}"
    printf 'kinship_missing_max_exclusive\t%s\n' "${KINSHIP_MISSING_MAX}"
    printf 'n_intersection_missingness_measured\t%s\n' "${n_vmiss}"
    printf 'n_intersection_and_missing_lt_%s\t%s\n' "${KINSHIP_MISSING_MAX}" "${n_final}"
    printf 'n_intersection_and_missing_lt_threshold\t%s\n' "${n_final}"
} > "${OUTDIR}/kinship_snp_subset_summary.tsv"

log "=== done: ${n_final}/${n_intersection} SNPs pass missingness < ${KINSHIP_MISSING_MAX} ==="
cat "${OUTDIR}/kinship_snp_subset_summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
