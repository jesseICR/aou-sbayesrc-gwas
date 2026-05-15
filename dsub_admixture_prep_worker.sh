#!/bin/bash
# dsub_admixture_prep_worker.sh - Build ADMIXTURE K=6 aligned bfile + .P.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${REF_TSV:?REF_TSV not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${ADMIXTURE_GENO_MAX:?ADMIXTURE_GENO_MAX not set}"
: "${EXPECTED_SOURCE_VARIANTS:?EXPECTED_SOURCE_VARIANTS not set}"
: "${EXPECTED_SOURCE_SAMPLES:?EXPECTED_SOURCE_SAMPLES not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[admixture-prep $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/admixture_prep
mkdir -p "${SCRATCH}" "${OUTDIR}"
ln -sf "${BED}" "${SCRATCH}/source.bed"
ln -sf "${BIM}" "${SCRATCH}/source.bim"
ln -sf "${FAM}" "${SCRATCH}/source.fam"
cp "${REF_TSV}" "${SCRATCH}/admixture_allele_freqs.tsv"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "ADMIXTURE_GENO_MAX=${ADMIXTURE_GENO_MAX}"
df -h /mnt/data | sed 's/^/  /'

source_variants=$(wc -l < "${SCRATCH}/source.bim")
source_samples=$(wc -l < "${SCRATCH}/source.fam")
if [[ "${source_variants}" -ne "${EXPECTED_SOURCE_VARIANTS}" ]]; then
    log "ERROR: source bfile has ${source_variants} variants, expected ${EXPECTED_SOURCE_VARIANTS}"
    exit 1
fi
if [[ "${source_samples}" -ne "${EXPECTED_SOURCE_SAMPLES}" ]]; then
    log "ERROR: source bfile has ${source_samples} samples, expected ${EXPECTED_SOURCE_SAMPLES}"
    exit 1
fi

ref_variants=$(awk 'NR > 1 {n++} END {print n + 0}' "${SCRATCH}/admixture_allele_freqs.tsv")
awk -F'\t' 'NR > 1 {print $2}' "${SCRATCH}/admixture_allele_freqs.tsv" > "${SCRATCH}/ref_snp_ids.txt"

log "computing all-sample variant missingness on source HQ direct bfile"
"${PLINK2}" \
    --bfile "${SCRATCH}/source" \
    --missing variant-only \
    --threads "$(nproc)" \
    --out "${SCRATCH}/source_all_samples"

awk -v max="${ADMIXTURE_GENO_MAX}" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "ID") id_col = i
            if ($i == "F_MISS") fmiss_col = i
        }
        next
    }
    fmiss_col > 0 && ($fmiss_col + 0) <= max { print $id_col }
' "${SCRATCH}/source_all_samples.vmiss" > "${SCRATCH}/geno_pass_snps.txt"
geno_pass_variants=$(wc -l < "${SCRATCH}/geno_pass_snps.txt")

awk 'NR == FNR {ref[$1] = 1; next} ($1 in ref) {print $1}' \
    "${SCRATCH}/ref_snp_ids.txt" \
    "${SCRATCH}/geno_pass_snps.txt" > "${SCRATCH}/geno_pass_ref_snps.txt"
geno_pass_ref_variants=$(wc -l < "${SCRATCH}/geno_pass_ref_snps.txt")

log "extracting reference-overlap SNPs after --geno ${ADMIXTURE_GENO_MAX}"
"${PLINK2}" \
    --bfile "${SCRATCH}/source" \
    --extract "${SCRATCH}/geno_pass_ref_snps.txt" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/aou_admixture_extracted"
extracted_variants=$(wc -l < "${SCRATCH}/aou_admixture_extracted.bim")

log "aligning reference frequencies to extracted BIM allele order"
rm -f "${SCRATCH}/snps_aligned.txt" "${SCRATCH}/ref_aligned.P" \
      "${SCRATCH}/admixture_align_log.txt" "${SCRATCH}/admixture_align_summary.tsv"
awk -v OFS="\t" '
    function comp(a) {
        if (a == "A") return "T"
        if (a == "T") return "A"
        if (a == "C") return "G"
        if (a == "G") return "C"
        return ""
    }
    function ambiguous(a, b) {
        return ((a == "A" && b == "T") || (a == "T" && b == "A") ||
                (a == "C" && b == "G") || (a == "G" && b == "C"))
    }
    function write_freqs(rsid, f1, f2, f3, f4, f5, f6) {
        print rsid >> snps_out
        printf "%.10g %.10g %.10g %.10g %.10g %.10g\n", f1, f2, f3, f4, f5, f6 >> p_out
    }
    BEGIN {
        snps_out = "'"${SCRATCH}"'/snps_aligned.txt"
        p_out = "'"${SCRATCH}"'/ref_aligned.P"
        log_out = "'"${SCRATCH}"'/admixture_align_log.txt"
        summary_out = "'"${SCRATCH}"'/admixture_align_summary.tsv"
    }
    NR == FNR && FNR > 1 {
        rsid = $2
        ref_a1[rsid] = $4
        ref_a2[rsid] = $5
        for (i = 1; i <= 6; i++) {
            ref_f[rsid, i] = $(5 + i) + 0
        }
        reference_variants++
        next
    }
    NR != FNR {
        bim_variants++
        rsid = $2
        bim_a1 = $5
        bim_a2 = $6
        if (!(rsid in ref_a1)) {
            missing_from_reference++
            next
        }
        if (ambiguous(bim_a1, bim_a2)) {
            strand_ambiguous++
            next
        }
        if (bim_a1 == ref_a1[rsid] && bim_a2 == ref_a2[rsid]) {
            same++
            retained++
            write_freqs(rsid, ref_f[rsid,1], ref_f[rsid,2], ref_f[rsid,3], ref_f[rsid,4], ref_f[rsid,5], ref_f[rsid,6])
            next
        }
        if (bim_a1 == ref_a2[rsid] && bim_a2 == ref_a1[rsid]) {
            swapped++
            retained++
            write_freqs(rsid, 1-ref_f[rsid,1], 1-ref_f[rsid,2], 1-ref_f[rsid,3], 1-ref_f[rsid,4], 1-ref_f[rsid,5], 1-ref_f[rsid,6])
            next
        }
        if (bim_a1 == comp(ref_a1[rsid]) && bim_a2 == comp(ref_a2[rsid])) {
            strand_flip++
            retained++
            write_freqs(rsid, ref_f[rsid,1], ref_f[rsid,2], ref_f[rsid,3], ref_f[rsid,4], ref_f[rsid,5], ref_f[rsid,6])
            next
        }
        if (bim_a1 == comp(ref_a2[rsid]) && bim_a2 == comp(ref_a1[rsid])) {
            swapped_strand_flip++
            retained++
            write_freqs(rsid, 1-ref_f[rsid,1], 1-ref_f[rsid,2], 1-ref_f[rsid,3], 1-ref_f[rsid,4], 1-ref_f[rsid,5], 1-ref_f[rsid,6])
            next
        }
        allele_mismatch++
    }
    END {
        excluded = bim_variants - retained
        printf "Allele alignment retained %d of %d extracted SNPs\n", retained, bim_variants > log_out
        printf "same=%d\n", same >> log_out
        printf "swapped=%d\n", swapped >> log_out
        printf "strand_flip=%d\n", strand_flip >> log_out
        printf "swapped_strand_flip=%d\n", swapped_strand_flip >> log_out
        printf "strand_ambiguous=%d\n", strand_ambiguous >> log_out
        printf "allele_mismatch=%d\n", allele_mismatch >> log_out
        printf "missing_from_reference=%d\n", missing_from_reference >> log_out
        printf "populations=European,East Asian,American,African,South Asian,Oceanian\n" >> log_out

        printf "metric\tvalue\n" > summary_out
        printf "bim_variants\t%d\n", bim_variants >> summary_out
        printf "reference_variants\t%d\n", reference_variants >> summary_out
        printf "same\t%d\n", same >> summary_out
        printf "swapped\t%d\n", swapped >> summary_out
        printf "strand_flip\t%d\n", strand_flip >> summary_out
        printf "swapped_strand_flip\t%d\n", swapped_strand_flip >> summary_out
        printf "strand_ambiguous\t%d\n", strand_ambiguous >> summary_out
        printf "allele_mismatch\t%d\n", allele_mismatch >> summary_out
        printf "missing_from_reference\t%d\n", missing_from_reference >> summary_out
        printf "retained\t%d\n", retained >> summary_out
        printf "excluded\t%d\n", excluded >> summary_out
        printf "populations\tEuropean,East Asian,American,African,South Asian,Oceanian\n" >> summary_out

        while ((getline line < log_out) > 0) print line
    }
' "${SCRATCH}/admixture_allele_freqs.tsv" "${SCRATCH}/aou_admixture_extracted.bim"
aligned_variants=$(wc -l < "${SCRATCH}/snps_aligned.txt")
p_rows=$(wc -l < "${SCRATCH}/ref_aligned.P")
if [[ "${aligned_variants}" -ne "${p_rows}" ]]; then
    log "ERROR: aligned SNP list has ${aligned_variants} rows but ref_aligned.P has ${p_rows}"
    exit 1
fi

log "building final aligned ADMIXTURE bfile"
"${PLINK2}" \
    --bfile "${SCRATCH}/aou_admixture_extracted" \
    --extract "${SCRATCH}/snps_aligned.txt" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${OUTDIR}/aou_admixture_aligned"

final_variants=$(wc -l < "${OUTDIR}/aou_admixture_aligned.bim")
final_samples=$(wc -l < "${OUTDIR}/aou_admixture_aligned.fam")
if [[ "${final_variants}" -ne "${aligned_variants}" ]]; then
    log "ERROR: final aligned bfile has ${final_variants} variants, expected ${aligned_variants}"
    exit 1
fi
if [[ "${final_samples}" -ne "${source_samples}" ]]; then
    log "ERROR: final aligned bfile has ${final_samples} samples, expected ${source_samples}"
    exit 1
fi

cp "${SCRATCH}/admixture_allele_freqs.tsv" "${OUTDIR}/admixture_allele_freqs.tsv"
cp "${SCRATCH}/ref_aligned.P" "${OUTDIR}/ref_aligned.P"
cp "${SCRATCH}/snps_aligned.txt" "${OUTDIR}/snps_aligned.txt"
cp "${SCRATCH}/admixture_align_log.txt" "${OUTDIR}/admixture_align_log.txt"
cp "${SCRATCH}/admixture_align_summary.tsv" "${OUTDIR}/admixture_align_summary.tsv"
cp "${SCRATCH}/source_all_samples.vmiss" "${OUTDIR}/aou_admixture_source_all_samples.vmiss"
cp "${SCRATCH}/source_all_samples.log" "${OUTDIR}/aou_admixture_source_all_samples.log"
cp "${PARAMS}" "${OUTDIR}/admixture_prep.params.tsv"

{
    printf 'metric\tvalue\n'
    printf 'source_hq_direct_variants\t%s\n' "${source_variants}"
    printf 'source_hq_direct_samples\t%s\n' "${source_samples}"
    printf 'reference_variants\t%s\n' "${ref_variants}"
    printf 'geno_missingness_max\t%s\n' "${ADMIXTURE_GENO_MAX}"
    printf 'variants_passing_all_sample_geno\t%s\n' "${geno_pass_variants}"
    printf 'variants_dropped_by_all_sample_geno\t%s\n' "$((source_variants - geno_pass_variants))"
    printf 'variants_after_reference_intersection\t%s\n' "${geno_pass_ref_variants}"
    printf 'variants_extracted_for_alignment\t%s\n' "${extracted_variants}"
    printf 'variants_dropped_by_allele_alignment\t%s\n' "$((extracted_variants - aligned_variants))"
    printf 'final_aligned_variants\t%s\n' "${final_variants}"
    printf 'final_samples\t%s\n' "${final_samples}"
    printf 'ref_aligned_p_rows\t%s\n' "${p_rows}"
} > "${OUTDIR}/admixture_prep_summary.tsv"

log "=== done: ${final_variants} aligned variants, ${final_samples} samples ==="
cat "${OUTDIR}/admixture_prep_summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
