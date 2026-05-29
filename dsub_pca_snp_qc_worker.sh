#!/bin/bash
# dsub_pca_snp_qc_worker.sh - Build PCA-ready SNP bfile.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${KEEP:?KEEP not set}"
: "${VARIANT_QC:?VARIANT_QC not set}"
: "${HIGH_LD_BED:?HIGH_LD_BED not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${PCA_AF_DIFF_MAX:?PCA_AF_DIFF_MAX not set}"
: "${PCA_MAF_MIN:?PCA_MAF_MIN not set}"
: "${PCA_GENO_MAX:?PCA_GENO_MAX not set}"
: "${PCA_MIND_MAX:?PCA_MIND_MAX not set}"
: "${PCA_LD_WINDOW:?PCA_LD_WINDOW not set}"
: "${PCA_LD_STEP:?PCA_LD_STEP not set}"
: "${PCA_LD_R2:?PCA_LD_R2 not set}"
: "${EXPECTED_SOURCE_VARIANTS:?EXPECTED_SOURCE_VARIANTS not set}"
: "${EXPECTED_SOURCE_SAMPLES:?EXPECTED_SOURCE_SAMPLES not set}"
: "${EXPECTED_FIT_SAMPLES:?EXPECTED_FIT_SAMPLES not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[pca-snp-qc $(ts)] $*"; }

count_variants() { wc -l < "${1}.bim"; }
count_samples() { wc -l < "${1}.fam"; }

append_step() {
    local step="$1"
    local filter="$2"
    local in_variants="$3"
    local in_samples="$4"
    local out_variants="$5"
    local out_samples="$6"
    local desc="$7"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${step}" "${filter}" "${in_variants}" "${in_samples}" \
        "$((in_variants - out_variants))" "$((in_samples - out_samples))" \
        "${out_variants}" "${out_samples}" "${desc}" >> "${FILTER_STEPS}"
}

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/pca_snp_qc
mkdir -p "${SCRATCH}" "${OUTDIR}" "${OUTDIR}/scrap"
ln -sf "${BED}" "${SCRATCH}/source.bed"
ln -sf "${BIM}" "${SCRATCH}/source.bim"
ln -sf "${FAM}" "${SCRATCH}/source.fam"

FILTER_STEPS="${OUTDIR}/pca_snp_qc.filter_steps.tsv"
printf 'step\tfilter\tinput_variants\tinput_samples\tdropped_variants\tdropped_samples\tremaining_variants\tremaining_samples\tdescription\n' > "${FILTER_STEPS}"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "PCA_AF_DIFF_MAX=${PCA_AF_DIFF_MAX}"
log "PCA_MAF_MIN=${PCA_MAF_MIN}"
log "PCA_GENO_MAX=${PCA_GENO_MAX}"
log "PCA_MIND_MAX=${PCA_MIND_MAX}"
log "PCA_LD_WINDOW=${PCA_LD_WINDOW}"
log "PCA_LD_STEP=${PCA_LD_STEP}"
log "PCA_LD_R2=${PCA_LD_R2}"
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
keep_samples=$(wc -l < "${KEEP}")
if [[ "${keep_samples}" -ne "${EXPECTED_FIT_SAMPLES}" ]]; then
    log "ERROR: keep file has ${keep_samples} samples, expected ${EXPECTED_FIT_SAMPLES}"
    exit 1
fi

log "building rsid extract for abs ALT-frequency difference <= ${PCA_AF_DIFF_MAX}"
awk -F'\t' -v max="${PCA_AF_DIFF_MAX}" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "rsid") rsid_col = i
            if ($i == "abs_alt_freq_diff") diff_col = i
            if ($i == "pass_hq_direct") pass_col = i
        }
        next
    }
    rsid_col > 0 && diff_col > 0 && pass_col > 0 {
        pass_value = tolower($pass_col)
        if ((pass_value == "true" || pass_value == "1") && ($diff_col + 0) <= max) {
            print $rsid_col
        }
    }
    END {
        if (rsid_col == 0 || diff_col == 0 || pass_col == 0) {
            exit 2
        }
    }
' "${VARIANT_QC}" > "${SCRATCH}/pca_afdiff_pass.extract.txt"
afdiff_pass=$(wc -l < "${SCRATCH}/pca_afdiff_pass.extract.txt")
if [[ "${afdiff_pass}" -le 0 ]]; then
    log "ERROR: AF-difference extract list is empty"
    exit 1
fi

log "converting high-LD BED to PLINK range format"
awk 'BEGIN {OFS="\t"} { chr=$1; sub(/^chr/, "", chr); print chr, $2 + 1, $3, $4 }' \
    "${HIGH_LD_BED}" > "${SCRATCH}/high_ld_ranges.txt"

printf '0\tsource_hq_direct\t%s\t%s\t0\t0\t%s\t%s\tsource direct_bfile_hq before PCA filtering\n' \
    "${source_variants}" "${source_samples}" "${source_variants}" "${source_samples}" >> "${FILTER_STEPS}"

log "Step 1: subset HQ direct bfile to PCA-fitting IIDs"
"${PLINK2}" \
    --bfile "${SCRATCH}/source" \
    --keep "${KEEP}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/01_fit_pca_subset"
v0="${source_variants}"
s0="${source_samples}"
v1=$(count_variants "${SCRATCH}/01_fit_pca_subset")
s1=$(count_samples "${SCRATCH}/01_fit_pca_subset")
append_step 1 "keep_fit_pca_iids" "${v0}" "${s0}" "${v1}" "${s1}" "keep final PCA-fitting IIDs"

log "Step 2: apply abs ALT-frequency difference <= ${PCA_AF_DIFF_MAX}"
"${PLINK2}" \
    --bfile "${SCRATCH}/01_fit_pca_subset" \
    --extract "${SCRATCH}/pca_afdiff_pass.extract.txt" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/02_afdiff_filtered"
v2=$(count_variants "${SCRATCH}/02_afdiff_filtered")
s2=$(count_samples "${SCRATCH}/02_afdiff_filtered")
append_step 2 "abs_alt_freq_diff_le_${PCA_AF_DIFF_MAX}" "${v1}" "${s1}" "${v2}" "${s2}" "tighten AoU EUR vs SBayesRC ALT-frequency agreement"
rm -f "${SCRATCH}/01_fit_pca_subset."{bed,bim,fam,log,nosex}

log "Step 3: apply MAF >= ${PCA_MAF_MIN}"
"${PLINK2}" \
    --bfile "${SCRATCH}/02_afdiff_filtered" \
    --maf "${PCA_MAF_MIN}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/03_maf_filtered"
v3=$(count_variants "${SCRATCH}/03_maf_filtered")
s3=$(count_samples "${SCRATCH}/03_maf_filtered")
append_step 3 "maf_ge_${PCA_MAF_MIN}" "${v2}" "${s2}" "${v3}" "${s3}" "drop variants with MAF below threshold in PCA-fitting samples"
rm -f "${SCRATCH}/02_afdiff_filtered."{bed,bim,fam,log,nosex}

log "Step 4: apply variant missingness <= ${PCA_GENO_MAX}"
"${PLINK2}" \
    --bfile "${SCRATCH}/03_maf_filtered" \
    --geno "${PCA_GENO_MAX}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/04_geno_filtered"
v4=$(count_variants "${SCRATCH}/04_geno_filtered")
s4=$(count_samples "${SCRATCH}/04_geno_filtered")
append_step 4 "geno_le_${PCA_GENO_MAX}" "${v3}" "${s3}" "${v4}" "${s4}" "drop variants with missingness above threshold"
rm -f "${SCRATCH}/03_maf_filtered."{bed,bim,fam,log,nosex}

log "Step 5: apply sample missingness <= ${PCA_MIND_MAX}"
"${PLINK2}" \
    --bfile "${SCRATCH}/04_geno_filtered" \
    --mind "${PCA_MIND_MAX}" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/05_mind_filtered"
v5=$(count_variants "${SCRATCH}/05_mind_filtered")
s5=$(count_samples "${SCRATCH}/05_mind_filtered")
append_step 5 "mind_le_${PCA_MIND_MAX}" "${v4}" "${s4}" "${v5}" "${s5}" "drop samples with missingness above threshold"
rm -f "${SCRATCH}/04_geno_filtered."{bed,bim,fam,log,nosex}

log "Step 6: exclude long-range LD regions"
"${PLINK2}" \
    --bfile "${SCRATCH}/05_mind_filtered" \
    --exclude range "${SCRATCH}/high_ld_ranges.txt" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${SCRATCH}/06_ldregion_excluded"
v6=$(count_variants "${SCRATCH}/06_ldregion_excluded")
s6=$(count_samples "${SCRATCH}/06_ldregion_excluded")
append_step 6 "exclude_long_range_ld_regions_hg38" "${v5}" "${s5}" "${v6}" "${s6}" "exclude Price et al. long-range LD regions, hg38"
rm -f "${SCRATCH}/05_mind_filtered."{bed,bim,fam,log,nosex}

log "Step 7: LD prune window=${PCA_LD_WINDOW}, step=${PCA_LD_STEP}, r2=${PCA_LD_R2}"
"${PLINK2}" \
    --bfile "${SCRATCH}/06_ldregion_excluded" \
    --indep-pairwise "${PCA_LD_WINDOW}" "${PCA_LD_STEP}" "${PCA_LD_R2}" \
    --threads "$(nproc)" \
    --out "${SCRATCH}/pca_ld_prune"
"${PLINK2}" \
    --bfile "${SCRATCH}/06_ldregion_excluded" \
    --extract "${SCRATCH}/pca_ld_prune.prune.in" \
    --make-bed \
    --threads "$(nproc)" \
    --out "${OUTDIR}/pca_ready"
v7=$(count_variants "${OUTDIR}/pca_ready")
s7=$(count_samples "${OUTDIR}/pca_ready")
append_step 7 "ld_prune_${PCA_LD_WINDOW}_${PCA_LD_STEP}_${PCA_LD_R2}" "${v6}" "${s6}" "${v7}" "${s7}" "LD prune retained variants"
rm -f "${SCRATCH}/06_ldregion_excluded."{bed,bim,fam,log,nosex}

cp "${SCRATCH}/pca_afdiff_pass.extract.txt" "${OUTDIR}/pca_afdiff_${PCA_AF_DIFF_MAX}.extract.txt"
cp "${SCRATCH}/pca_ld_prune.prune.in" "${OUTDIR}/pca_ld_prune.prune.in"
cp "${SCRATCH}/pca_ld_prune.prune.out" "${OUTDIR}/pca_ld_prune.prune.out"
cp "${PARAMS}" "${OUTDIR}/pca_snp_qc.params.tsv"

{
    printf 'metric\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'source_variants\t%s\n' "${source_variants}"
    printf 'source_samples\t%s\n' "${source_samples}"
    printf 'fit_pca_iids\t%s\n' "${keep_samples}"
    printf 'pca_af_diff_max\t%s\n' "${PCA_AF_DIFF_MAX}"
    printf 'pca_maf_min\t%s\n' "${PCA_MAF_MIN}"
    printf 'pca_geno_max\t%s\n' "${PCA_GENO_MAX}"
    printf 'pca_mind_max\t%s\n' "${PCA_MIND_MAX}"
    printf 'pca_ld_window\t%s\n' "${PCA_LD_WINDOW}"
    printf 'pca_ld_step\t%s\n' "${PCA_LD_STEP}"
    printf 'pca_ld_r2\t%s\n' "${PCA_LD_R2}"
    printf 'afdiff_extract_snps\t%s\n' "${afdiff_pass}"
    printf 'final_variants\t%s\n' "${v7}"
    printf 'final_samples\t%s\n' "${s7}"
} > "${OUTDIR}/pca_snp_qc.summary.tsv"

log "=== done: ${v7} PCA-ready variants, ${s7} samples ==="
cat "${OUTDIR}/pca_snp_qc.filter_steps.tsv" | sed 's/^/  /'
cat "${OUTDIR}/pca_snp_qc.summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
