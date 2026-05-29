#!/bin/bash
# dsub_fit_project_pca_worker.sh - Fit PCA and project allele weights to all samples.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${PCA_BED:?PCA_BED not set}"
: "${PCA_BIM:?PCA_BIM not set}"
: "${PCA_FAM:?PCA_FAM not set}"
: "${DIRECT_BED:?DIRECT_BED not set}"
: "${DIRECT_BIM:?DIRECT_BIM not set}"
: "${DIRECT_FAM:?DIRECT_FAM not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${PCA_NPCS:?PCA_NPCS not set}"
: "${PCA_SEED:?PCA_SEED not set}"
: "${EXPECTED_PCA_VARIANTS:?EXPECTED_PCA_VARIANTS not set}"
: "${EXPECTED_PCA_SAMPLES:?EXPECTED_PCA_SAMPLES not set}"
: "${EXPECTED_PROJECTED_SAMPLES:?EXPECTED_PROJECTED_SAMPLES not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[fit-project-pca $(ts)] $*"; }

count_rows_after_header() {
    local file="$1"
    local lines
    lines=$(wc -l < "${file}")
    echo $((lines - 1))
}

count_unique_score_ids() {
    awk 'NR > 1 {ids[$2] = 1} END {print length(ids)}' "$1"
}

chmod +x "${PLINK2}"

SCRATCH=/mnt/data/scratch/fit_project_pca
mkdir -p "${SCRATCH}" "${OUTDIR}" "${OUTDIR}/scrap"
ln -sf "${PCA_BED}" "${SCRATCH}/pca_ready.bed"
ln -sf "${PCA_BIM}" "${SCRATCH}/pca_ready.bim"
ln -sf "${PCA_FAM}" "${SCRATCH}/pca_ready.fam"
ln -sf "${DIRECT_BED}" "${SCRATCH}/direct_hq.bed"
ln -sf "${DIRECT_BIM}" "${SCRATCH}/direct_hq.bim"
ln -sf "${DIRECT_FAM}" "${SCRATCH}/direct_hq.fam"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "PCA_NPCS=${PCA_NPCS}"
log "PCA_SEED=${PCA_SEED}"
df -h /mnt/data | sed 's/^/  /'

pca_variants=$(wc -l < "${SCRATCH}/pca_ready.bim")
pca_samples=$(wc -l < "${SCRATCH}/pca_ready.fam")
direct_samples=$(wc -l < "${SCRATCH}/direct_hq.fam")
if [[ "${pca_variants}" -ne "${EXPECTED_PCA_VARIANTS}" ]]; then
    log "ERROR: pca_ready has ${pca_variants} variants, expected ${EXPECTED_PCA_VARIANTS}"
    exit 1
fi
if [[ "${pca_samples}" -ne "${EXPECTED_PCA_SAMPLES}" ]]; then
    log "ERROR: pca_ready has ${pca_samples} samples, expected ${EXPECTED_PCA_SAMPLES}"
    exit 1
fi
if [[ "${direct_samples}" -ne "${EXPECTED_PROJECTED_SAMPLES}" ]]; then
    log "ERROR: projection bfile has ${direct_samples} samples, expected ${EXPECTED_PROJECTED_SAMPLES}"
    exit 1
fi

log "Step 1: fit PCA on unrelated Europeans"
"${PLINK2}" \
    --bfile "${SCRATCH}/pca_ready" \
    --pca allele-wts "${PCA_NPCS}" approx \
    --seed "${PCA_SEED}" \
    --threads "$(nproc)" \
    --out "${OUTDIR}/aou_pcs"

if [[ ! -s "${OUTDIR}/aou_pcs.eigenvec.allele" || ! -s "${OUTDIR}/aou_pcs.eigenvec" || ! -s "${OUTDIR}/aou_pcs.eigenval" ]]; then
    log "ERROR: missing PCA output files"
    exit 1
fi

allele_weight_rows=$(count_rows_after_header "${OUTDIR}/aou_pcs.eigenvec.allele")
fit_variant_weights=$(count_unique_score_ids "${OUTDIR}/aou_pcs.eigenvec.allele")
fit_samples=$(count_rows_after_header "${OUTDIR}/aou_pcs.eigenvec")
eigenvalues=$(wc -l < "${OUTDIR}/aou_pcs.eigenval")
if [[ "${fit_variant_weights}" -ne "${pca_variants}" ]]; then
    log "ERROR: PCA allele weights have ${fit_variant_weights} unique variant IDs, expected ${pca_variants}"
    exit 1
fi
if [[ "${fit_samples}" -ne "${pca_samples}" ]]; then
    log "ERROR: PCA eigenvectors have ${fit_samples} sample rows, expected ${pca_samples}"
    exit 1
fi
if [[ "${eigenvalues}" -ne "${PCA_NPCS}" ]]; then
    log "ERROR: PCA eigenvalue count is ${eigenvalues}, expected ${PCA_NPCS}"
    exit 1
fi

log "Eigenvalues:"
sed 's/^/  /' "${OUTDIR}/aou_pcs.eigenval"

log "Step 2: compute allele frequency counts in the PCA fit set"
"${PLINK2}" \
    --bfile "${SCRATCH}/pca_ready" \
    --freq counts \
    --threads "$(nproc)" \
    --out "${OUTDIR}/pca_eur_counts"
if [[ ! -s "${OUTDIR}/pca_eur_counts.acount" ]]; then
    log "ERROR: missing allele count output"
    exit 1
fi

log "Step 3: project PCs onto all samples in direct_bfile_hq"
awk 'NR > 1 {print $2}' "${OUTDIR}/aou_pcs.eigenvec.allele" | sort -u > "${SCRATCH}/pca_snps.txt"
pca_score_snps=$(wc -l < "${SCRATCH}/pca_snps.txt")

LC_ALL=C comm -23 \
    <(LC_ALL=C sort "${SCRATCH}/pca_snps.txt") \
    <(awk '{print $2}' "${SCRATCH}/direct_hq.bim" | LC_ALL=C sort) \
    > "${OUTDIR}/scrap/pca_projection_missing_snps.txt"
missing_projection_snps=$(wc -l < "${OUTDIR}/scrap/pca_projection_missing_snps.txt")
if [[ "${missing_projection_snps}" -ne 0 ]]; then
    log "ERROR: ${missing_projection_snps} PCA SNPs are absent from the projection bfile"
    head "${OUTDIR}/scrap/pca_projection_missing_snps.txt" | sed 's/^/  /'
    exit 1
fi

read -r a1_col first_pc_col last_pc_col < <(
    awk -v last_pc="PC${PCA_NPCS}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "A1") a1 = i
                if ($i == "PC1") pc1 = i
                if ($i == last_pc) pclast = i
            }
            print a1, pc1, pclast
            exit
        }
    ' "${OUTDIR}/aou_pcs.eigenvec.allele"
)
if [[ -z "${a1_col}" || -z "${first_pc_col}" || -z "${last_pc_col}" ||
      "${a1_col}" -le 0 || "${first_pc_col}" -le 0 || "${last_pc_col}" -le 0 ]]; then
    log "ERROR: could not parse A1/PC columns from aou_pcs.eigenvec.allele"
    head -1 "${OUTDIR}/aou_pcs.eigenvec.allele" | sed 's/^/  /'
    exit 1
fi
log "A1 column=${a1_col}; score columns=${first_pc_col}-${last_pc_col}; score SNPs=${pca_score_snps}"

"${PLINK2}" \
    --bfile "${SCRATCH}/direct_hq" \
    --extract "${SCRATCH}/pca_snps.txt" \
    --read-freq "${OUTDIR}/pca_eur_counts.acount" \
    --score "${OUTDIR}/aou_pcs.eigenvec.allele" 2 "${a1_col}" header-read no-mean-imputation variance-standardize \
    --score-col-nums "${first_pc_col}-${last_pc_col}" \
    --threads "$(nproc)" \
    --out "${OUTDIR}/aou_projected"

if [[ ! -s "${OUTDIR}/aou_projected.sscore" ]]; then
    log "ERROR: missing projected score output"
    exit 1
fi
projected_samples=$(count_rows_after_header "${OUTDIR}/aou_projected.sscore")
if [[ "${projected_samples}" -ne "${EXPECTED_PROJECTED_SAMPLES}" ]]; then
    log "ERROR: projected scores have ${projected_samples} rows, expected ${EXPECTED_PROJECTED_SAMPLES}"
    exit 1
fi

cp "${PARAMS}" "${OUTDIR}/fit_project_pca.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'pca_bfile\t%s\n' "pca_eur/pca_ready"
    printf 'projection_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'pca_npcs\t%s\n' "${PCA_NPCS}"
    printf 'pca_seed\t%s\n' "${PCA_SEED}"
    printf 'fit_variants\t%s\n' "${pca_variants}"
    printf 'fit_samples\t%s\n' "${pca_samples}"
    printf 'allele_weight_rows\t%s\n' "${allele_weight_rows}"
    printf 'allele_weight_unique_variant_ids\t%s\n' "${fit_variant_weights}"
    printf 'eigenvalues\t%s\n' "${eigenvalues}"
    printf 'pca_score_snps\t%s\n' "${pca_score_snps}"
    printf 'missing_projection_snps\t%s\n' "${missing_projection_snps}"
    printf 'projected_samples\t%s\n' "${projected_samples}"
    printf 'a1_col\t%s\n' "${a1_col}"
    printf 'first_pc_col\t%s\n' "${first_pc_col}"
    printf 'last_pc_col\t%s\n' "${last_pc_col}"
} > "${OUTDIR}/fit_project_pca.summary.tsv"

log "=== done: ${pca_variants} variants, ${pca_samples} fit samples, ${projected_samples} projected samples ==="
cat "${OUTDIR}/fit_project_pca.summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
