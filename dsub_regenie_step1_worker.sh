#!/bin/bash
# dsub_regenie_step1_worker.sh - Fit REGENIE Step 1 null model.

set -euo pipefail

: "${REGENIE_BUNDLE:?REGENIE_BUNDLE not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${PHEN:?PHEN not set}"
: "${COVAR:?COVAR not set}"
: "${KEEP:?KEEP not set}"
: "${PARAMS:?PARAMS not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${PHENO_COL:?PHENO_COL not set}"
: "${COVAR_COLS:?COVAR_COLS not set}"
: "${APPLY_RINT:?APPLY_RINT not set}"
: "${STEP1_BLOCK_SIZE:?STEP1_BLOCK_SIZE not set}"
: "${EXPECTED_KEEP_SAMPLES:?EXPECTED_KEEP_SAMPLES not set}"
: "${EXPECTED_BFILE_VARIANTS:?EXPECTED_BFILE_VARIANTS not set}"
: "${EXPECTED_BFILE_SAMPLES:?EXPECTED_BFILE_SAMPLES not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[regenie-step1 $(ts)] $*"; }

SCRATCH=/mnt/data/scratch/regenie_step1
mkdir -p "${SCRATCH}/regenie" "${OUTDIR}"
tar -xzf "${REGENIE_BUNDLE}" -C "${SCRATCH}/regenie"
export LD_LIBRARY_PATH="${SCRATCH}/regenie:${LD_LIBRARY_PATH:-}"
REGENIE_BIN="${SCRATCH}/regenie/regenie"
chmod +x "${REGENIE_BIN}"

ln -sf "${BED}" "${SCRATCH}/direct_hq.bed"
ln -sf "${BIM}" "${SCRATCH}/direct_hq.bim"
ln -sf "${FAM}" "${SCRATCH}/direct_hq.fam"

log "=== starting on $(hostname) ==="
log "regenie = $("${REGENIE_BIN}" --help 2>&1 | head -2 | tail -1 | sed 's/^ *//' || true)"
log "PHENO_COL=${PHENO_COL}"
log "COVAR_COLS=${COVAR_COLS}"
log "APPLY_RINT=${APPLY_RINT}"
log "STEP1_BLOCK_SIZE=${STEP1_BLOCK_SIZE}"
df -h /mnt/data | sed 's/^/  /'

keep_samples=$(wc -l < "${KEEP}")
bfile_variants=$(wc -l < "${SCRATCH}/direct_hq.bim")
bfile_samples=$(wc -l < "${SCRATCH}/direct_hq.fam")
if [[ "${keep_samples}" -ne "${EXPECTED_KEEP_SAMPLES}" ]]; then
    log "ERROR: keep file has ${keep_samples} samples, expected ${EXPECTED_KEEP_SAMPLES}"
    exit 1
fi
if [[ "${bfile_variants}" -ne "${EXPECTED_BFILE_VARIANTS}" || "${bfile_samples}" -ne "${EXPECTED_BFILE_SAMPLES}" ]]; then
    log "ERROR: bfile counts variants/samples=${bfile_variants}/${bfile_samples}, expected ${EXPECTED_BFILE_VARIANTS}/${EXPECTED_BFILE_SAMPLES}"
    exit 1
fi

cmd=(
    "${REGENIE_BIN}"
    --step 1
    --bed "${SCRATCH}/direct_hq"
    --phenoFile "${PHEN}"
    --phenoCol "${PHENO_COL}"
    --covarFile "${COVAR}"
    --covarColList "${COVAR_COLS}"
    --keep "${KEEP}"
    --qt
    --bsize "${STEP1_BLOCK_SIZE}"
    --threads "$(nproc)"
    --lowmem
    --lowmem-prefix "${SCRATCH}/height_step1_l0"
    --use-relative-path
    --out "${OUTDIR}/height_step1"
)
if [[ "${APPLY_RINT}" == "1" ]]; then
    cmd+=(--apply-rint)
fi

log "running REGENIE Step 1"
"${cmd[@]}"

pred_list="${OUTDIR}/height_step1_pred.list"
if [[ ! -s "${pred_list}" ]]; then
    log "ERROR: missing REGENIE prediction list ${pred_list}"
    ls -lh "${OUTDIR}" | sed 's/^/  /'
    exit 1
fi
pred_files=$(wc -l < "${pred_list}")
if [[ "${pred_files}" -le 0 ]]; then
    log "ERROR: prediction list is empty"
    exit 1
fi

cp "${PARAMS}" "${OUTDIR}/regenie_step1.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'step1_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'step1_bfile_variants\t%s\n' "${bfile_variants}"
    printf 'step1_bfile_samples\t%s\n' "${bfile_samples}"
    printf 'keep_samples\t%s\n' "${keep_samples}"
    printf 'pheno_col\t%s\n' "${PHENO_COL}"
    printf 'covar_cols\t%s\n' "${COVAR_COLS}"
    printf 'apply_rint\t%s\n' "${APPLY_RINT}"
    printf 'step1_block_size\t%s\n' "${STEP1_BLOCK_SIZE}"
    printf 'prediction_files\t%s\n' "${pred_files}"
} > "${OUTDIR}/regenie_step1.summary.tsv"

log "=== done Step 1: ${keep_samples} samples, ${bfile_variants} variants ==="
cat "${OUTDIR}/regenie_step1.summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
