#!/bin/bash
# dsub_regenie_step2_worker.sh - Run REGENIE Step 2 for one chromosome.

set -euo pipefail

: "${REGENIE_BUNDLE:?REGENIE_BUNDLE not set}"
: "${PGEN:?PGEN not set}"
: "${PVAR:?PVAR not set}"
: "${PSAM:?PSAM not set}"
: "${PHEN:?PHEN not set}"
: "${COVAR:?COVAR not set}"
: "${KEEP:?KEEP not set}"
: "${PARAMS:?PARAMS not set}"
: "${STEP1_DIR:?STEP1_DIR not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${CHROM:?CHROM not set}"
: "${PHENO_COL:?PHENO_COL not set}"
: "${COVAR_COLS:?COVAR_COLS not set}"
: "${APPLY_RINT:?APPLY_RINT not set}"
: "${STEP2_BLOCK_SIZE:?STEP2_BLOCK_SIZE not set}"
: "${EXPECTED_KEEP_SAMPLES:?EXPECTED_KEEP_SAMPLES not set}"
: "${EXPECTED_VARIANTS:?EXPECTED_VARIANTS not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[regenie-step2 ${CHROM} $(ts)] $*"; }

SCRATCH=/mnt/data/scratch/regenie_step2_${CHROM}
mkdir -p "${SCRATCH}/regenie" "${OUTDIR}"
tar -xzf "${REGENIE_BUNDLE}" -C "${SCRATCH}/regenie"
export LD_LIBRARY_PATH="${SCRATCH}/regenie:${LD_LIBRARY_PATH:-}"
REGENIE_BIN="${SCRATCH}/regenie/regenie"
chmod +x "${REGENIE_BIN}"

ln -sf "${PGEN}" "${SCRATCH}/${CHROM}.pgen"
ln -sf "${PVAR}" "${SCRATCH}/${CHROM}.pvar"

# AoU pfiles ship with #IID-only psam headers. REGENIE requires #FID IID...
# so create a local psam with FID=0, matching the fam/phen/covar files.
awk '
BEGIN { OFS = "\t" }
NR == 1 {
    if ($1 == "#FID") {
        print
        mode = "copy"
        next
    }
    if ($1 == "#IID") {
        $1 = "IID"
        print "#FID", $0
        mode = "add_fid"
        next
    }
    print "ERROR: unsupported psam header: " $0 > "/dev/stderr"
    exit 1
}
mode == "add_fid" {
    print "0", $0
    next
}
{
    print
}
' "${PSAM}" > "${SCRATCH}/${CHROM}.psam"

log "=== starting on $(hostname) ==="
log "regenie = $("${REGENIE_BIN}" --help 2>&1 | head -2 | tail -1 | sed 's/^ *//' || true)"
log "PHENO_COL=${PHENO_COL}"
log "COVAR_COLS=${COVAR_COLS}"
log "APPLY_RINT=${APPLY_RINT}"
log "STEP2_BLOCK_SIZE=${STEP2_BLOCK_SIZE}"
df -h /mnt/data | sed 's/^/  /'

keep_samples=$(wc -l < "${KEEP}")
variants=$(grep -vc '^#' "${SCRATCH}/${CHROM}.pvar")
if [[ "${keep_samples}" -ne "${EXPECTED_KEEP_SAMPLES}" ]]; then
    log "ERROR: keep file has ${keep_samples} samples, expected ${EXPECTED_KEEP_SAMPLES}"
    exit 1
fi
if [[ "${variants}" -ne "${EXPECTED_VARIANTS}" ]]; then
    log "ERROR: ${CHROM} has ${variants} variants, expected ${EXPECTED_VARIANTS}"
    exit 1
fi

pred_list=$(find "${STEP1_DIR}" -maxdepth 1 -name '*_pred.list' | head -1)
if [[ -z "${pred_list}" || ! -s "${pred_list}" ]]; then
    log "ERROR: could not find Step 1 *_pred.list in ${STEP1_DIR}"
    find "${STEP1_DIR}" -maxdepth 1 -type f -print | sed 's/^/  /'
    exit 1
fi
pred_base=$(basename "${pred_list}")
pred_dir=$(dirname "${pred_list}")
pred_rewritten="${SCRATCH}/${pred_base%.list}.localized.list"
while read -r pred_pheno pred_path extra; do
    if [[ -z "${pred_pheno:-}" ]]; then
        continue
    fi
    if [[ -n "${extra:-}" ]]; then
        log "ERROR: unexpected extra columns in ${pred_list}: ${pred_pheno} ${pred_path} ${extra}"
        exit 1
    fi
    pred_file="${pred_dir}/$(basename "${pred_path}")"
    if [[ ! -s "${pred_file}" ]]; then
        log "ERROR: prediction file listed in ${pred_list} is missing after localization: ${pred_file}"
        find "${pred_dir}" -maxdepth 1 -type f -print | sed 's/^/  /'
        exit 1
    fi
    printf '%s %s\n' "${pred_pheno}" "${pred_file}" >> "${pred_rewritten}"
done < "${pred_list}"
if [[ ! -s "${pred_rewritten}" ]]; then
    log "ERROR: localized prediction list is empty"
    exit 1
fi

cmd=(
    "${REGENIE_BIN}"
    --step 2
    --pgen "${SCRATCH}/${CHROM}"
    --phenoFile "${PHEN}"
    --phenoCol "${PHENO_COL}"
    --covarFile "${COVAR}"
    --covarColList "${COVAR_COLS}"
    --keep "${KEEP}"
    --qt
    --bsize "${STEP2_BLOCK_SIZE}"
    --threads "$(nproc)"
    --pred "${pred_rewritten}"
    --gz
    --no-split
    --out "${OUTDIR}/${CHROM}_height"
)
if [[ "${APPLY_RINT}" == "1" ]]; then
    cmd+=(--apply-rint)
fi

log "running REGENIE Step 2 with localized prediction list ${pred_rewritten}"
"${cmd[@]}"

result_file=$(find "${OUTDIR}" -maxdepth 1 -type f \
    \( -name "${CHROM}_height.regenie.gz" -o -name "${CHROM}_height.regenie" \) |
    head -1)
if [[ -z "${result_file}" || ! -s "${result_file}" ]]; then
    log "ERROR: missing REGENIE Step 2 result file in ${OUTDIR}"
    ls -lh "${OUTDIR}" | sed 's/^/  /'
    exit 1
fi

if [[ "${result_file}" == *.gz ]]; then
    result_lines=$(gzip -cd "${result_file}" | wc -l)
else
    result_lines=$(wc -l < "${result_file}")
fi
tested_variants=$((result_lines > 0 ? result_lines - 1 : 0))

cp "${PARAMS}" "${OUTDIR}/${CHROM}_height.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'chrom\t%s\n' "${CHROM}"
    printf 'source_pfile\t%s\n' "gwas_genotypes/step2_wgs_pfiles/${CHROM}"
    printf 'source_variants\t%s\n' "${variants}"
    printf 'keep_samples\t%s\n' "${keep_samples}"
    printf 'pheno_col\t%s\n' "${PHENO_COL}"
    printf 'covar_cols\t%s\n' "${COVAR_COLS}"
    printf 'apply_rint\t%s\n' "${APPLY_RINT}"
    printf 'step2_block_size\t%s\n' "${STEP2_BLOCK_SIZE}"
    printf 'result_file\t%s\n' "$(basename "${result_file}")"
    printf 'result_lines\t%s\n' "${result_lines}"
    printf 'tested_variants\t%s\n' "${tested_variants}"
} > "${OUTDIR}/${CHROM}_height.summary.tsv"

log "=== done ${CHROM}: ${tested_variants} tested variants ==="
cat "${OUTDIR}/${CHROM}_height.summary.tsv" | sed 's/^/  /'
ls -lh "${OUTDIR}" | sed 's/^/  /'
