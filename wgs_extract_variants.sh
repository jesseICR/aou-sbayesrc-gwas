#!/bin/bash
# wgs_extract_variants.sh — Submit per-chromosome SBayesRC extraction jobs via dsub.
#
# For each autosome (chr1..chr22), launches a Google Batch job that runs:
#   plink2 \
#     --bed/--bim/--fam <AoU acaf_threshold/plink_bed/chr{N}.{bed,bim,fam}> \
#     --extract <chr{N}.extract.txt>     # SBayesRC chr:pos:REF:ALT IDs
#     --update-name <chr{N}.idmap.txt>   # remap IDs to rsids in one shot
#     --no-pheno --make-pgen --out chr{N}
# and uploads the resulting chr{N}.{pgen,pvar,psam} + chr{N}.summary.tsv to
# ${DX_WGS_PFILE_DIR}.
#
# Idempotent: skips a chromosome if chr{N}.pgen already exists at the
# destination (cheap metadata-only `gcloud storage ls`).
#
# Expects env vars (set by get_genotypes.sh):
#   GOOGLE_PROJECT          billing + Batch project
#   DSUB_REGION             e.g. us-central1
#   PLINK2_IMAGE            e.g. quay.io/biocontainers/plink2:...
#   DSUB_MACHINE_TYPE       e.g. n2-standard-8
#   AOU_PLINK_BED_DIR       gs://vwb-aou-datasets-controlled/v8/wgs/.../plink_bed
#   DX_SBAYESRC_ID_DIR      gs://...  (dest for chr{N}.{extract,idmap}.txt)
#   DX_WGS_PFILE_DIR        gs://...  (dest for chr{N}.{pgen,pvar,psam})
#   DX_LOGS_DIR             gs://...  (dsub --logging prefix)
#   SCRIPT_DIR              local script dir (logs/dsub/ written here)
#
# Optional smoke-test override:
#   SBAYESRC_TEST_CHROM     If set (e.g. "22"), submit only that chromosome
#                           and skip the others. Useful for first-run
#                           validation before launching all 22 in parallel.

set -euo pipefail

mkdir -p "${SCRIPT_DIR}/logs/dsub"

submit_extract_job() {
    local n=$1

    local disk
    case "${n}" in
        1|2)         disk=800 ;;
        3|4|5|6)     disk=600 ;;
        *)           disk=400 ;;
    esac

    # Use single-quoted heredoc so dsub gets the literal command and expands
    # the env vars (BED, BIM, ..., CHR) on the worker.
    local cmd
    cmd=$(cat <<'EOF'
set -euo pipefail

REQUESTED=$(wc -l < "${EXTRACT}")
BIM_TOTAL=$(wc -l < "${BIM}")
SAMPLES=$(wc -l < "${FAM}")
echo "[chr${CHR}] requested SBayesRC variants: ${REQUESTED}"
echo "[chr${CHR}] AoU bim variants: ${BIM_TOTAL}"
echo "[chr${CHR}] AoU samples: ${SAMPLES}"

mkdir -p "${OUT_DIR}"
PREFIX="${OUT_DIR}/chr${CHR}"

plink2 \
    --bed "${BED}" \
    --bim "${BIM}" \
    --fam "${FAM}" \
    --extract "${EXTRACT}" \
    --update-name "${IDMAP}" \
    --no-pheno \
    --make-pgen \
    --out "${PREFIX}"

EXTRACTED=$(grep -vc '^#' "${PREFIX}.pvar")
OUT_SAMPLES=$(grep -vc '^#' "${PREFIX}.psam")
MISSING=$(( REQUESTED - EXTRACTED ))
echo "[chr${CHR}] extracted variants: ${EXTRACTED}"
echo "[chr${CHR}] output samples: ${OUT_SAMPLES}"
echo "[chr${CHR}] SBayesRC variants not found in AoU: ${MISSING}"

# One-line summary uploaded alongside the pfile triplet.
{
    printf 'chrom\trequested\tbim_total\tsamples\textracted\tout_samples\tmissing\n'
    printf 'chr%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${CHR}" "${REQUESTED}" "${BIM_TOTAL}" "${SAMPLES}" \
        "${EXTRACTED}" "${OUT_SAMPLES}" "${MISSING}"
} > "${PREFIX}.summary.tsv"
EOF
)

    dsub \
        --provider google-batch \
        --project "${GOOGLE_PROJECT}" \
        --regions "${DSUB_REGION}" \
        --user-project "${GOOGLE_PROJECT}" \
        --logging "${DX_LOGS_DIR}/wgs_extract/chr${n}/" \
        --image "${PLINK2_IMAGE}" \
        --machine-type "${DSUB_MACHINE_TYPE}" \
        --boot-disk-size 50 \
        --disk-size "${disk}" \
        --env CHR="${n}" \
        --input BED="${AOU_PLINK_BED_DIR}/chr${n}.bed" \
        --input BIM="${AOU_PLINK_BED_DIR}/chr${n}.bim" \
        --input FAM="${AOU_PLINK_BED_DIR}/chr${n}.fam" \
        --input EXTRACT="${DX_SBAYESRC_ID_DIR}/chr${n}.extract.txt" \
        --input IDMAP="${DX_SBAYESRC_ID_DIR}/chr${n}.idmap.txt" \
        --output-recursive OUT_DIR="${DX_WGS_PFILE_DIR}/" \
        --name "wgs_extract_chr${n}" \
        --command "${cmd}" \
        --wait
}

if [[ -n "${SBAYESRC_TEST_CHROM:-}" ]]; then
    chrom_list="${SBAYESRC_TEST_CHROM}"
    echo "SBAYESRC_TEST_CHROM=${SBAYESRC_TEST_CHROM} — submitting only this chromosome."
else
    chrom_list=$(seq 1 22)
fi

submitted=0
skipped=0
pids=()

for chrom in ${chrom_list}; do
    if gcloud storage ls --billing-project="${GOOGLE_PROJECT}" \
            "${DX_WGS_PFILE_DIR}/chr${chrom}.pgen" &>/dev/null; then
        echo "chr${chrom}: skipping — pfiles already exist"
        skipped=$((skipped + 1))
        continue
    fi
    echo "chr${chrom}: submitting dsub job (logs/dsub/chr${chrom}.log)"
    submit_extract_job "${chrom}" \
        > "${SCRIPT_DIR}/logs/dsub/chr${chrom}.log" 2>&1 &
    pids+=("$!")
    submitted=$((submitted + 1))
done

echo ""
echo "Submitted: ${submitted}, Skipped (already exist): ${skipped}"

if [[ ${#pids[@]} -gt 0 ]]; then
    echo "Waiting for ${#pids[@]} dsub jobs ..."
    fail=0
    for pid in "${pids[@]}"; do
        if ! wait "${pid}"; then
            fail=$((fail + 1))
        fi
    done
    if [[ ${fail} -gt 0 ]]; then
        echo "ERROR: ${fail} dsub job(s) failed. Inspect logs/dsub/chr*.log."
        exit 1
    fi
    echo "All extraction jobs complete."
fi

# Print combined summary by streaming the per-chromosome summary.tsv files
# back from the bucket. These are tiny (1 data line each), so streaming is
# both cheap and consistent with the no-bulk-download rule.
echo ""
echo "=== SBayesRC extraction summary (all chromosomes) ==="
combined="${SCRIPT_DIR}/logs/sbayesrc_extract_summary.tsv"
{
    header_written=0
    for chrom in ${chrom_list}; do
        uri="${DX_WGS_PFILE_DIR}/chr${chrom}.summary.tsv"
        if ! gcloud storage ls --billing-project="${GOOGLE_PROJECT}" \
                "${uri}" &>/dev/null; then
            continue
        fi
        if [[ ${header_written} -eq 0 ]]; then
            gcloud storage cat --billing-project="${GOOGLE_PROJECT}" "${uri}" \
                | head -1
            header_written=1
        fi
        gcloud storage cat --billing-project="${GOOGLE_PROJECT}" "${uri}" \
            | tail -n +2
    done
} | tee "${combined}"

echo ""
echo "Summary written to ${combined}"
