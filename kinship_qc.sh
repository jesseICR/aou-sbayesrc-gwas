#!/bin/bash
# kinship_qc.sh - Compare our KING kinship output against AoU relatedness.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${DX_KINSHIP_DIR:?DX_KINSHIP_DIR not set}"
: "${AOU_RELATEDNESS_FILE:?AOU_RELATEDNESS_FILE not set}"

KING_TABLE_FILTER="${KING_TABLE_FILTER:-0.035}"
KINSHIP_CLOSE_LOWER="${KINSHIP_CLOSE_LOWER:-0.1767}"

kin0="${DX_KINSHIP_DIR}/aou_hq_direct_rel.kin0"
king_summary="${DX_KINSHIP_DIR}/aou_hq_direct_rel.summary.tsv"
outdir="${DX_KINSHIP_DIR}/qc"
mkdir -p "${outdir}"

if [[ ! -s "${kin0}" || ! -s "${king_summary}" ]]; then
    echo "ERROR: missing KING outputs; run run_king_kinship.sh first" >&2
    exit 1
fi
if [[ ! -s "${AOU_RELATEDNESS_FILE}" ]]; then
    echo "ERROR: missing AoU relatedness file ${AOU_RELATEDNESS_FILE}" >&2
    exit 1
fi

desired_params="${outdir}/kinship_qc.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'ours_kin0\t%s\n' "${kin0}"
    printf 'ours_kin0_size\t%s\n' "$(stat -c%s "${kin0}")"
    printf 'aou_relatedness_file\t%s\n' "${AOU_RELATEDNESS_FILE}"
    printf 'aou_relatedness_file_size\t%s\n' "$(stat -c%s "${AOU_RELATEDNESS_FILE}")"
    printf 'king_table_filter\t%s\n' "${KING_TABLE_FILTER}"
    printf 'kinship_close_lower\t%s\n' "${KINSHIP_CLOSE_LOWER}"
} > "${desired_params}"

params="${outdir}/kinship_qc.params.tsv"
required_outputs=(
    "${outdir}/kinship_comparison_summary.tsv"
    "${outdir}/kinship_comparison_summary.txt"
    "${outdir}/kinship_comparison_pairs.tsv"
    "${outdir}/kinship_comparison_plots.png"
    "${outdir}/kinship_bland_altman.png"
)

outputs_present=1
for output in "${required_outputs[@]}"; do
    if [[ ! -s "${output}" ]]; then
        outputs_present=0
        break
    fi
done

if [[ "${outputs_present}" -eq 1 ]] && diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
    overlap=$(awk -F'\t' '$1 == "overlapping_pairs" {print $2; exit}' "${outdir}/kinship_comparison_summary.tsv")
    echo "  Kinship QC already exists (${overlap:-unknown} overlapping pairs) — skipping"
    exit 0
fi

echo "  Comparing our KING output to AoU relatedness ..."
MPLBACKEND=Agg python3 "${SCRIPT_DIR}/kinship_qc.py" \
    --ours-kin0 "${kin0}" \
    --aou-relatedness "${AOU_RELATEDNESS_FILE}" \
    --output-dir "${outdir}" \
    --king-table-filter "${KING_TABLE_FILTER}"
cp "${desired_params}" "${params}"

echo "  Kinship QC summary:"
awk -F'\t' '
    $1 ~ /^(ours_total_pairs|aou_total_pairs|overlapping_pairs|ours_only_pairs|aou_only_pairs|pearson_r|mean_abs_diff|median_abs_diff)$/ {
        printf "    %s = %s\n", $1, $2
    }
' "${outdir}/kinship_comparison_summary.tsv"
