#!/bin/bash
# classify_relations.sh - Classify close relationships from KING kinship output.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${DX_KINSHIP_DIR:?DX_KINSHIP_DIR not set}"

KINSHIP_CLOSE_LOWER="${KINSHIP_CLOSE_LOWER:-0.1767}"
KINSHIP_FIRST_DEGREE_UPPER="${KINSHIP_FIRST_DEGREE_UPPER:-0.3535}"
KINSHIP_IBS0_CUTOFF="${KINSHIP_IBS0_CUTOFF:-0.0012}"

kin0="${DX_KINSHIP_DIR}/aou_hq_direct_rel.kin0"
outdir="${DX_KINSHIP_DIR}"
if [[ ! -s "${kin0}" ]]; then
    echo "ERROR: missing KING output ${kin0}" >&2
    exit 1
fi

desired_params="${outdir}/close_relations.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'kin0\t%s\n' "${kin0}"
    printf 'kin0_size\t%s\n' "$(stat -c%s "${kin0}")"
    printf 'kinship_close_lower\t%s\n' "${KINSHIP_CLOSE_LOWER}"
    printf 'kinship_first_degree_upper\t%s\n' "${KINSHIP_FIRST_DEGREE_UPPER}"
    printf 'kinship_ibs0_cutoff\t%s\n' "${KINSHIP_IBS0_CUTOFF}"
    printf 'age_gap_filter\t%s\n' "not_applied"
} > "${desired_params}"

params="${outdir}/close_relations.params.tsv"
close_csv="${outdir}/close_relations.csv"
summary="${outdir}/close_relations.summary.tsv"
if [[ -s "${close_csv}" && -s "${summary}" && -s "${params}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        total=$(awk -F'\t' '$1 == "total_close_relationships" {print $2; exit}' "${summary}")
        echo "  Close relationship classification already exists (${total:-unknown} close relationships) — skipping"
        exit 0
    fi
    echo "  Close relationship outputs exist but params do not match — rebuilding"
fi

echo "  Classifying close relationships from KING output ..."
python3 "${SCRIPT_DIR}/classify_relations.py" \
    --kin0 "${kin0}" \
    --output-dir "${outdir}" \
    --close-lower "${KINSHIP_CLOSE_LOWER}" \
    --first-degree-upper "${KINSHIP_FIRST_DEGREE_UPPER}" \
    --ibs0-cutoff "${KINSHIP_IBS0_CUTOFF}"
cp "${desired_params}" "${params}"

echo "  Close relationship summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
