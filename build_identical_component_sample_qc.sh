#!/bin/bash
# build_identical_component_sample_qc.sh - Exclude anomalous identical-genotype components.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${DX_KINSHIP_DIR:?DX_KINSHIP_DIR not set}"
: "${DX_SAMPLE_QC_DIR:?DX_SAMPLE_QC_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"

IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE="${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE:-3}"

close_relations="${DX_KINSHIP_DIR}/close_relations.csv"
fam="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.fam"
outdir="${DX_SAMPLE_QC_DIR}"
mkdir -p "${outdir}"

for f in "${close_relations}" "${fam}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing required input ${f}" >&2
        exit 1
    fi
done

desired_params="${outdir}/identical_component_sample_qc.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'close_relations\t%s\n' "${close_relations}"
    printf 'close_relations_size\t%s\n' "$(stat -c%s "${close_relations}")"
    printf 'fam\t%s\n' "${fam}"
    printf 'fam_size\t%s\n' "$(stat -c%s "${fam}")"
    printf 'exclude_min_component_size\t%s\n' "${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}"
    printf 'build_identical_component_sample_qc_py_sha256\t%s\n' "$(sha256sum "${SCRIPT_DIR}/build_identical_component_sample_qc.py" | awk '{print $1}')"
} > "${desired_params}"

params="${outdir}/identical_component_sample_qc.params.tsv"
summary="${outdir}/identical_component_sample_qc.summary.tsv"
members="${outdir}/identical_components.tsv"
component_summary="${outdir}/identical_component_summary.tsv"
exclude_iids="${outdir}/exclude_identical_component_size_ge${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}_iids.txt"
log_file="${outdir}/identical_component_sample_qc.log"

if [[ -s "${params}" && -s "${summary}" && -s "${members}" &&
      -s "${component_summary}" && -s "${exclude_iids}" && -s "${log_file}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' -v k="iids_in_components_size_ge_${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}" '$1 == k {gsub(/\r/, "", $2); print $2; exit}' "${summary}")
        observed=$(wc -l < "${exclude_iids}")
        if [[ -n "${expected}" && "${observed}" == "${expected}" ]]; then
            echo "  Identical-component sample QC already exists (${observed} excluded IIDs) — skipping"
            exit 0
        fi
    fi
    echo "  Identical-component sample QC outputs exist but params/counts do not match — rebuilding"
fi

echo "  Building identical-genotype component sample QC ..."
python3 "${SCRIPT_DIR}/build_identical_component_sample_qc.py" \
    --close-relations "${close_relations}" \
    --fam "${fam}" \
    --out-dir "${outdir}" \
    --exclude-min-component-size "${IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE}"
cp "${desired_params}" "${params}"

echo "  Identical-component sample QC summary:"
awk -F'\t' 'NR > 1 {gsub(/\r/, "", $2); printf "    %s = %s\n", $1, $2}' "${summary}"
