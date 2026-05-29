#!/bin/bash
# select_pca_europeans.sh - Select unrelated European IIDs for fitting PCA.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_EUROPEANS_DIR:?DX_EUROPEANS_DIR not set}"
: "${DX_KINSHIP_DIR:?DX_KINSHIP_DIR not set}"
: "${PLINK2:?PLINK2 not set}"

PCA_KINSHIP_THRESHOLD="${PCA_KINSHIP_THRESHOLD:-0.0441941}"
PCA_SEED_RELATIONSHIPS="${PCA_SEED_RELATIONSHIPS:-sibling,identical}"

mkdir -p "${DX_PCA_EUR_DIR}" "${DX_PCA_EUR_DIR}/scrap"

europeans="${DX_EUROPEANS_DIR}/classified_european_iids.txt"
close_relations="${DX_KINSHIP_DIR}/close_relations.csv"
kin0="${DX_KINSHIP_DIR}/aou_hq_direct_rel.kin0"

for input in "${europeans}" "${close_relations}" "${kin0}"; do
    if [[ ! -s "${input}" ]]; then
        echo "ERROR: missing PCA European selector input ${input}" >&2
        exit 1
    fi
done

desired_params="${DX_PCA_EUR_DIR}/select_pca_europeans.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'europeans_file\t%s\n' "${europeans}"
    printf 'europeans_file_size\t%s\n' "$(stat -c%s "${europeans}")"
    printf 'close_relations_file\t%s\n' "${close_relations}"
    printf 'close_relations_file_size\t%s\n' "$(stat -c%s "${close_relations}")"
    printf 'kin0_file\t%s\n' "${kin0}"
    printf 'kin0_file_size\t%s\n' "$(stat -c%s "${kin0}")"
    printf 'pca_kinship_threshold\t%s\n' "${PCA_KINSHIP_THRESHOLD}"
    printf 'pca_seed_relationships\t%s\n' "${PCA_SEED_RELATIONSHIPS}"
    printf 'plink2_version\t%s\n' "$("${PLINK2}" --version 2>&1 | head -1)"
} > "${desired_params}"

params="${DX_PCA_EUR_DIR}/select_pca_europeans.params.tsv"
summary="${DX_PCA_EUR_DIR}/select_pca_europeans.summary.tsv"
fit_iids="${DX_PCA_EUR_DIR}/fit_pca_iids.txt"
if [[ -s "${fit_iids}" && -s "${summary}" && -s "${params}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' '$1 == "fit_pca_iids" {print $2; exit}' "${summary}")
        observed=$(wc -l < "${fit_iids}")
        if [[ -n "${expected}" && "${observed}" -eq "${expected}" ]]; then
            echo "  PCA European fitting IIDs already exist (${observed} samples) — skipping"
            exit 0
        fi
    fi
    echo "  PCA European selector outputs exist but params/counts do not match — rebuilding"
fi

echo "  Selecting unrelated European IIDs for PCA fitting ..."
python3 "${SCRIPT_DIR}/select_pca_europeans.py" \
    --europeans "${europeans}" \
    --close-relations "${close_relations}" \
    --kin0 "${kin0}" \
    --output-dir "${DX_PCA_EUR_DIR}" \
    --plink2 "${PLINK2}" \
    --kinship-threshold "${PCA_KINSHIP_THRESHOLD}" \
    --seed-relationships "${PCA_SEED_RELATIONSHIPS}"
cp "${desired_params}" "${params}"

echo "  PCA European selector summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
