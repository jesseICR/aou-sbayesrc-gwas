#!/bin/bash
# compare_aou_ancestry.sh - Compare AoU-provided ancestry fractions to ours.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${AOU_ANCESTRY_PRED_FILE:?AOU_ANCESTRY_PRED_FILE not set}"
: "${AOU_ADMIXTURE_Q_FILE:?AOU_ADMIXTURE_Q_FILE not set}"
: "${DX_STATGEN_DIR:?DX_STATGEN_DIR not set}"
: "${DX_AOU_VS_OURS_DIR:?DX_AOU_VS_OURS_DIR not set}"
: "${DX_EUROPEANS_DIR:?DX_EUROPEANS_DIR not set}"
: "${LOCAL_ANCESTRY_COMPARE_DIR:?LOCAL_ANCESTRY_COMPARE_DIR not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"

OURS_EUR_MIN="${OURS_EUR_MIN:-0.8}"
OURS_AFR_MAX="${OURS_AFR_MAX:-0.1}"
OURS_AMR_MAX="${OURS_AMR_MAX:-0.1}"
OURS_EAS_MAX="${OURS_EAS_MAX:-0.1}"
OURS_OCE_MAX="${OURS_OCE_MAX:-0.1}"
AOU_MID_THRESHOLDS="${AOU_MID_THRESHOLDS:-0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2}"

mkdir -p "${DX_AOU_VS_OURS_DIR}" "${DX_AOU_VS_OURS_DIR}/plots" "${DX_EUROPEANS_DIR}" "${LOCAL_ANCESTRY_COMPARE_DIR}"

ours_admix="${DX_STATGEN_DIR}/aou_admixture_k${ADMIXTURE_K}.tsv"
for path in "${AOU_ANCESTRY_PRED_FILE}" "${AOU_ADMIXTURE_Q_FILE}" "${ours_admix}"; do
    if [[ ! -s "${path}" ]]; then
        echo "ERROR: missing required ancestry input ${path}" >&2
        exit 1
    fi
done

desired_params="${LOCAL_ANCESTRY_COMPARE_DIR}/aou_vs_ours.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'admixture_k\t%s\n' "${ADMIXTURE_K}"
    printf 'ours_admixture_file\t%s\n' "${ours_admix}"
    printf 'ours_admixture_size\t%s\n' "$(stat -c%s "${ours_admix}")"
    printf 'aou_admixture_q_file\t%s\n' "${AOU_ADMIXTURE_Q_FILE}"
    printf 'aou_admixture_q_size\t%s\n' "$(stat -c%s "${AOU_ADMIXTURE_Q_FILE}")"
    printf 'aou_ancestry_pred_file\t%s\n' "${AOU_ANCESTRY_PRED_FILE}"
    printf 'aou_ancestry_pred_size\t%s\n' "$(stat -c%s "${AOU_ANCESTRY_PRED_FILE}")"
    printf 'ours_eur_min\t%s\n' "${OURS_EUR_MIN}"
    printf 'ours_afr_max\t%s\n' "${OURS_AFR_MAX}"
    printf 'ours_amr_max\t%s\n' "${OURS_AMR_MAX}"
    printf 'ours_eas_max\t%s\n' "${OURS_EAS_MAX}"
    printf 'ours_oce_max\t%s\n' "${OURS_OCE_MAX}"
    printf 'aou_mid_thresholds\t%s\n' "${AOU_MID_THRESHOLDS}"
} > "${desired_params}"

params="${DX_AOU_VS_OURS_DIR}/aou_vs_ours.params.tsv"
summary="${DX_AOU_VS_OURS_DIR}/aou_vs_ours.summary.tsv"
eur_keep="${DX_EUROPEANS_DIR}/classified_european_iids.txt"
required_outputs=(
    "${summary}"
    "${params}"
    "${eur_keep}"
    "${DX_AOU_VS_OURS_DIR}/component_pair_metrics.tsv"
    "${DX_AOU_VS_OURS_DIR}/discordant_set_component_summary.tsv"
    "${DX_AOU_VS_OURS_DIR}/european_set_overlap_summary.tsv"
    "${DX_AOU_VS_OURS_DIR}/aou_mid_threshold_summary.tsv"
    "${DX_AOU_VS_OURS_DIR}/plots/scatter_aou_vs_ours_European.png"
    "${DX_AOU_VS_OURS_DIR}/plots/aou_mid_threshold_ours_means.png"
    "${DX_AOU_VS_OURS_DIR}/plots/discordant_mean_composition.png"
    "${DX_AOU_VS_OURS_DIR}/plots/discordant_component_boxplots.png"
    "${DX_AOU_VS_OURS_DIR}/plots/discordant_european_fraction_histograms.png"
)

outputs_present=1
for output in "${required_outputs[@]}"; do
    if [[ ! -s "${output}" ]]; then
        outputs_present=0
        break
    fi
done

if [[ "${outputs_present}" -eq 1 ]] && diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
    samples=$(awk -F'\t' '$1 == "samples_in_all_sources" {print $2; exit}' "${summary}")
    ours_eur=$(awk -F'\t' '$1 == "ours_european_count" {print $2; exit}' "${summary}")
    echo "  AoU-vs-ours ancestry comparison already exists (${samples:-unknown} samples; ${ours_eur:-unknown} ours-European) — skipping"
    exit 0
fi

echo "  Building AoU-vs-ours ancestry comparison outputs ..."
MPLBACKEND=Agg python3 "${SCRIPT_DIR}/compare_aou_ancestry.py" \
    --ours-admixture "${ours_admix}" \
    --aou-admixture-q "${AOU_ADMIXTURE_Q_FILE}" \
    --aou-ancestry-pred "${AOU_ANCESTRY_PRED_FILE}" \
    --output-dir "${DX_AOU_VS_OURS_DIR}" \
    --europeans-dir "${DX_EUROPEANS_DIR}" \
    --eur-min "${OURS_EUR_MIN}" \
    --afr-max "${OURS_AFR_MAX}" \
    --amr-max "${OURS_AMR_MAX}" \
    --eas-max "${OURS_EAS_MAX}" \
    --oce-max "${OURS_OCE_MAX}" \
    --mid-thresholds "${AOU_MID_THRESHOLDS}"

cp "${desired_params}" "${params}"

echo "  Wrote ancestry comparison to ${DX_AOU_VS_OURS_DIR}"
echo "  Wrote European keep-list to ${eur_keep}"
if [[ -s "${summary}" ]]; then
    echo "  Key counts:"
    awk -F'\t' '
        $1 ~ /^(samples_in_all_sources|aou_eur_pred_count|ours_european_count|european_both|aou_eur_only|ours_eur_only|neither_european)$/ {
            printf "    %s = %s\n", $1, $2
        }
    ' "${summary}"
fi
