#!/bin/bash
# fit_project_pca.sh - Fit PCA on unrelated Europeans and project onto all samples.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_URI:?DX_HQ_DIRECT_BFILE_URI not set}"
: "${DX_PCA_EUR_DIR:?DX_PCA_EUR_DIR not set}"
: "${DX_PCA_EUR_URI:?DX_PCA_EUR_URI not set}"
: "${LOCAL_PCA_QC_DIR:?LOCAL_PCA_QC_DIR not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

PCA_NPCS="${PCA_NPCS:-20}"
PCA_SEED="${PCA_SEED:-0}"
PCA_PROJECT_DSUB_MIN_CORES="${PCA_PROJECT_DSUB_MIN_CORES:-16}"
PCA_PROJECT_DSUB_MIN_RAM="${PCA_PROJECT_DSUB_MIN_RAM:-64}"
PCA_PROJECT_DSUB_DISK_SIZE="${PCA_PROJECT_DSUB_DISK_SIZE:-300}"
PCA_PROJECT_DSUB_DISK_TYPE="${PCA_PROJECT_DSUB_DISK_TYPE:-pd-ssd}"

mkdir -p "${LOCAL_PCA_QC_DIR}" "${DX_PCA_EUR_DIR}" "${DX_PCA_EUR_DIR}/scrap"

pca_prefix="${DX_PCA_EUR_DIR}/pca_ready"
pca_uri_prefix="${DX_PCA_EUR_URI}/pca_ready"
source_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
source_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"

for ext in bed bim fam; do
    if [[ ! -s "${pca_prefix}.${ext}" ]]; then
        echo "ERROR: missing PCA-ready bfile input ${pca_prefix}.${ext}" >&2
        exit 1
    fi
    if [[ ! -s "${source_prefix}.${ext}" ]]; then
        echo "ERROR: missing HQ direct bfile input ${source_prefix}.${ext}" >&2
        exit 1
    fi
done

pca_variants=$(wc -l < "${pca_prefix}.bim")
pca_samples=$(wc -l < "${pca_prefix}.fam")
source_variants=$(wc -l < "${source_prefix}.bim")
source_samples=$(wc -l < "${source_prefix}.fam")
if [[ "${pca_variants}" -le 0 || "${pca_samples}" -le 0 ]]; then
    echo "ERROR: PCA-ready bfile has invalid counts: ${pca_variants} variants, ${pca_samples} samples" >&2
    exit 1
fi
if [[ "${source_samples}" -le 0 ]]; then
    echo "ERROR: HQ direct bfile has no samples" >&2
    exit 1
fi

desired_params="${LOCAL_PCA_QC_DIR}/fit_project_pca.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'pca_bfile\t%s\n' "pca_eur/pca_ready"
    printf 'pca_variants\t%s\n' "${pca_variants}"
    printf 'pca_samples\t%s\n' "${pca_samples}"
    printf 'pca_bed_size\t%s\n' "$(stat -c%s "${pca_prefix}.bed")"
    printf 'pca_bim_size\t%s\n' "$(stat -c%s "${pca_prefix}.bim")"
    printf 'pca_fam_size\t%s\n' "$(stat -c%s "${pca_prefix}.fam")"
    printf 'projection_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'projection_variants\t%s\n' "${source_variants}"
    printf 'projection_samples\t%s\n' "${source_samples}"
    printf 'projection_bed_size\t%s\n' "$(stat -c%s "${source_prefix}.bed")"
    printf 'projection_bim_size\t%s\n' "$(stat -c%s "${source_prefix}.bim")"
    printf 'projection_fam_size\t%s\n' "$(stat -c%s "${source_prefix}.fam")"
    printf 'pca_npcs\t%s\n' "${PCA_NPCS}"
    printf 'pca_seed\t%s\n' "${PCA_SEED}"
    printf 'plink2_version\t%s\n' "$("${PLINK2}" --version 2>&1 | head -1)"
} > "${desired_params}"

params="${DX_PCA_EUR_DIR}/fit_project_pca.params.tsv"
summary="${DX_PCA_EUR_DIR}/fit_project_pca.summary.tsv"
projected="${DX_PCA_EUR_DIR}/aou_projected.sscore"
eigenvec="${DX_PCA_EUR_DIR}/aou_pcs.eigenvec"
eigenvec_allele="${DX_PCA_EUR_DIR}/aou_pcs.eigenvec.allele"
eigenval="${DX_PCA_EUR_DIR}/aou_pcs.eigenval"
acount="${DX_PCA_EUR_DIR}/pca_eur_counts.acount"

if [[ -s "${projected}" && -s "${eigenvec}" && -s "${eigenvec_allele}" &&
      -s "${eigenval}" && -s "${acount}" && -s "${summary}" && -s "${params}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected_fit_variants=$(awk -F'\t' '$1 == "fit_variants" {print $2; exit}' "${summary}")
        expected_fit_samples=$(awk -F'\t' '$1 == "fit_samples" {print $2; exit}' "${summary}")
        expected_projected_samples=$(awk -F'\t' '$1 == "projected_samples" {print $2; exit}' "${summary}")
        expected_pcs=$(awk -F'\t' '$1 == "pca_npcs" {print $2; exit}' "${summary}")
        observed_weight_ids=$(awk 'NR > 1 {ids[$2] = 1} END {print length(ids)}' "${eigenvec_allele}")
        observed_fit_samples=$(( $(wc -l < "${eigenvec}") - 1 ))
        observed_projected=$(( $(wc -l < "${projected}") - 1 ))
        observed_pcs=$(wc -l < "${eigenval}")
        if [[ -n "${expected_fit_variants}" && -n "${expected_fit_samples}" &&
              -n "${expected_projected_samples}" && -n "${expected_pcs}" &&
              "${observed_weight_ids}" -eq "${expected_fit_variants}" &&
              "${observed_fit_samples}" -eq "${expected_fit_samples}" &&
              "${observed_projected}" -eq "${expected_projected_samples}" &&
              "${observed_pcs}" -eq "${expected_pcs}" ]]; then
            echo "  PCA fit/projection already exists (${observed_weight_ids} variants, ${observed_fit_samples} fit samples, ${observed_projected} projected samples) — skipping"
            exit 0
        fi
    fi
    echo "  PCA fit/projection outputs exist but params/counts do not match — rebuilding"
fi

echo "  Staging plink2 and PCA fit/projection params ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
export DSUB_PCA_QC_PREP_URI="${DSUB_PCA_QC_PREP_URI:-${DX_PCA_EUR_URI}/_prep}"
gcloud storage cp "${desired_params}" "${DSUB_PCA_QC_PREP_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

pca_log="${SCRIPT_DIR}/logs/dsub_fit_project_pca_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${pca_log}")"

echo "  Submitting dsub job to fit ${PCA_NPCS} PCs and project ${source_samples} samples ..."
dsub \
    --provider "${DSUB_PROVIDER}" \
    --project "${GOOGLE_PROJECT}" \
    --regions "${DSUB_REGION}" \
    --service-account "${DSUB_PET_SA}" \
    --use-private-address \
    --network "${DSUB_NETWORK}" \
    --subnetwork "${DSUB_SUBNETWORK}" \
    --user-project "${GOOGLE_PROJECT}" \
    --logging "${DSUB_LOG_URI}" \
    --name "sbayesrc-fit-project-pca" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_fit_project_pca_worker.sh" \
    --env PCA_NPCS="${PCA_NPCS}" \
    --env PCA_SEED="${PCA_SEED}" \
    --env EXPECTED_PCA_VARIANTS="${pca_variants}" \
    --env EXPECTED_PCA_SAMPLES="${pca_samples}" \
    --env EXPECTED_PROJECTED_SAMPLES="${source_samples}" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input PCA_BED="${pca_uri_prefix}.bed" \
    --input PCA_BIM="${pca_uri_prefix}.bim" \
    --input PCA_FAM="${pca_uri_prefix}.fam" \
    --input DIRECT_BED="${source_uri_prefix}.bed" \
    --input DIRECT_BIM="${source_uri_prefix}.bim" \
    --input DIRECT_FAM="${source_uri_prefix}.fam" \
    --input PARAMS="${DSUB_PCA_QC_PREP_URI}/fit_project_pca.desired_params.tsv" \
    --output-recursive OUTDIR="${DX_PCA_EUR_URI}/" \
    --min-cores "${PCA_PROJECT_DSUB_MIN_CORES}" \
    --min-ram "${PCA_PROJECT_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${PCA_PROJECT_DSUB_DISK_SIZE}" \
    --disk-type "${PCA_PROJECT_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${pca_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub PCA fit/projection job returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

echo "  Verifying PCA fit/projection outputs ..."
for attempt in $(seq 1 60); do
    if [[ -s "${projected}" && -s "${eigenvec}" && -s "${eigenvec_allele}" &&
          -s "${eigenval}" && -s "${acount}" && -s "${summary}" && -s "${params}" ]]; then
        expected_fit_variants=$(awk -F'\t' '$1 == "fit_variants" {print $2; exit}' "${summary}")
        expected_fit_samples=$(awk -F'\t' '$1 == "fit_samples" {print $2; exit}' "${summary}")
        expected_projected_samples=$(awk -F'\t' '$1 == "projected_samples" {print $2; exit}' "${summary}")
        expected_pcs=$(awk -F'\t' '$1 == "pca_npcs" {print $2; exit}' "${summary}")
        observed_weight_ids=$(awk 'NR > 1 {ids[$2] = 1} END {print length(ids)}' "${eigenvec_allele}")
        observed_fit_samples=$(( $(wc -l < "${eigenvec}") - 1 ))
        observed_projected=$(( $(wc -l < "${projected}") - 1 ))
        observed_pcs=$(wc -l < "${eigenval}")
        if [[ "${observed_weight_ids}" -eq "${expected_fit_variants}" &&
              "${observed_fit_samples}" -eq "${expected_fit_samples}" &&
              "${observed_projected}" -eq "${expected_projected_samples}" &&
              "${observed_pcs}" -eq "${expected_pcs}" ]]; then
            echo "  Done: ${expected_fit_variants} variants, ${expected_fit_samples} fit samples, ${expected_projected_samples} projected samples"
            echo "  PCA fit/projection summary:"
            sed 's/^/    /' "${summary}"
            exit 0
        fi
        echo "  Outputs visible but counts do not match yet; waiting ..."
    else
        echo "  Waiting for PCA fit/projection outputs to appear (${attempt}/60) ..."
    fi
    sleep 10
done

echo "ERROR: PCA fit/projection outputs were not visible with expected counts after dsub completed" >&2
exit 1
