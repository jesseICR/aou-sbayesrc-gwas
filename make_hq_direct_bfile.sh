#!/bin/bash
# make_hq_direct_bfile.sh - Build high-quality direct-SNP bfile for REGENIE step 1.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${LOCAL_DIRECT_SNPS_FILE:?LOCAL_DIRECT_SNPS_FILE not set}"
: "${LOCAL_HQ_DIRECT_DIR:?LOCAL_HQ_DIRECT_DIR not set}"
: "${SBAYESRC_LIFTOVER_FILE:?SBAYESRC_LIFTOVER_FILE not set}"
: "${AOU_ANCESTRY_PRED_FILE:?AOU_ANCESTRY_PRED_FILE not set}"
: "${DX_DIRECT_BFILE_DIR:?DX_DIRECT_BFILE_DIR not set}"
: "${DX_DIRECT_BFILE_URI:?DX_DIRECT_BFILE_URI not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_URI:?DX_HQ_DIRECT_BFILE_URI not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

mkdir -p "${LOCAL_HQ_DIRECT_DIR}" "${DX_HQ_DIRECT_BFILE_DIR}"

HQ_AF_DIFF_MAX="${HQ_AF_DIFF_MAX:-0.04}"
HQ_EUR_MAF_MIN="${HQ_EUR_MAF_MIN:-0.007}"
HQ_EUR_MISSING_MAX="${HQ_EUR_MISSING_MAX:-0.05}"
HQ_DIRECT_DSUB_MIN_CORES="${HQ_DIRECT_DSUB_MIN_CORES:-8}"
HQ_DIRECT_DSUB_MIN_RAM="${HQ_DIRECT_DSUB_MIN_RAM:-32}"
HQ_DIRECT_DSUB_DISK_SIZE="${HQ_DIRECT_DSUB_DISK_SIZE:-300}"
HQ_DIRECT_DSUB_DISK_TYPE="${HQ_DIRECT_DSUB_DISK_TYPE:-pd-ssd}"

direct_prefix="${DX_DIRECT_BFILE_DIR}/chr1_22_merged"
direct_uri_prefix="${DX_DIRECT_BFILE_URI}/chr1_22_merged"
hq_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
hq_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"

for ext in bed bim fam; do
    if [[ ! -s "${direct_prefix}.${ext}" ]]; then
        echo "ERROR: missing direct bfile input ${direct_prefix}.${ext}" >&2
        exit 1
    fi
done
if [[ ! -s "${LOCAL_DIRECT_SNPS_FILE}" ]]; then
    echo "ERROR: missing original direct SNP list ${LOCAL_DIRECT_SNPS_FILE}" >&2
    exit 1
fi
if [[ ! -s "${SBAYESRC_LIFTOVER_FILE}" ]]; then
    echo "ERROR: missing SBayesRC liftover file ${SBAYESRC_LIFTOVER_FILE}" >&2
    exit 1
fi
if [[ ! -s "${AOU_ANCESTRY_PRED_FILE}" ]]; then
    echo "ERROR: missing AoU ancestry predictions ${AOU_ANCESTRY_PRED_FILE}" >&2
    exit 1
fi

requested_direct=$(wc -l < "${LOCAL_DIRECT_SNPS_FILE}")
direct_variants=$(wc -l < "${direct_prefix}.bim")

desired_params="${LOCAL_HQ_DIRECT_DIR}/chr1_22_merged_hq.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'af_diff_max\t%s\n' "${HQ_AF_DIFF_MAX}"
    printf 'aou_eur_maf_min\t%s\n' "${HQ_EUR_MAF_MIN}"
    printf 'aou_eur_missing_rate_max\t%s\n' "${HQ_EUR_MISSING_MAX}"
    printf 'requested_direct_snps\t%s\n' "${requested_direct}"
    printf 'direct_bfile_variants\t%s\n' "${direct_variants}"
    printf 'liftover_file_size\t%s\n' "$(stat -c%s "${SBAYESRC_LIFTOVER_FILE}")"
} > "${desired_params}"

summary="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.filter_summary.tsv"
params="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.params.tsv"
sample_summary="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq.sample_missingness_summary.tsv"
if [[ -s "${hq_prefix}.bed" && -s "${hq_prefix}.bim" && -s "${hq_prefix}.fam" && -s "${summary}" && -s "${params}" && -s "${sample_summary}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected_final=$(awk -F'\t' '$1 == "final_hq_direct_snps" {print $2; exit}' "${summary}")
        observed_final=$(wc -l < "${hq_prefix}.bim")
        if [[ -n "${expected_final}" && "${observed_final}" -eq "${expected_final}" ]]; then
            echo "  High-quality direct bfile already exists (${observed_final} variants) — skipping"
            exit 0
        fi
    fi
    echo "  High-quality direct bfile exists but params/counts do not match — rebuilding"
fi

echo "  Building AoU EUR keep-list from ancestry_pred == eur ..."
eur_keep="${LOCAL_HQ_DIRECT_DIR}/aou_eur.keep"
awk -F'\t' '
    NR == FNR {
        if (FNR > 1 && $2 == "eur") {
            eur[$1] = 1
        }
        next
    }
    ($2 in eur) {
        print $1, $2
    }
' "${AOU_ANCESTRY_PRED_FILE}" "${direct_prefix}.fam" > "${eur_keep}"
eur_samples=$(wc -l < "${eur_keep}")
if [[ "${eur_samples}" -le 0 ]]; then
    echo "ERROR: EUR keep-list is empty: ${eur_keep}" >&2
    exit 1
fi
echo "  EUR samples in direct bfile: ${eur_samples}"

echo "  Staging plink2 + EUR keep-list to workspace bucket ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
export DSUB_HQ_DIRECT_PREP_URI="${DSUB_HQ_DIRECT_PREP_URI:-${DX_HQ_DIRECT_BFILE_URI}/_prep}"
gcloud storage cp "${eur_keep}" "${DSUB_HQ_DIRECT_PREP_URI}/aou_eur.keep" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

metrics_afreq="${DX_HQ_DIRECT_BFILE_DIR}/eur_direct_qc.afreq"
metrics_vmiss="${DX_HQ_DIRECT_BFILE_DIR}/eur_direct_qc.vmiss"
metrics_expected_lines=$((direct_variants + 1))
if [[ -s "${metrics_afreq}" && -s "${metrics_vmiss}" ]]; then
    afreq_lines=$(wc -l < "${metrics_afreq}")
    vmiss_lines=$(wc -l < "${metrics_vmiss}")
else
    afreq_lines=0
    vmiss_lines=0
fi

if [[ "${afreq_lines}" -eq "${metrics_expected_lines}" && "${vmiss_lines}" -eq "${metrics_expected_lines}" ]]; then
    echo "  EUR direct SNP metrics already exist (${direct_variants} variants) — skipping metrics dsub"
else
    metrics_log="${SCRIPT_DIR}/logs/dsub_hq_direct_metrics_$(date +%Y%m%d_%H%M%S).dsub.out"
    mkdir -p "$(dirname "${metrics_log}")"
    echo "  Submitting dsub job for EUR direct SNP frequency/missingness metrics ..."
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
        --name "sbayesrc-hq-direct-metrics" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_hq_direct_metrics_worker.sh" \
        --input PLINK2="${DSUB_PLINK2_GS}" \
        --input BED="${direct_uri_prefix}.bed" \
        --input BIM="${direct_uri_prefix}.bim" \
        --input FAM="${direct_uri_prefix}.fam" \
        --input KEEP="${DSUB_HQ_DIRECT_PREP_URI}/aou_eur.keep" \
        --output-recursive OUTDIR="${DX_HQ_DIRECT_BFILE_URI}/" \
        --min-cores "${HQ_DIRECT_DSUB_MIN_CORES}" \
        --min-ram "${HQ_DIRECT_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${HQ_DIRECT_DSUB_DISK_SIZE}" \
        --disk-type "${HQ_DIRECT_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${metrics_log}"
    dsub_rc=${PIPESTATUS[0]}
    if [[ "${dsub_rc}" -ne 0 ]]; then
        echo "ERROR: dsub high-quality direct metrics job returned ${dsub_rc}" >&2
        echo "Check logs at ${DSUB_LOG_URI}/" >&2
        exit "${dsub_rc}"
    fi
fi

if [[ ! -s "${metrics_afreq}" || ! -s "${metrics_vmiss}" ]]; then
    echo "ERROR: missing high-quality direct metrics outputs" >&2
    exit 1
fi

echo "  Applying high-quality direct SNP filters ..."
python3 "${SCRIPT_DIR}/filter_hq_direct_snps.py" \
    --direct-snps "${LOCAL_DIRECT_SNPS_FILE}" \
    --afreq "${metrics_afreq}" \
    --vmiss "${metrics_vmiss}" \
    --liftover "${SBAYESRC_LIFTOVER_FILE}" \
    --output-dir "${LOCAL_HQ_DIRECT_DIR}" \
    --af-diff-max "${HQ_AF_DIFF_MAX}" \
    --maf-min "${HQ_EUR_MAF_MIN}" \
    --missing-max "${HQ_EUR_MISSING_MAX}"

final_extract="${LOCAL_HQ_DIRECT_DIR}/chr1_22_merged_hq.extract.txt"
final_summary="${LOCAL_HQ_DIRECT_DIR}/chr1_22_merged_hq.filter_summary.tsv"
final_variant_qc="${LOCAL_HQ_DIRECT_DIR}/chr1_22_merged_hq.variant_qc.tsv"
final_params="${LOCAL_HQ_DIRECT_DIR}/chr1_22_merged_hq.params.tsv"
expected_final=$(wc -l < "${final_extract}")
if [[ "${expected_final}" -le 0 ]]; then
    echo "ERROR: high-quality direct extract list is empty" >&2
    exit 1
fi

echo "  Staging high-quality filter outputs to workspace bucket ..."
gcloud storage cp \
    "${final_extract}" \
    "${final_summary}" \
    "${final_variant_qc}" \
    "${final_params}" \
    "${DX_HQ_DIRECT_BFILE_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

if [[ -s "${hq_prefix}.bed" && -s "${hq_prefix}.bim" && -s "${hq_prefix}.fam" ]]; then
    observed_final=$(wc -l < "${hq_prefix}.bim")
else
    observed_final=0
fi

if [[ "${observed_final}" -eq "${expected_final}" && -s "${sample_summary}" ]]; then
    echo "  High-quality direct bfile already has expected ${expected_final} variants — skipping bfile dsub"
else
    bfile_log="${SCRIPT_DIR}/logs/dsub_hq_direct_bfile_$(date +%Y%m%d_%H%M%S).dsub.out"
    mkdir -p "$(dirname "${bfile_log}")"
    echo "  Submitting dsub job to build high-quality direct bfile ..."
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
        --name "sbayesrc-hq-direct-bfile" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_hq_direct_bfile_worker.sh" \
        --env EXPECTED_VARIANTS="${expected_final}" \
        --input PLINK2="${DSUB_PLINK2_GS}" \
        --input BED="${direct_uri_prefix}.bed" \
        --input BIM="${direct_uri_prefix}.bim" \
        --input FAM="${direct_uri_prefix}.fam" \
        --input EXTRACT="${hq_uri_prefix}.extract.txt" \
        --input KEEP="${DSUB_HQ_DIRECT_PREP_URI}/aou_eur.keep" \
        --output-recursive OUTDIR="${DX_HQ_DIRECT_BFILE_URI}/" \
        --min-cores "${HQ_DIRECT_DSUB_MIN_CORES}" \
        --min-ram "${HQ_DIRECT_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${HQ_DIRECT_DSUB_DISK_SIZE}" \
        --disk-type "${HQ_DIRECT_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${bfile_log}"
    dsub_rc=${PIPESTATUS[0]}
    if [[ "${dsub_rc}" -ne 0 ]]; then
        echo "ERROR: dsub high-quality direct bfile job returned ${dsub_rc}" >&2
        echo "Check logs at ${DSUB_LOG_URI}/" >&2
        exit "${dsub_rc}"
    fi
fi

echo "  Verifying high-quality direct bfile ..."
if [[ ! -s "${hq_prefix}.bed" || ! -s "${hq_prefix}.bim" || ! -s "${hq_prefix}.fam" ]]; then
    echo "ERROR: missing high-quality direct bfile outputs under ${DX_HQ_DIRECT_BFILE_DIR}" >&2
    exit 1
fi
observed_final=$(wc -l < "${hq_prefix}.bim")
if [[ "${observed_final}" -ne "${expected_final}" ]]; then
    echo "ERROR: high-quality direct bfile has ${observed_final}/${expected_final} variants" >&2
    exit 1
fi

echo "  Done: ${hq_prefix}.{bed,bim,fam} (${observed_final} variants)"
echo "  Filter summary: ${summary}"
echo "  Sample missingness summary: ${sample_summary}"
