#!/bin/bash
# run_king_kinship.sh - Run KING kinship on the HQ direct bfile.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_HQ_DIRECT_BFILE_DIR:?DX_HQ_DIRECT_BFILE_DIR not set}"
: "${DX_HQ_DIRECT_BFILE_URI:?DX_HQ_DIRECT_BFILE_URI not set}"
: "${DX_KINSHIP_DIR:?DX_KINSHIP_DIR not set}"
: "${DX_KINSHIP_URI:?DX_KINSHIP_URI not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

KINSHIP_MISSING_MAX="${KINSHIP_MISSING_MAX:-0.01}"
KING_TABLE_FILTER="${KING_TABLE_FILTER:-0.035}"
KING_DSUB_MIN_CORES="${KING_DSUB_MIN_CORES:-32}"
KING_DSUB_MIN_RAM="${KING_DSUB_MIN_RAM:-256}"
KING_DSUB_DISK_SIZE="${KING_DSUB_DISK_SIZE:-300}"
KING_DSUB_DISK_TYPE="${KING_DSUB_DISK_TYPE:-pd-ssd}"

mkdir -p "${DX_KINSHIP_DIR}"

hq_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
hq_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"
extract="${DX_KINSHIP_DIR}/ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt"
subset_params="${DX_KINSHIP_DIR}/kinship_snp_subset.params.tsv"
subset_summary="${DX_KINSHIP_DIR}/kinship_snp_subset_summary.tsv"
for ext in bed bim fam; do
    if [[ ! -s "${hq_prefix}.${ext}" ]]; then
        echo "ERROR: missing HQ direct bfile input ${hq_prefix}.${ext}" >&2
        exit 1
    fi
done
if [[ ! -s "${extract}" || ! -s "${subset_params}" || ! -s "${subset_summary}" ]]; then
    echo "ERROR: missing kinship SNP subset outputs; run subset_kinship_snps.sh first" >&2
    exit 1
fi

kinship_snps=$(wc -l < "${extract}")
hq_variants=$(wc -l < "${hq_prefix}.bim")
hq_samples=$(wc -l < "${hq_prefix}.fam")
if [[ "${kinship_snps}" -le 0 ]]; then
    echo "ERROR: final kinship SNP extract is empty: ${extract}" >&2
    exit 1
fi

desired_params="${DX_KINSHIP_DIR}/aou_hq_direct_rel.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'source_bfile_variants\t%s\n' "${hq_variants}"
    printf 'source_bfile_samples\t%s\n' "${hq_samples}"
    printf 'kinship_snps\t%s\n' "${kinship_snps}"
    printf 'kinship_missing_max_exclusive\t%s\n' "${KINSHIP_MISSING_MAX}"
    printf 'king_table_filter\t%s\n' "${KING_TABLE_FILTER}"
    printf 'subset_params_size\t%s\n' "$(stat -c%s "${subset_params}")"
    printf 'subset_summary_size\t%s\n' "$(stat -c%s "${subset_summary}")"
} > "${desired_params}"

kin0="${DX_KINSHIP_DIR}/aou_hq_direct_rel.kin0"
params="${DX_KINSHIP_DIR}/aou_hq_direct_rel.params.tsv"
summary="${DX_KINSHIP_DIR}/aou_hq_direct_rel.summary.tsv"
if [[ -s "${kin0}" && -s "${params}" && -s "${summary}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected_pairs=$(awk -F'\t' '$1 == "king_pairs" {print $2; exit}' "${summary}")
        observed_pairs=$(( $(wc -l < "${kin0}") - 1 ))
        if [[ -n "${expected_pairs}" && "${observed_pairs}" -eq "${expected_pairs}" ]]; then
            echo "  KING kinship already exists (${observed_pairs} pairs; filter >= ${KING_TABLE_FILTER}) — skipping"
            exit 0
        fi
    fi
    echo "  KING kinship exists but params/counts do not match — rebuilding"
fi

echo "  Staging plink2 for KING ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

king_log="${SCRIPT_DIR}/logs/dsub_king_kinship_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${king_log}")"
echo "  Submitting dsub job for KING kinship on ${kinship_snps} SNPs ..."
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
    --name "sbayesrc-king-kinship" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_king_kinship_worker.sh" \
    --env KING_TABLE_FILTER="${KING_TABLE_FILTER}" \
    --env EXPECTED_SNPS="${kinship_snps}" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input BED="${hq_uri_prefix}.bed" \
    --input BIM="${hq_uri_prefix}.bim" \
    --input FAM="${hq_uri_prefix}.fam" \
    --input EXTRACT="${DX_KINSHIP_URI}/ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt" \
    --input PARAMS="${DX_KINSHIP_URI}/aou_hq_direct_rel.desired_params.tsv" \
    --output-recursive OUTDIR="${DX_KINSHIP_URI}/" \
    --min-cores "${KING_DSUB_MIN_CORES}" \
    --min-ram "${KING_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${KING_DSUB_DISK_SIZE}" \
    --disk-type "${KING_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${king_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub KING kinship job returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

if [[ ! -s "${kin0}" || ! -s "${summary}" ]]; then
    echo "ERROR: missing KING outputs under ${DX_KINSHIP_DIR}" >&2
    exit 1
fi

echo "  KING kinship summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
