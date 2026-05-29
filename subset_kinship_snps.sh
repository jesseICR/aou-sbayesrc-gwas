#!/bin/bash
# subset_kinship_snps.sh - Build HQ-direct SNP subset for KING kinship.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${LOCAL_SNP_QC_FILE:?LOCAL_SNP_QC_FILE not set}"
: "${LOCAL_KINSHIP_DIR:?LOCAL_KINSHIP_DIR not set}"
: "${UKB_SNP_QC_URL:?UKB_SNP_QC_URL not set}"
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
KINSHIP_SUBSET_DSUB_MIN_CORES="${KINSHIP_SUBSET_DSUB_MIN_CORES:-8}"
KINSHIP_SUBSET_DSUB_MIN_RAM="${KINSHIP_SUBSET_DSUB_MIN_RAM:-32}"
KINSHIP_SUBSET_DSUB_DISK_SIZE="${KINSHIP_SUBSET_DSUB_DISK_SIZE:-300}"
KINSHIP_SUBSET_DSUB_DISK_TYPE="${KINSHIP_SUBSET_DSUB_DISK_TYPE:-pd-ssd}"

mkdir -p "${LOCAL_KINSHIP_DIR}" "${DX_KINSHIP_DIR}"

hq_prefix="${DX_HQ_DIRECT_BFILE_DIR}/chr1_22_merged_hq"
hq_uri_prefix="${DX_HQ_DIRECT_BFILE_URI}/chr1_22_merged_hq"
params="${DX_KINSHIP_DIR}/kinship_snp_subset.params.tsv"
summary="${DX_KINSHIP_DIR}/kinship_snp_subset_summary.tsv"
final_extract="${DX_KINSHIP_DIR}/ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt"
vmiss="${DX_KINSHIP_DIR}/kinship_snp_subset_all_sample_missingness.vmiss"
threshold_counts="${DX_KINSHIP_DIR}/kinship_snp_missingness_threshold_counts.tsv"

if [[ -s "${summary}" && -s "${params}" && -s "${final_extract}" && -s "${vmiss}" ]]; then
    existing_source=$(awk -F'\t' '$1 == "source_bfile" {print $2; exit}' "${params}")
    existing_missing=$(awk -F'\t' '$1 == "kinship_missing_max_exclusive" {print $2; exit}' "${params}")
    expected=$(awk -F'\t' '$1 == "n_intersection_and_missing_lt_threshold" {print $2; exit}' "${summary}")
    observed=$(wc -l < "${final_extract}")
    if [[ "${existing_source}" == "direct_bfile_hq/chr1_22_merged_hq" \
          && "${existing_missing}" == "${KINSHIP_MISSING_MAX}" \
          && -n "${expected}" \
          && "${observed}" -eq "${expected}" ]]; then
        echo "  Kinship SNP subset already exists (${observed} SNPs passing missingness < ${KINSHIP_MISSING_MAX}) — skipping"
        exit 0
    fi
    echo "  Kinship SNP subset exists but params/counts do not match — rebuilding"
fi

for ext in bed bim fam; do
    if [[ ! -s "${hq_prefix}.${ext}" ]]; then
        echo "ERROR: missing HQ direct bfile input ${hq_prefix}.${ext}" >&2
        exit 1
    fi
done

if [[ -s "${LOCAL_SNP_QC_FILE}" ]]; then
    echo "  ukb_snp_qc.txt already cached locally — skipping download"
else
    echo "  Downloading ukb_snp_qc.txt ..."
    tmp_snp_qc="${LOCAL_SNP_QC_FILE}.tmp"
    rm -f "${tmp_snp_qc}"
    curl -fsSL -o "${tmp_snp_qc}" "${UKB_SNP_QC_URL}"
    mv "${tmp_snp_qc}" "${LOCAL_SNP_QC_FILE}"
    echo "  Downloaded ($(wc -l < "${LOCAL_SNP_QC_FILE}") lines)"
fi

header=$(head -1 "${LOCAL_SNP_QC_FILE}")
rel_col=$(echo "${header}" | tr ' ' '\n' | grep -n '^in_Relatedness$' | cut -d: -f1)
if [[ -z "${rel_col}" ]]; then
    echo "ERROR: in_Relatedness column not found in ${LOCAL_SNP_QC_FILE}" >&2
    exit 1
fi
echo "  in_Relatedness is column ${rel_col}"

related_rsids="${LOCAL_KINSHIP_DIR}/ukb_in_relatedness_rsids.txt"
intersection="${LOCAL_KINSHIP_DIR}/ukbb_relatedness_snps_in_hq_direct.txt"
awk -v col="${rel_col}" 'NR > 1 && $col == 1 {print $1}' "${LOCAL_SNP_QC_FILE}" \
    | sort -u > "${related_rsids}"
ukb_relatedness_snps=$(wc -l < "${related_rsids}")
hq_variants=$(wc -l < "${hq_prefix}.bim")
hq_samples=$(wc -l < "${hq_prefix}.fam")
awk 'NR == FNR {keep[$1] = 1; next} ($2 in keep) {print $2}' \
    "${related_rsids}" "${hq_prefix}.bim" > "${intersection}"
n_intersection=$(wc -l < "${intersection}")

echo "  UKB in_Relatedness SNPs: ${ukb_relatedness_snps}"
echo "  HQ direct bfile variants: ${hq_variants}"
echo "  n_intersection_hq_direct: ${n_intersection}"
if [[ "${n_intersection}" -le 0 ]]; then
    echo "ERROR: kinship SNP intersection is empty" >&2
    exit 1
fi

desired_params="${LOCAL_KINSHIP_DIR}/kinship_snp_subset.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
    printf 'source_bfile_variants\t%s\n' "${hq_variants}"
    printf 'source_bfile_samples\t%s\n' "${hq_samples}"
    printf 'ukb_snp_qc_url\t%s\n' "${UKB_SNP_QC_URL}"
    printf 'ukb_snp_qc_file_size\t%s\n' "$(stat -c%s "${LOCAL_SNP_QC_FILE}")"
    printf 'ukb_relatedness_snps\t%s\n' "${ukb_relatedness_snps}"
    printf 'n_intersection_hq_direct\t%s\n' "${n_intersection}"
    printf 'kinship_missing_max_exclusive\t%s\n' "${KINSHIP_MISSING_MAX}"
} > "${desired_params}"

write_subset_from_existing_vmiss() {
    local tmp_extract="${final_extract}.tmp"
    local tmp_counts="${threshold_counts}.tmp"
    local stats

    rm -f "${tmp_extract}" "${tmp_counts}"
    if ! stats=$(awk \
        -v max_missing="${KINSHIP_MISSING_MAX}" \
        -v intersection="${intersection}" \
        -v out_extract="${tmp_extract}" \
        -v out_counts="${tmp_counts}" '
        BEGIN {
            while ((getline id < intersection) > 0) {
                if (id != "") {
                    keep[id] = 1
                    n_intersection++
                }
            }
            close(intersection)
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "ID") id_col = i
                if ($i == "F_MISS") fmiss_col = i
            }
            next
        }
        id_col > 0 && fmiss_col > 0 {
            id = $id_col
            fmiss = $fmiss_col + 0
            measured++
            if (!(id in keep)) {
                unexpected_ids++
                next
            }
            if (fmiss < max_missing) {
                print id > out_extract
                pass_threshold++
            }
            if (fmiss < 0.05) pass_005++
            if (fmiss < 0.04) pass_004++
            if (fmiss < 0.03) pass_003++
            if (fmiss < 0.02) pass_002++
            if (fmiss < 0.01) pass_001++
        }
        END {
            if (id_col == 0 || fmiss_col == 0) {
                exit 2
            }
            if (unexpected_ids > 0 || measured != n_intersection) {
                exit 3
            }
            print "threshold\tpassing_variants\tfailing_variants" > out_counts
            printf "missingness_lt_0.05\t%d\t%d\n", pass_005, measured - pass_005 > out_counts
            printf "missingness_lt_0.04\t%d\t%d\n", pass_004, measured - pass_004 > out_counts
            printf "missingness_lt_0.03\t%d\t%d\n", pass_003, measured - pass_003 > out_counts
            printf "missingness_lt_0.02\t%d\t%d\n", pass_002, measured - pass_002 > out_counts
            printf "missingness_lt_0.01\t%d\t%d\n", pass_001, measured - pass_001 > out_counts
            printf "total_measured\t%d\t0\n", measured > out_counts
            printf "%d\t%d\n", measured, pass_threshold
        }
    ' "${vmiss}"); then
        rm -f "${tmp_extract}" "${tmp_counts}"
        return 1
    fi

    local n_vmiss n_final
    IFS=$'\t' read -r n_vmiss n_final <<< "${stats}"
    if [[ -z "${n_vmiss}" || -z "${n_final}" || "${n_final}" -le 0 ]]; then
        rm -f "${tmp_extract}" "${tmp_counts}"
        return 1
    fi

    mv "${tmp_extract}" "${final_extract}"
    mv "${tmp_counts}" "${threshold_counts}"
    cp "${desired_params}" "${params}"
    {
        printf 'metric\tvalue\n'
        printf 'source_bfile\t%s\n' "direct_bfile_hq/chr1_22_merged_hq"
        printf 'source_bfile_variants\t%s\n' "${hq_variants}"
        printf 'source_bfile_samples\t%s\n' "${hq_samples}"
        printf 'ukb_relatedness_snps\t%s\n' "${ukb_relatedness_snps}"
        printf 'n_intersection_hq_direct\t%s\n' "${n_intersection}"
        printf 'kinship_missing_max_exclusive\t%s\n' "${KINSHIP_MISSING_MAX}"
        printf 'n_intersection_missingness_measured\t%s\n' "${n_vmiss}"
        printf 'n_intersection_and_missing_lt_%s\t%s\n' "${KINSHIP_MISSING_MAX}" "${n_final}"
        printf 'n_intersection_and_missing_lt_threshold\t%s\n' "${n_final}"
    } > "${summary}"
    echo "  Reused existing all-sample missingness report for threshold < ${KINSHIP_MISSING_MAX} (${n_final} SNPs)"
    return 0
}

if [[ -s "${summary}" && -s "${params}" && -s "${final_extract}" && -s "${vmiss}" ]]; then
    if diff -q "${desired_params}" "${params}" >/dev/null 2>&1; then
        expected=$(awk -F'\t' '$1 == "n_intersection_and_missing_lt_threshold" {print $2; exit}' "${summary}")
        observed=$(wc -l < "${final_extract}")
        if [[ -n "${expected}" && "${observed}" -eq "${expected}" ]]; then
            echo "  Kinship SNP subset already exists (${observed} SNPs passing missingness < ${KINSHIP_MISSING_MAX}) — skipping"
            exit 0
        fi
    fi
    echo "  Kinship SNP subset exists but params/counts do not match — rebuilding"
fi

if [[ -s "${vmiss}" ]]; then
    echo "  Existing all-sample missingness report found — trying local threshold re-filter"
    if write_subset_from_existing_vmiss; then
        echo "  Kinship SNP missingness threshold counts:"
        awk -F'\t' 'NR > 1 {printf "    %s = %s pass, %s fail\n", $1, $2, $3}' "${threshold_counts}"
        echo "  Kinship SNP subset summary:"
        awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
        exit 0
    fi
    echo "  Existing missingness report did not validate against the current SNP intersection — recomputing with dsub"
fi

echo "  Staging plink2 + kinship SNP inputs to workspace bucket ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
export DSUB_KINSHIP_PREP_URI="${DSUB_KINSHIP_PREP_URI:-${DX_KINSHIP_URI}/_prep}"
gcloud storage cp "${intersection}" "${desired_params}" "${DSUB_KINSHIP_PREP_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

subset_log="${SCRIPT_DIR}/logs/dsub_kinship_subset_$(date +%Y%m%d_%H%M%S).dsub.out"
mkdir -p "$(dirname "${subset_log}")"
echo "  Submitting dsub job for all-sample kinship SNP missingness ..."
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
    --name "sbayesrc-kinship-snp-subset" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_kinship_subset_worker.sh" \
    --env KINSHIP_MISSING_MAX="${KINSHIP_MISSING_MAX}" \
    --input PLINK2="${DSUB_PLINK2_GS}" \
    --input BED="${hq_uri_prefix}.bed" \
    --input BIM="${hq_uri_prefix}.bim" \
    --input FAM="${hq_uri_prefix}.fam" \
    --input EXTRACT="${DSUB_KINSHIP_PREP_URI}/ukbb_relatedness_snps_in_hq_direct.txt" \
    --input PARAMS="${DSUB_KINSHIP_PREP_URI}/kinship_snp_subset.desired_params.tsv" \
    --output-recursive OUTDIR="${DX_KINSHIP_URI}/" \
    --min-cores "${KINSHIP_SUBSET_DSUB_MIN_CORES}" \
    --min-ram "${KINSHIP_SUBSET_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${KINSHIP_SUBSET_DSUB_DISK_SIZE}" \
    --disk-type "${KINSHIP_SUBSET_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${subset_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub kinship SNP subset job returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

if [[ ! -s "${summary}" || ! -s "${final_extract}" ]]; then
    echo "ERROR: missing kinship SNP subset outputs under ${DX_KINSHIP_DIR}" >&2
    exit 1
fi

echo "  Kinship SNP subset summary:"
awk -F'\t' 'NR > 1 {printf "    %s = %s\n", $1, $2}' "${summary}"
