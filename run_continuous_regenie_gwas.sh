#!/bin/bash
# run_continuous_regenie_gwas.sh - Run a continuous-trait REGENIE GWAS on AoU.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage: bash run_continuous_regenie_gwas.sh <input_name> <output_name> [OPTIONS]

Runs REGENIE Step 1 on gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1
and Step 2 on gwas_genotypes/step2_wgs_pfiles/chr{1..22}.pgen/.pvar/.psam.

Options:
  --apply-rint            Apply rank-inverse normal transform (default)
  --no-apply-rint         Disable rank-inverse normal transform
  --chroms LIST           Chromosomes to run, e.g. 22, 1,2,3, or 1-22
  --step1-block-size N    REGENIE Step 1 block size
  --step2-block-size N    REGENIE Step 2 block size
  -h, --help              Show this help
EOF
}

if [[ $# -lt 1 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if [[ $# -lt 2 ]]; then
    usage >&2
    exit 1
fi

INPUT_NAME="$1"
OUTPUT_NAME="$2"
shift 2

APPLY_RINT="${REGENIE_APPLY_RINT:-1}"
STEP1_BLOCK_SIZE="${REGENIE_STEP1_BLOCK_SIZE:-1000}"
STEP2_BLOCK_SIZE="${REGENIE_STEP2_BLOCK_SIZE:-200}"
CHROMS_SPEC="${REGENIE_CHROMS:-1-22}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply-rint)
            APPLY_RINT=1
            shift
            ;;
        --no-apply-rint)
            APPLY_RINT=0
            shift
            ;;
        --chroms)
            CHROMS_SPEC="${2:?--chroms requires a value}"
            shift 2
            ;;
        --step1-block-size)
            STEP1_BLOCK_SIZE="${2:?--step1-block-size requires a value}"
            shift 2
            ;;
        --step2-block-size)
            STEP2_BLOCK_SIZE="${2:?--step2-block-size requires a value}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${REGENIE:?REGENIE not set}"
: "${DX_GWAS_STEP1_BFILE_DIR:?DX_GWAS_STEP1_BFILE_DIR not set}"
: "${DX_GWAS_STEP1_BFILE_URI:?DX_GWAS_STEP1_BFILE_URI not set}"
: "${DX_GWAS_STEP2_PFILE_DIR:?DX_GWAS_STEP2_PFILE_DIR not set}"
: "${DX_GWAS_STEP2_PFILE_URI:?DX_GWAS_STEP2_PFILE_URI not set}"
: "${DX_REGENIE_INPUT_DIR:?DX_REGENIE_INPUT_DIR not set}"
: "${DX_REGENIE_INPUT_URI:?DX_REGENIE_INPUT_URI not set}"
: "${DX_REGENIE_OUTPUT_DIR:?DX_REGENIE_OUTPUT_DIR not set}"
: "${DX_REGENIE_OUTPUT_URI:?DX_REGENIE_OUTPUT_URI not set}"
: "${LOCAL_REGENIE_DIR:?LOCAL_REGENIE_DIR not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"
: "${DSUB_REGENIE_BUNDLE_URI:?DSUB_REGENIE_BUNDLE_URI not set}"

if [[ ! -x "${REGENIE}" ]]; then
    echo "ERROR: REGENIE is not executable: ${REGENIE}" >&2
    exit 1
fi

is_uint() {
    [[ "${1:-}" =~ ^[0-9]+$ ]]
}

sanitize_prefix() {
    local value="$1" sanitized
    sanitized="$(printf '%s' "${value}" | tr -c 'A-Za-z0-9_.-' '_' | sed 's/^_*//; s/_*$//')"
    if [[ -z "${sanitized}" ]]; then
        echo "ERROR: output name ${value} cannot be converted to a valid file prefix" >&2
        return 1
    fi
    printf '%s\n' "${sanitized}"
}

regenie_version_text() {
    local help_text
    help_text="$("$1" --help 2>&1 || true)"
    printf '%s\n' "${help_text}" | sed -n '2{s/^ *//;p;q;}'
}

expand_chroms() {
    local spec="$1"
    local -a out=()
    local token start end c
    IFS=',' read -ra tokens <<< "${spec}"
    for token in "${tokens[@]}"; do
        if [[ "${token}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            start="${BASH_REMATCH[1]}"
            end="${BASH_REMATCH[2]}"
            if (( start > end )); then
                echo "ERROR: invalid chromosome range ${token}" >&2
                return 1
            fi
            for ((c = start; c <= end; c++)); do
                out+=("${c}")
            done
        elif [[ "${token}" =~ ^[0-9]+$ ]]; then
            out+=("${token}")
        else
            echo "ERROR: invalid chromosome token ${token}" >&2
            return 1
        fi
    done
    for c in "${out[@]}"; do
        if (( c < 1 || c > 22 )); then
            echo "ERROR: chromosome ${c} outside 1..22" >&2
            return 1
        fi
    done
    printf '%s\n' "${out[@]}" | awk '!seen[$1]++'
}

mapfile -t CHROMS < <(expand_chroms "${CHROMS_SPEC}")
if [[ "${#CHROMS[@]}" -eq 0 ]]; then
    echo "ERROR: no chromosomes selected by ${CHROMS_SPEC}" >&2
    exit 1
fi

stage_regenie_bundle() {
    local regenie_sha bundle local_bundle tmp bundle_uri version regenie_lib_dir
    mkdir -p "${SCRIPT_DIR}/tools" "${LOCAL_REGENIE_DIR}"
    regenie_sha=$(sha256sum "${REGENIE}" | awk '{print $1}')
    local_bundle="${SCRIPT_DIR}/tools/regenie_${regenie_sha:0:16}_mkl.tar.gz"
    bundle_uri="${DSUB_REGENIE_BUNDLE_URI}/$(basename "${local_bundle}")"
    if [[ ! -s "${local_bundle}" ]]; then
        echo "  Creating minimal REGENIE runtime bundle ..." >&2
        tmp=$(mktemp -d)
        cp "${REGENIE}" "${tmp}/regenie"
        ldd "${REGENIE}" | awk '/=>/ {print $3}' | while read -r lib; do
            if [[ -f "${lib}" ]]; then
                cp -L "${lib}" "${tmp}/"
            fi
        done
        regenie_lib_dir="$(cd "$(dirname "${REGENIE}")/../lib" && pwd)"
        find "${regenie_lib_dir}" -maxdepth 1 -type f -name 'libmkl*.so*' -print0 |
            xargs -0 -r -I{} cp -L "{}" "${tmp}/"
        version=$(regenie_version_text "${REGENIE}")
        {
            printf 'regenie\t%s\n' "${version}"
            printf 'regenie_sha256\t%s\n' "${regenie_sha}"
            printf 'bundle_contents\t%s\n' "ldd libraries plus Workbench libmkl*.so* runtime libraries"
            printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        } > "${tmp}/manifest.tsv"
        tar -C "${tmp}" -czf "${local_bundle}" .
        rm -rf "${tmp}"
    fi
    if ! gcloud storage ls "${bundle_uri}" --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1; then
        echo "  Uploading REGENIE bundle to ${bundle_uri} ..." >&2
        gcloud storage cp "${local_bundle}" "${bundle_uri}" \
            --billing-project="${GOOGLE_PROJECT}" >&2
    else
        echo "  REGENIE bundle already staged: ${bundle_uri}" >&2
    fi
    printf '%s\n' "${bundle_uri}"
}

input_dir="${DX_REGENIE_INPUT_DIR}/${INPUT_NAME}"
input_uri="${DX_REGENIE_INPUT_URI}/${INPUT_NAME}"
output_dir="${DX_REGENIE_OUTPUT_DIR}/${OUTPUT_NAME}"
output_uri="${DX_REGENIE_OUTPUT_URI}/${OUTPUT_NAME}"
step1_dir="${output_dir}/step1"
step1_uri="${output_uri}/step1"
step2_dir="${output_dir}/step2"
step2_uri="${output_uri}/step2"
result_prefix="$(sanitize_prefix "${OUTPUT_NAME}")"
step1_output_prefix="${result_prefix}_step1"
mkdir -p "${output_dir}" "${step1_dir}" "${step2_dir}" "${LOCAL_REGENIE_DIR}"

phen="${input_dir}/phen.txt"
covar="${input_dir}/covar.txt"
keep="${input_dir}/training_iids.txt"

find_input_metadata_file() {
    local suffix="$1" candidate matches
    for candidate in \
        "${input_dir}/${INPUT_NAME}.${suffix}" \
        "${input_dir}/${OUTPUT_NAME}.${suffix}" \
        "${input_dir}/height_gwas.${suffix}" \
        "${input_dir}/ea_gwas.${suffix}" \
        "${input_dir}/income_gwas.${suffix}"; do
        if [[ -s "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done
    mapfile -t matches < <(find "${input_dir}" -maxdepth 1 -type f -name "*.${suffix}" -print | sort)
    if [[ "${#matches[@]}" -eq 1 && -s "${matches[0]}" ]]; then
        printf '%s\n' "${matches[0]}"
        return 0
    fi
    echo "ERROR: could not identify unique input metadata file *.${suffix} in ${input_dir}" >&2
    if [[ "${#matches[@]}" -gt 1 ]]; then
        printf '  candidate: %s\n' "${matches[@]}" >&2
    fi
    return 1
}

input_params="$(find_input_metadata_file params.tsv)"
input_summary="$(find_input_metadata_file summary.tsv)"
for f in "${phen}" "${covar}" "${keep}" "${input_params}" "${input_summary}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing REGENIE input file ${f}" >&2
        exit 1
    fi
done

step1_prefix="${DX_GWAS_STEP1_BFILE_DIR}/chr1_22_merged_gwas_step1"
step1_uri_prefix="${DX_GWAS_STEP1_BFILE_URI}/chr1_22_merged_gwas_step1"
for ext in bed bim fam; do
    if [[ ! -s "${step1_prefix}.${ext}" ]]; then
        echo "ERROR: missing final GWAS Step 1 bfile input ${step1_prefix}.${ext}" >&2
        exit 1
    fi
done

pheno_col="${REGENIE_PHENO_COL:-height}"
if [[ -n "${REGENIE_COVAR_COLS:-}" ]]; then
    covar_cols="${REGENIE_COVAR_COLS}"
else
    n_pcs=$(awk -F'\t' '$1 == "n_pcs" {print $2; exit}' "${input_summary}")
    n_pcs="${n_pcs:-10}"
    covar_cols="age_c,sex_c,age_c_sex_c_inter"
    for i in $(seq 1 "${n_pcs}"); do
        covar_cols+=",PC${i}_AVG"
    done
fi

keep_samples=$(wc -l < "${keep}")
step1_variants=$(wc -l < "${step1_prefix}.bim")
step1_samples=$(wc -l < "${step1_prefix}.fam")
regenie_version=$(regenie_version_text "${REGENIE}")
regenie_sha=$(sha256sum "${REGENIE}" | awk '{print $1}')

chroms_joined=$(IFS=,; echo "${CHROMS[*]}")
desired_params="${LOCAL_REGENIE_DIR}/${OUTPUT_NAME}.regenie_gwas.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'input_name\t%s\n' "${INPUT_NAME}"
    printf 'output_name\t%s\n' "${OUTPUT_NAME}"
    printf 'result_prefix\t%s\n' "${result_prefix}"
    printf 'step1_bfile\t%s\n' "gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1"
    printf 'step1_bfile_variants\t%s\n' "${step1_variants}"
    printf 'step1_bfile_samples\t%s\n' "${step1_samples}"
    printf 'keep_samples\t%s\n' "${keep_samples}"
    printf 'pheno_col\t%s\n' "${pheno_col}"
    printf 'covar_cols\t%s\n' "${covar_cols}"
    printf 'apply_rint\t%s\n' "${APPLY_RINT}"
    printf 'step1_block_size\t%s\n' "${STEP1_BLOCK_SIZE}"
    printf 'step2_block_size\t%s\n' "${STEP2_BLOCK_SIZE}"
    printf 'chroms\t%s\n' "${chroms_joined}"
    printf 'regenie_version\t%s\n' "${regenie_version}"
    printf 'regenie_sha256\t%s\n' "${regenie_sha}"
    printf 'phen_size\t%s\n' "$(stat -c%s "${phen}")"
    printf 'covar_size\t%s\n' "$(stat -c%s "${covar}")"
    printf 'keep_size\t%s\n' "$(stat -c%s "${keep}")"
    printf 'input_params_size\t%s\n' "$(stat -c%s "${input_params}")"
    for c in "${CHROMS[@]}"; do
        summary="${DX_GWAS_STEP2_PFILE_DIR}/chr${c}.summary.tsv"
        if [[ ! -s "${summary}" ]]; then
            echo "ERROR: missing final GWAS Step 2 pfile summary ${summary}" >&2
            exit 1
        fi
        variants=$(awk -F'\t' '$1 == "final_variants" {print $2; exit}' "${summary}")
        if ! is_uint "${variants}" || [[ "${variants}" -le 0 ]]; then
            echo "ERROR: could not read extracted variant count from ${summary}" >&2
            exit 1
        fi
        printf 'chr%s_variants\t%s\n' "${c}" "${variants}"
        printf 'chr%s_pgen_size\t%s\n' "${c}" "$(stat -c%s "${DX_GWAS_STEP2_PFILE_DIR}/chr${c}.pgen")"
        printf 'chr%s_pvar_size\t%s\n' "${c}" "$(stat -c%s "${DX_GWAS_STEP2_PFILE_DIR}/chr${c}.pvar")"
        printf 'chr%s_psam_size\t%s\n' "${c}" "$(stat -c%s "${DX_GWAS_STEP2_PFILE_DIR}/chr${c}.psam")"
    done
} > "${desired_params}"

echo "=== REGENIE continuous GWAS ==="
echo "  input        = ${input_dir}"
echo "  output       = ${output_dir}"
echo "  file prefix  = ${result_prefix}"
echo "  chroms       = ${chroms_joined}"
echo "  keep samples = ${keep_samples}"
echo "  Step 1 bfile = ${step1_variants} variants, ${step1_samples} samples"
echo "  RINT         = ${APPLY_RINT}"
echo "  covars       = ${covar_cols}"

bundle_uri=$(stage_regenie_bundle)
gcloud storage cp "${desired_params}" "${output_uri}/regenie_gwas.desired_params.tsv" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

step1_summary="${step1_dir}/regenie_step1.summary.tsv"
step1_params="${step1_dir}/regenie_step1.params.tsv"
step1_pred="${step1_dir}/${step1_output_prefix}_pred.list"
run_step1=1
if [[ -s "${step1_summary}" && -s "${step1_params}" && -s "${step1_pred}" ]]; then
    if diff -q "${desired_params}" "${step1_params}" >/dev/null 2>&1; then
        observed_keep=$(awk -F'\t' '$1 == "keep_samples" {print $2; exit}' "${step1_summary}")
        if [[ "${observed_keep}" == "${keep_samples}" ]]; then
            run_step1=0
        fi
    fi
fi

if [[ "${run_step1}" -eq 1 ]]; then
    echo "  Submitting REGENIE Step 1 ..."
    step1_log="${SCRIPT_DIR}/logs/dsub_regenie_step1_${OUTPUT_NAME}_$(date +%Y%m%d_%H%M%S).dsub.out"
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
        --name "sbayesrc-regenie-step1" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_regenie_step1_worker.sh" \
        --env PHENO_COL="${pheno_col}" \
        --env COVAR_COLS="${covar_cols}" \
        --env APPLY_RINT="${APPLY_RINT}" \
        --env STEP1_BLOCK_SIZE="${STEP1_BLOCK_SIZE}" \
        --env RESULT_PREFIX="${result_prefix}" \
        --env EXPECTED_KEEP_SAMPLES="${keep_samples}" \
        --env EXPECTED_BFILE_VARIANTS="${step1_variants}" \
        --env EXPECTED_BFILE_SAMPLES="${step1_samples}" \
        --input REGENIE_BUNDLE="${bundle_uri}" \
        --input BED="${step1_uri_prefix}.bed" \
        --input BIM="${step1_uri_prefix}.bim" \
        --input FAM="${step1_uri_prefix}.fam" \
        --input PHEN="${input_uri}/phen.txt" \
        --input COVAR="${input_uri}/covar.txt" \
        --input KEEP="${input_uri}/training_iids.txt" \
        --input PARAMS="${output_uri}/regenie_gwas.desired_params.tsv" \
        --output-recursive OUTDIR="${step1_uri}/" \
        --min-cores "${REGENIE_STEP1_DSUB_MIN_CORES:-16}" \
        --min-ram "${REGENIE_STEP1_DSUB_MIN_RAM:-64}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${REGENIE_STEP1_DSUB_DISK_SIZE:-300}" \
        --disk-type "${REGENIE_STEP1_DSUB_DISK_TYPE:-pd-ssd}" \
        --wait \
        --summary 2>&1 | tee "${step1_log}"
    dsub_rc=${PIPESTATUS[0]}
    if [[ "${dsub_rc}" -ne 0 ]]; then
        echo "ERROR: REGENIE Step 1 dsub job returned ${dsub_rc}" >&2
        exit "${dsub_rc}"
    fi
else
    echo "  REGENIE Step 1 already exists (${keep_samples} samples) — skipping"
fi

if [[ ! -s "${step1_pred}" ]]; then
    echo "ERROR: missing Step 1 prediction list after Step 1: ${step1_pred}" >&2
    exit 1
fi

tasks_tsv="${SCRIPT_DIR}/logs/dsub_regenie_step2_${OUTPUT_NAME}_$(date +%Y%m%d_%H%M%S).tsv"
{
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '--env CHROM' '--env EXPECTED_VARIANTS' \
        '--input PGEN' '--input PVAR' '--input PSAM' \
        '--output-recursive OUTDIR' '--env PHENO_COL' '--env COVAR_COLS'
    for c in "${CHROMS[@]}"; do
        chrom="chr${c}"
        chrom_result_prefix="${chrom}_${result_prefix}"
        chrom_summary="${step2_dir}/${chrom}/${chrom_result_prefix}.summary.tsv"
        submit=1
        if [[ -s "${chrom_summary}" && -s "${step2_dir}/${chrom}/${chrom_result_prefix}.params.tsv" ]]; then
            if diff -q "${desired_params}" "${step2_dir}/${chrom}/${chrom_result_prefix}.params.tsv" >/dev/null 2>&1; then
                submit=0
            fi
        fi
        if [[ "${submit}" -eq 1 ]]; then
            variants=$(awk -F'\t' -v key="chr${c}_variants" '$1 == key {print $2; exit}' "${desired_params}")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${chrom}" "${variants}" \
                "${DX_GWAS_STEP2_PFILE_URI}/${chrom}.pgen" \
                "${DX_GWAS_STEP2_PFILE_URI}/${chrom}.pvar" \
                "${DX_GWAS_STEP2_PFILE_URI}/${chrom}.psam" \
                "${step2_uri}/${chrom}/" \
                "${pheno_col}" "${covar_cols}"
        fi
    done
} > "${tasks_tsv}"

to_submit=$(( $(wc -l < "${tasks_tsv}") - 1 ))
echo "  REGENIE Step 2 chromosomes selected: ${#CHROMS[@]}; to submit: ${to_submit}"

if [[ "${to_submit}" -gt 0 ]]; then
    step2_log="${tasks_tsv%.tsv}.dsub.out"
    dsub_rc=0
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
        --name "sbayesrc-regenie-step2" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_regenie_step2_worker.sh" \
        --tasks "${tasks_tsv}" \
        --env APPLY_RINT="${APPLY_RINT}" \
        --env STEP2_BLOCK_SIZE="${STEP2_BLOCK_SIZE}" \
        --env RESULT_PREFIX="${result_prefix}" \
        --env EXPECTED_KEEP_SAMPLES="${keep_samples}" \
        --input REGENIE_BUNDLE="${bundle_uri}" \
        --input PHEN="${input_uri}/phen.txt" \
        --input COVAR="${input_uri}/covar.txt" \
        --input KEEP="${input_uri}/training_iids.txt" \
        --input PARAMS="${output_uri}/regenie_gwas.desired_params.tsv" \
        --input-recursive STEP1_DIR="${step1_uri}/" \
        --min-cores "${REGENIE_STEP2_DSUB_MIN_CORES:-8}" \
        --min-ram "${REGENIE_STEP2_DSUB_MIN_RAM:-32}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${REGENIE_STEP2_DSUB_DISK_SIZE:-300}" \
        --disk-type "${REGENIE_STEP2_DSUB_DISK_TYPE:-pd-ssd}" \
        --wait \
        --summary 2>&1 | tee "${step2_log}"
    dsub_rc=${PIPESTATUS[0]}

    dsub_job_id="$(awk '/^Launched job-id:/ {print $NF; exit}' "${step2_log}")"
    if [[ -n "${dsub_job_id}" ]]; then
        echo "  Polling dstat for Step 2 job ${dsub_job_id} until all ${to_submit} tasks terminal ..."
        while true; do
            terminal_count=$(dstat --provider "${DSUB_PROVIDER}" \
                                   --project "${GOOGLE_PROJECT}" \
                                   --location "${DSUB_REGION}" \
                                   --jobs "${dsub_job_id}" \
                                   --users jupyter \
                                   --status '*' 2>/dev/null \
                             | awk 'NR>2 && /SUCCESS|FAILURE|CANCEL/ {c++} END {print c+0}')
            if (( terminal_count >= to_submit )); then
                echo "  ${terminal_count}/${to_submit} Step 2 tasks terminal — proceeding"
                break
            fi
            echo "  $(date -u +%H:%M:%SZ) ${terminal_count}/${to_submit} terminal — waiting 30s ..."
            sleep 30
        done
    fi

    if [[ "${dsub_rc}" -ne 0 ]]; then
        echo "ERROR: REGENIE Step 2 dsub job returned ${dsub_rc}" >&2
        exit "${dsub_rc}"
    fi
fi

echo "  Verifying REGENIE Step 2 outputs ..."
total_tested=0
for c in "${CHROMS[@]}"; do
    chrom="chr${c}"
    chrom_result_prefix="${chrom}_${result_prefix}"
    chrom_summary="${step2_dir}/${chrom}/${chrom_result_prefix}.summary.tsv"
    chrom_params="${step2_dir}/${chrom}/${chrom_result_prefix}.params.tsv"
    if [[ ! -s "${chrom_summary}" || ! -s "${chrom_params}" ]]; then
        echo "ERROR: missing Step 2 summary/params for ${chrom}" >&2
        exit 1
    fi
    if ! diff -q "${desired_params}" "${chrom_params}" >/dev/null 2>&1; then
        echo "ERROR: Step 2 params mismatch for ${chrom}" >&2
        exit 1
    fi
    tested=$(awk -F'\t' '$1 == "tested_variants" {print $2; exit}' "${chrom_summary}")
    if ! is_uint "${tested}"; then
        echo "ERROR: could not read tested variant count from ${chrom_summary}" >&2
        exit 1
    fi
    total_tested=$((total_tested + tested))
done

cp "${desired_params}" "${output_dir}/regenie_gwas.params.tsv"
{
    printf 'metric\tvalue\n'
    printf 'input_name\t%s\n' "${INPUT_NAME}"
    printf 'output_name\t%s\n' "${OUTPUT_NAME}"
    printf 'chromosomes\t%s\n' "${chroms_joined}"
    printf 'step1_samples\t%s\n' "${keep_samples}"
    printf 'step1_variants\t%s\n' "${step1_variants}"
    printf 'step2_chromosomes\t%s\n' "${#CHROMS[@]}"
    printf 'step2_total_tested_variants\t%s\n' "${total_tested}"
    printf 'apply_rint\t%s\n' "${APPLY_RINT}"
    printf 'pheno_col\t%s\n' "${pheno_col}"
    printf 'covar_cols\t%s\n' "${covar_cols}"
} > "${output_dir}/regenie_gwas.summary.tsv"

echo "  Done: Step 1 samples=${keep_samples}; Step 2 tested variants=${total_tested}"
