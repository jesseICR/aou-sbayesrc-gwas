#!/bin/bash
# extract_direct_snps.sh — Build per-chromosome direct-SNP pfiles via dsub.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${PLINK2:?PLINK2 not set}"
: "${DX_WGS_PFILE_DIR:?DX_WGS_PFILE_DIR not set}"
: "${DX_WGS_PFILE_URI:?DX_WGS_PFILE_URI not set}"
: "${DX_DIRECT_PFILE_DIR:?DX_DIRECT_PFILE_DIR not set}"
: "${DX_DIRECT_PFILE_URI:?DX_DIRECT_PFILE_URI not set}"
: "${LOCAL_DIRECT_PREP_DIR:?LOCAL_DIRECT_PREP_DIR not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_PLINK2_GS:?DSUB_PLINK2_GS not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

mkdir -p "${DX_DIRECT_PFILE_DIR}"

summary="${LOCAL_DIRECT_PREP_DIR}/summary.tsv"
if [[ ! -s "${summary}" ]]; then
    echo "ERROR: missing direct SNP summary: ${summary}" >&2
    echo "Run prepare_direct_snps.py before extract_direct_snps.sh." >&2
    exit 1
fi

echo "  Staging direct SNP metadata + plink2 to workspace bucket ..."
gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

export DSUB_DIRECT_PREP_URI="${DSUB_DIRECT_PREP_URI:-${DX_DIRECT_PFILE_URI}/_prep}"
gcloud storage cp \
    "${LOCAL_DIRECT_PREP_DIR}"/chr*.extract.txt \
    "${LOCAL_DIRECT_PREP_DIR}"/summary.tsv \
    "${LOCAL_DIRECT_PREP_DIR}"/missing_direct_snps.tsv \
    "${DSUB_DIRECT_PREP_URI}/" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

to_submit=()
already_done=()

for chrom in $(seq 1 22); do
    chrom_name="chr${chrom}"
    final_prefix="${DX_DIRECT_PFILE_DIR}/${chrom_name}"
    desired=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $2 }' "${summary}")
    available=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $3 }' "${summary}")
    missing=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $4 }' "${summary}")

    if [[ -z "${desired}" || -z "${available}" || -z "${missing}" ]]; then
        echo "ERROR: ${summary} has no complete row for ${chrom_name}" >&2
        exit 1
    fi

    if [[ -s "${final_prefix}.pgen" && -s "${final_prefix}.pvar" && -s "${final_prefix}.psam" ]]; then
        existing=$(grep -vc '^#' "${final_prefix}.pvar")
        if [[ "${existing}" -eq "${available}" ]]; then
            already_done+=("${chrom_name}")
            continue
        fi
        echo "  ${chrom_name}: direct pfile exists but has ${existing}/${available} present variants — rebuilding"
    fi
    to_submit+=("${chrom}")
done

echo "  Already on bucket: ${#already_done[@]} chrom(s) (${already_done[*]:-none})"
if [[ ${#to_submit[@]} -gt 0 ]]; then
    printf -v submit_names 'chr%s ' "${to_submit[@]}"
else
    submit_names="none"
fi
echo "  To submit:         ${#to_submit[@]} chrom(s) (${submit_names})"

if [[ ${#to_submit[@]} -gt 0 ]]; then
    tasks_tsv="${SCRIPT_DIR}/logs/dsub_direct_tasks_$(date +%Y%m%d_%H%M%S).tsv"
    mkdir -p "$(dirname "${tasks_tsv}")"
    {
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            '--env CHROM' \
            '--env DESIRED' '--env AVAILABLE' '--env MISSING' \
            '--input PLINK2' \
            '--input PGEN' '--input PVAR' '--input PSAM' \
            '--input EXTRACT' \
            '--output-recursive OUTDIR'
        for c in "${to_submit[@]}"; do
            chrom_name="chr${c}"
            desired=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $2 }' "${summary}")
            available=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $3 }' "${summary}")
            missing=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $4 }' "${summary}")
            printf 'chr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${c}" "${desired}" "${available}" "${missing}" \
                "${DSUB_PLINK2_GS}" \
                "${DX_WGS_PFILE_URI}/chr${c}.pgen" \
                "${DX_WGS_PFILE_URI}/chr${c}.pvar" \
                "${DX_WGS_PFILE_URI}/chr${c}.psam" \
                "${DSUB_DIRECT_PREP_URI}/chr${c}.extract.txt" \
                "${DX_DIRECT_PFILE_URI}/"
        done
    } > "${tasks_tsv}"

    echo "  Submitting dsub direct-SNP extraction tasks ..."
    dsub_out="${tasks_tsv%.tsv}.dsub.out"
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
        --name "sbayesrc-direct" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_direct_worker.sh" \
        --tasks "${tasks_tsv}" \
        --min-cores "${DSUB_MIN_CORES:-4}" \
        --min-ram "${DSUB_MIN_RAM:-32}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${DSUB_DISK_SIZE:-300}" \
        --wait \
        --summary 2>&1 | tee "${dsub_out}"
    dsub_rc=${PIPESTATUS[0]}

    dsub_job_id="$(awk '/^Launched job-id:/ {print $NF; exit}' "${dsub_out}")"
    if [[ -n "${dsub_job_id}" ]]; then
        expected=${#to_submit[@]}
        echo "  Polling dstat for job ${dsub_job_id} until all ${expected} tasks terminal ..."
        while true; do
            terminal_count=$(dstat --provider "${DSUB_PROVIDER}" \
                                   --project "${GOOGLE_PROJECT}" \
                                   --location "${DSUB_REGION}" \
                                   --jobs "${dsub_job_id}" \
                                   --users jupyter \
                                   --status '*' 2>/dev/null \
                             | awk 'NR>2 && /SUCCESS|FAILURE|CANCEL/ {c++} END {print c+0}')
            if (( terminal_count >= expected )); then
                echo "  ${terminal_count}/${expected} tasks terminal — proceeding to verification"
                break
            fi
            echo "  $(date -u +%H:%M:%SZ) ${terminal_count}/${expected} terminal — waiting 30s ..."
            sleep 30
        done
    fi

    if [[ "${dsub_rc}" -ne 0 ]]; then
        echo "ERROR: dsub direct-SNP extraction returned ${dsub_rc}" >&2
        echo "Check logs at ${DSUB_LOG_URI}/" >&2
        exit "${dsub_rc}"
    fi
fi

echo "  Verifying direct pfiles ..."
failed=()
for chrom in $(seq 1 22); do
    chrom_name="chr${chrom}"
    final_prefix="${DX_DIRECT_PFILE_DIR}/${chrom_name}"
    desired=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $2 }' "${summary}")
    available=$(awk -F'\t' -v c="${chrom_name}" '$1 == c { print $3 }' "${summary}")
    if [[ ! -s "${final_prefix}.pgen" || ! -s "${final_prefix}.pvar" || ! -s "${final_prefix}.psam" ]]; then
        failed+=("${chrom_name}:missing_output")
        continue
    fi
    observed=$(grep -vc '^#' "${final_prefix}.pvar")
    if [[ "${observed}" -ne "${available}" ]]; then
        failed+=("${chrom_name}:${observed}/${available}")
    fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "ERROR: direct pfile verification failed: ${failed[*]}" >&2
    exit 1
fi

echo "  Direct per-chromosome pfiles are in ${DX_DIRECT_PFILE_DIR}"
