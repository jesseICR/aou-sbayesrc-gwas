#!/bin/bash
# admixture_run_projection.sh - Run ADMIXTURE projection and concatenate results.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set}"
: "${DX_STATGEN_DIR:?DX_STATGEN_DIR not set}"
: "${DX_STATGEN_URI:?DX_STATGEN_URI not set}"
: "${DX_ADMIXTURE_SCRAP_DIR:?DX_ADMIXTURE_SCRAP_DIR not set}"
: "${DX_ADMIXTURE_SCRAP_URI:?DX_ADMIXTURE_SCRAP_URI not set}"
: "${DX_ADMIXTURE_BATCH_DIR:?DX_ADMIXTURE_BATCH_DIR not set}"
: "${DX_ADMIXTURE_BATCH_URI:?DX_ADMIXTURE_BATCH_URI not set}"
: "${DX_ADMIXTURE_Q_DIR:?DX_ADMIXTURE_Q_DIR not set}"
: "${DX_ADMIXTURE_Q_URI:?DX_ADMIXTURE_Q_URI not set}"
: "${LOCAL_ADMIXTURE_DIR:?LOCAL_ADMIXTURE_DIR not set}"
: "${ADMIXTURE_K:?ADMIXTURE_K not set}"
: "${DSUB_PROVIDER:?DSUB_PROVIDER not set}"
: "${DSUB_REGION:?DSUB_REGION not set}"
: "${DSUB_NETWORK:?DSUB_NETWORK not set}"
: "${DSUB_SUBNETWORK:?DSUB_SUBNETWORK not set}"
: "${DSUB_PET_SA:?DSUB_PET_SA not set}"
: "${DSUB_IMAGE:?DSUB_IMAGE not set}"
: "${DSUB_LOG_URI:?DSUB_LOG_URI not set}"

mkdir -p "${LOCAL_ADMIXTURE_DIR}" "${DX_STATGEN_DIR}" "${DX_ADMIXTURE_Q_DIR}"

ADMIXTURE_PROJECT_DSUB_MIN_CORES="${ADMIXTURE_PROJECT_DSUB_MIN_CORES:-2}"
ADMIXTURE_PROJECT_DSUB_MIN_RAM="${ADMIXTURE_PROJECT_DSUB_MIN_RAM:-8}"
ADMIXTURE_PROJECT_DSUB_DISK_SIZE="${ADMIXTURE_PROJECT_DSUB_DISK_SIZE:-100}"
ADMIXTURE_PROJECT_DSUB_DISK_TYPE="${ADMIXTURE_PROJECT_DSUB_DISK_TYPE:-pd-ssd}"
ADMIXTURE_CONCAT_DSUB_MIN_CORES="${ADMIXTURE_CONCAT_DSUB_MIN_CORES:-2}"
ADMIXTURE_CONCAT_DSUB_MIN_RAM="${ADMIXTURE_CONCAT_DSUB_MIN_RAM:-8}"
ADMIXTURE_CONCAT_DSUB_DISK_SIZE="${ADMIXTURE_CONCAT_DSUB_DISK_SIZE:-50}"
ADMIXTURE_CONCAT_DSUB_DISK_TYPE="${ADMIXTURE_CONCAT_DSUB_DISK_TYPE:-pd-ssd}"

manifest="${DX_ADMIXTURE_BATCH_DIR}/batch_manifest.tsv"
split_params="${DX_ADMIXTURE_BATCH_DIR}/admixture_split.params.tsv"
aligned_fam="${DX_ADMIXTURE_SCRAP_DIR}/aou_admixture_aligned.fam"
admixture_bin="${DX_ADMIXTURE_SCRAP_DIR}/_prep/admixture"

if [[ ! -s "${manifest}" || ! -s "${split_params}" || ! -s "${aligned_fam}" || ! -s "${admixture_bin}" ]]; then
    echo "ERROR: missing ADMIXTURE split/prep outputs required for projection" >&2
    exit 1
fi

total_samples=$(wc -l < "${aligned_fam}")
n_batches=$(awk 'NR > 1 {n++} END {print n + 0}' "${manifest}")
split_sha256=$(sha256sum "${split_params}" | awk '{print $1}')
manifest_sha256=$(sha256sum "${manifest}" | awk '{print $1}')
admixture_sha256=$(sha256sum "${admixture_bin}" | awk '{print $1}')

desired_params="${LOCAL_ADMIXTURE_DIR}/admixture_projection.desired_params.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'admixture_k\t%s\n' "${ADMIXTURE_K}"
    printf 'samples\t%s\n' "${total_samples}"
    printf 'batches\t%s\n' "${n_batches}"
    printf 'split_params_sha256\t%s\n' "${split_sha256}"
    printf 'batch_manifest_sha256\t%s\n' "${manifest_sha256}"
    printf 'admixture_binary_sha256\t%s\n' "${admixture_sha256}"
} > "${desired_params}"

final_tsv="${DX_STATGEN_DIR}/aou_admixture_k${ADMIXTURE_K}.tsv"
final_summary="${DX_STATGEN_DIR}/aou_admixture_k${ADMIXTURE_K}.summary.tsv"
final_params="${DX_STATGEN_DIR}/aou_admixture_k${ADMIXTURE_K}.params.tsv"
if [[ -s "${final_tsv}" && -s "${final_summary}" && -s "${final_params}" ]]; then
    if diff -q "${desired_params}" "${final_params}" >/dev/null 2>&1; then
        final_lines=$(wc -l < "${final_tsv}")
        if [[ "${final_lines}" -eq $((total_samples + 1)) ]]; then
            echo "  ADMIXTURE K=${ADMIXTURE_K} result already exists (${total_samples} samples) — skipping"
            exit 0
        fi
    fi
    echo "  Final ADMIXTURE output exists but params/counts do not match — rebuilding projection outputs"
fi

q_params="${DX_ADMIXTURE_Q_DIR}/admixture_projection.params.tsv"
q_params_match=0
if [[ -s "${q_params}" ]] && diff -q "${desired_params}" "${q_params}" >/dev/null 2>&1; then
    q_params_match=1
elif [[ -s "${q_params}" ]]; then
    echo "  Existing ADMIXTURE Q params do not match — stale batch Q files will be overwritten"
    gcloud storage rm "${DX_ADMIXTURE_Q_URI}/batch_*.${ADMIXTURE_K}.Q" \
        --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true
fi

tasks_tsv="${SCRIPT_DIR}/logs/dsub_admixture_projection_tasks_$(date +%Y%m%d_%H%M%S).tsv"
mkdir -p "$(dirname "${tasks_tsv}")"
{
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '--env BATCH' '--env ADMIXTURE_K' '--env EXPECTED_SAMPLES' '--env EXPECTED_VARIANTS' \
        '--input ADMIXTURE' '--input BED' '--input BIM' '--input FAM' '--input PFILE' '--output Q'
    while IFS=$'\t' read -r batch start end samples variants p_rows; do
        [[ "${batch}" == "batch" ]] && continue
        q_file="${DX_ADMIXTURE_Q_DIR}/${batch}.${ADMIXTURE_K}.Q"
        submit=1
        if [[ "${q_params_match}" -eq 1 && -s "${q_file}" ]]; then
            q_lines=$(wc -l < "${q_file}")
            if [[ "${q_lines}" -eq "${samples}" ]]; then
                submit=0
            fi
        fi
        if [[ "${submit}" -eq 1 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${batch}" "${ADMIXTURE_K}" "${samples}" "${variants}" \
                "${DX_ADMIXTURE_SCRAP_URI}/_prep/admixture" \
                "${DX_ADMIXTURE_BATCH_URI}/${batch}.bed" \
                "${DX_ADMIXTURE_BATCH_URI}/${batch}.bim" \
                "${DX_ADMIXTURE_BATCH_URI}/${batch}.fam" \
                "${DX_ADMIXTURE_BATCH_URI}/${batch}.${ADMIXTURE_K}.P.in" \
                "${DX_ADMIXTURE_Q_URI}/${batch}.${ADMIXTURE_K}.Q"
        fi
    done < "${manifest}"
} > "${tasks_tsv}"

to_submit=$(( $(wc -l < "${tasks_tsv}") - 1 ))
echo "  ADMIXTURE projection batches: ${n_batches}; to submit: ${to_submit}"

if [[ "${to_submit}" -gt 0 ]]; then
    projection_log="${tasks_tsv%.tsv}.dsub.out"
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
        --name "sbayesrc-admixture-project" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_admixture_project_worker.sh" \
        --tasks "${tasks_tsv}" \
        --min-cores "${ADMIXTURE_PROJECT_DSUB_MIN_CORES}" \
        --min-ram "${ADMIXTURE_PROJECT_DSUB_MIN_RAM}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${ADMIXTURE_PROJECT_DSUB_DISK_SIZE}" \
        --disk-type "${ADMIXTURE_PROJECT_DSUB_DISK_TYPE}" \
        --wait \
        --summary 2>&1 | tee "${projection_log}"
    dsub_rc=${PIPESTATUS[0]}

    dsub_job_id="$(awk '/^Launched job-id:/ {print $NF; exit}' "${projection_log}")"
    if [[ -n "${dsub_job_id}" ]]; then
        echo "  Polling dstat for job ${dsub_job_id} until all ${to_submit} tasks terminal ..."
        while true; do
            terminal_count=$(dstat --provider "${DSUB_PROVIDER}" \
                                   --project "${GOOGLE_PROJECT}" \
                                   --location "${DSUB_REGION}" \
                                   --jobs "${dsub_job_id}" \
                                   --users jupyter \
                                   --status '*' 2>/dev/null \
                             | awk 'NR>2 && /SUCCESS|FAILURE|CANCEL/ {c++} END {print c+0}')
            if (( terminal_count >= to_submit )); then
                echo "  ${terminal_count}/${to_submit} tasks terminal — proceeding to verification"
                break
            fi
            echo "  $(date -u +%H:%M:%SZ) ${terminal_count}/${to_submit} terminal — waiting 30s ..."
            sleep 30
        done
    fi

    if [[ "${dsub_rc}" -ne 0 ]]; then
        echo "ERROR: dsub ADMIXTURE projection returned ${dsub_rc}" >&2
        echo "Check logs at ${DSUB_LOG_URI}/" >&2
        exit "${dsub_rc}"
    fi
fi

echo "  Verifying ADMIXTURE Q files ..."
while IFS=$'\t' read -r batch start end samples variants p_rows; do
    [[ "${batch}" == "batch" ]] && continue
    q_file="${DX_ADMIXTURE_Q_DIR}/${batch}.${ADMIXTURE_K}.Q"
    if [[ ! -s "${q_file}" ]]; then
        echo "ERROR: missing ${q_file}" >&2
        exit 1
    fi
    q_lines=$(wc -l < "${q_file}")
    if [[ "${q_lines}" -ne "${samples}" ]]; then
        echo "ERROR: ${q_file} has ${q_lines}/${samples} Q rows" >&2
        exit 1
    fi
done < "${manifest}"

gcloud storage cp "${desired_params}" "${DX_ADMIXTURE_Q_URI}/admixture_projection.params.tsv" \
    --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

concat_log="${SCRIPT_DIR}/logs/dsub_admixture_concat_$(date +%Y%m%d_%H%M%S).dsub.out"
echo "  Submitting dsub concat job for final ADMIXTURE TSV ..."
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
    --name "sbayesrc-admixture-concat" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_admixture_concat_worker.sh" \
    --env ADMIXTURE_K="${ADMIXTURE_K}" \
    --env EXPECTED_SAMPLES="${total_samples}" \
    --env EXPECTED_BATCHES="${n_batches}" \
    --input FAM="${DX_ADMIXTURE_SCRAP_URI}/aou_admixture_aligned.fam" \
    --input MANIFEST="${DX_ADMIXTURE_BATCH_URI}/batch_manifest.tsv" \
    --input PARAMS="${DX_ADMIXTURE_Q_URI}/admixture_projection.params.tsv" \
    --input-recursive QDIR="${DX_ADMIXTURE_Q_URI}/" \
    --output-recursive OUTDIR="${DX_STATGEN_URI}/" \
    --min-cores "${ADMIXTURE_CONCAT_DSUB_MIN_CORES}" \
    --min-ram "${ADMIXTURE_CONCAT_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
    --disk-size "${ADMIXTURE_CONCAT_DSUB_DISK_SIZE}" \
    --disk-type "${ADMIXTURE_CONCAT_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${concat_log}"
dsub_rc=${PIPESTATUS[0]}
if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: dsub ADMIXTURE concat returned ${dsub_rc}" >&2
    echo "Check logs at ${DSUB_LOG_URI}/" >&2
    exit "${dsub_rc}"
fi

echo "  Verifying final ADMIXTURE TSV ..."
for attempt in $(seq 1 60); do
    if [[ -s "${final_tsv}" && -s "${final_summary}" && -s "${final_params}" ]]; then
        final_lines=$(wc -l < "${final_tsv}")
        if [[ "${final_lines}" -eq $((total_samples + 1)) ]] && diff -q "${desired_params}" "${final_params}" >/dev/null 2>&1; then
            echo "  Done: ${final_tsv} (${total_samples} samples)"
            exit 0
        fi
        echo "  Final TSV visible but counts/params do not match yet; waiting ..."
    else
        echo "  Waiting for final ADMIXTURE TSV to appear (${attempt}/60) ..."
    fi
    sleep 10
done

echo "ERROR: final ADMIXTURE TSV was not visible with expected counts after dsub completed" >&2
exit 1
