#!/bin/bash
# wgs_extract_variants.sh — SBayesRC variant extraction via dsub fan-out.
#
# Submits one dsub task per chromosome (parallel across ~21 Google Batch
# workers in us-central1). Each worker runs dsub_extract_worker.sh:
#   Pass 1 (split + ID) → Pass 1.5 (normalize pvar) → Pass 2 (extract)
#   → Pass 3 (rsid remap). Outputs delocalized to DX_WGS_PFILE_URI.
#
# Why dsub-fanout instead of running plink2 sequentially on this Jupyter pod:
#   - 22× compute in parallel instead of 1 chrom at a time
#   - Workers run in us-central1 (same region as workspace bucket), avoiding
#     the cross-region upload ceiling that capped this pod at ~7 MB/s for the
#     32 GB chr2 pgen (~75 min just to upload). Intra-region same-bucket
#     writes are essentially line-rate.
#   - Frees the Jupyter pod for interactive work while extracts run.
#
# Idempotent: skips a chromosome if chrN.pgen already exists on the bucket
# (cheap gcsfuse metadata lookup; no submission cost).
#
# Required env vars (set by get_genotypes.sh):
#   GOOGLE_PROJECT
#   PLINK2                   local plink2 binary (staged to bucket once)
#   SCRIPT_DIR               repo root
#   LOCAL_SBAYESRC_ID_DIR    local dir with chr{N}.{extract,idmap}.txt
#   DX_WGS_PFILE_DIR         workspace-bucket-mount path (idempotency check)
#   DX_WGS_PFILE_URI         gs:// destination for delocalized pfiles
#   AOU_PGEN_GS_DIR          gs://vwb-aou-datasets-controlled/v8/.../pgen
#   DSUB_PROVIDER            "google-batch"
#   DSUB_REGION              "us-central1"
#   DSUB_NETWORK             projects/$PROJECT/global/networks/network
#   DSUB_SUBNETWORK          projects/$PROJECT/regions/$REGION/subnetworks/subnetwork
#   DSUB_PET_SA              this pod's pet service account (for --service-account)
#   DSUB_IMAGE               Docker image (default: ubuntu2204 marketplace image)
#   DSUB_PLINK2_GS           gs:// path of staged plink2 binary
#   DSUB_SBAYESRC_ID_URI     gs:// path of staged chr{N}.{extract,idmap}.txt
#   DSUB_LOG_URI             gs:// path for dsub task logs
#
# Optional:
#   SBAYESRC_TEST_CHROM      e.g. "21" — submit only this chromosome
#   DSUB_MIN_CORES           default 4
#   DSUB_MIN_RAM             default 32 (GB)
#   DSUB_DISK_SIZE           default 300 (GB) — peak chr1 needs ~250 GB
#   DSUB_BOOT_DISK_SIZE      default 50 (GB)

set -euo pipefail

# ---------------------------------------------------------------------------
# Build chrom list (idempotency: skip ones already on bucket)
# ---------------------------------------------------------------------------
if [[ -n "${SBAYESRC_TEST_CHROM:-}" ]]; then
    requested_chroms=("${SBAYESRC_TEST_CHROM}")
else
    requested_chroms=($(seq 1 22))
fi

to_submit=()
already_done=()
for c in "${requested_chroms[@]}"; do
    if [[ -f "${DX_WGS_PFILE_DIR}/chr${c}.pgen" ]]; then
        already_done+=("chr${c}")
    else
        to_submit+=("${c}")
    fi
done

echo "Already on bucket: ${#already_done[@]} chrom(s) (${already_done[*]:-none})"
echo "To submit:         ${#to_submit[@]} chrom(s) ($(printf 'chr%s ' "${to_submit[@]}" 2>/dev/null))"

if [[ ${#to_submit[@]} -eq 0 ]]; then
    echo "All requested chromosomes already extracted. Nothing to submit."
    # Still write combined summary below.
else
    # -----------------------------------------------------------------------
    # Stage worker binary + reference files to workspace bucket (idempotent)
    # -----------------------------------------------------------------------
    # The pvar-allele normalization logic from normalize_pvar_alleles.py is
    # inlined as awk in dsub_extract_worker.sh, so we don't stage the .py.
    echo ""
    echo "Staging worker binary + IDs to workspace bucket ..."
    gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" \
        --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1
    gcloud storage cp \
        "${LOCAL_SBAYESRC_ID_DIR}"/chr*.extract.txt \
        "${LOCAL_SBAYESRC_ID_DIR}"/chr*.idmap.txt \
        "${DSUB_SBAYESRC_ID_URI}/" \
        --billing-project="${GOOGLE_PROJECT}" 2>&1 | tail -1

    # -----------------------------------------------------------------------
    # Build --tasks TSV
    # -----------------------------------------------------------------------
    tasks_tsv="${SCRIPT_DIR}/logs/dsub_tasks_$(date +%Y%m%d_%H%M%S).tsv"
    mkdir -p "$(dirname "${tasks_tsv}")"
    {
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            '--env CHROM' \
            '--input PLINK2' \
            '--input PGEN' '--input PVAR' '--input PSAM' \
            '--input EXTRACT' '--input IDMAP' \
            '--output-recursive OUTDIR'
        for c in "${to_submit[@]}"; do
            printf 'chr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${c}" \
                "${DSUB_PLINK2_GS}" \
                "${AOU_PGEN_GS_DIR}/acaf_threshold.chr${c}.pgen" \
                "${AOU_PGEN_GS_DIR}/acaf_threshold.chr${c}.pvar" \
                "${AOU_PGEN_GS_DIR}/acaf_threshold.chr${c}.psam" \
                "${DSUB_SBAYESRC_ID_URI}/chr${c}.extract.txt" \
                "${DSUB_SBAYESRC_ID_URI}/chr${c}.idmap.txt" \
                "${DX_WGS_PFILE_URI}/"
        done
    } > "${tasks_tsv}"

    echo ""
    echo "Tasks TSV (${tasks_tsv}):"
    head -1 "${tasks_tsv}" | tr '\t' '\n' | nl | sed "s/^/  /"
    echo "  → ${#to_submit[@]} task rows"

    # -----------------------------------------------------------------------
    # Submit dsub
    # -----------------------------------------------------------------------
    echo ""
    echo "Submitting dsub: provider=${DSUB_PROVIDER}, region=${DSUB_REGION}, image=${DSUB_IMAGE}"
    echo "  cores=${DSUB_MIN_CORES:-4} ram=${DSUB_MIN_RAM:-32}G boot-disk=${DSUB_BOOT_DISK_SIZE:-50}G data-disk=${DSUB_DISK_SIZE:-300}G"
    echo "  logging=${DSUB_LOG_URI}"

    # Capture dsub stdout to extract the job-id. `dsub --wait` exits on the
    # FIRST task failure; the remaining tasks keep running on Batch but the
    # orchestrator returns prematurely. We then need to poll dstat ourselves
    # until every task reaches a terminal state, otherwise the bucket-output
    # verification below would false-fail on tasks that simply haven't
    # delocalized yet.
    dsub_out="${tasks_tsv%.tsv}.dsub.out"
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
        --name "sbayesrc-extract" \
        --image "${DSUB_IMAGE}" \
        --script "${SCRIPT_DIR}/dsub_extract_worker.sh" \
        --tasks "${tasks_tsv}" \
        --min-cores "${DSUB_MIN_CORES:-4}" \
        --min-ram "${DSUB_MIN_RAM:-32}" \
        --boot-disk-size "${DSUB_BOOT_DISK_SIZE:-50}" \
        --disk-size "${DSUB_DISK_SIZE:-300}" \
        --wait \
        --summary 2>&1 | tee "${dsub_out}"
    dsub_rc=${PIPESTATUS[0]}

    # -----------------------------------------------------------------------
    # If dsub exited early on a failed task, wait for the remaining tasks to
    # reach terminal state on Batch before declaring overall success/failure.
    # -----------------------------------------------------------------------
    dsub_job_id="$(awk '/^Launched job-id:/ {print $NF; exit}' "${dsub_out}")"
    if [[ -n "${dsub_job_id}" ]]; then
        expected=${#to_submit[@]}
        echo ""
        echo "Polling dstat for job ${dsub_job_id} until all ${expected} tasks terminal ..."
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
    else
        echo "WARN: could not parse dsub job-id from '${dsub_out}'; skipping post-wait poll"
    fi

    # -----------------------------------------------------------------------
    # Per-chrom verification (now that every task is actually terminal)
    # -----------------------------------------------------------------------
    echo ""
    echo "=== Verifying outputs on workspace bucket ==="
    succeeded=()
    failed=()
    for c in "${to_submit[@]}"; do
        if [[ -f "${DX_WGS_PFILE_DIR}/chr${c}.pgen" \
              && -f "${DX_WGS_PFILE_DIR}/chr${c}.pvar" \
              && -f "${DX_WGS_PFILE_DIR}/chr${c}.psam" \
              && -f "${DX_WGS_PFILE_DIR}/chr${c}.summary.tsv" ]]; then
            succeeded+=("chr${c}")
        else
            failed+=("chr${c}")
        fi
    done

    echo "Succeeded: ${#succeeded[@]}/${#to_submit[@]} (${succeeded[*]:-none})"
    if [[ ${#failed[@]} -gt 0 ]]; then
        echo "Failed:    ${#failed[@]}/${#to_submit[@]} (${failed[*]})"
        echo "Check dsub logs at: ${DSUB_LOG_URI}/"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Combined summary (always — useful even on idempotent re-runs)
# ---------------------------------------------------------------------------
echo ""
echo "=== SBayesRC extraction summary (chromosomes present on bucket) ==="
combined="${SCRIPT_DIR}/logs/sbayesrc_extract_summary.tsv"
{
    header_written=0
    for c in $(seq 1 22); do
        s="${DX_WGS_PFILE_DIR}/chr${c}.summary.tsv"
        [[ -f "${s}" ]] || continue
        if [[ ${header_written} -eq 0 ]]; then
            head -1 "${s}"
            header_written=1
        fi
        tail -n +2 "${s}"
    done
} | tee "${combined}"

echo ""
echo "Combined summary at ${combined}"
