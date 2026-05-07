#!/bin/bash
# upload_sbayesrc_ids.sh — Upload per-chromosome SBayesRC ID + idmap files to
# the workspace bucket on GCS (requester-pays).
#
# Expects env vars:
#   LOCAL_SBAYESRC_ID_DIR   local source dir (chr{N}.extract.txt, chr{N}.idmap.txt)
#   DX_SBAYESRC_ID_DIR      gs:// destination prefix
#   GOOGLE_PROJECT          billing project for requester-pays reads
#
# Idempotent: skips upload if the destination object already exists.
# Never downloads bytes; only `gcloud storage ls` for metadata checks.

set -euo pipefail

uploaded=0
skipped=0

for chrom in $(seq 1 22); do
    for kind in extract idmap; do
        filename="chr${chrom}.${kind}.txt"
        local_path="${LOCAL_SBAYESRC_ID_DIR}/${filename}"
        remote_path="${DX_SBAYESRC_ID_DIR}/${filename}"

        if ! [[ -f "${local_path}" ]]; then
            echo "ERROR: ${local_path} does not exist. Run step 1 first."
            exit 1
        fi

        if gcloud storage ls --billing-project="${GOOGLE_PROJECT}" \
                "${remote_path}" &>/dev/null; then
            skipped=$((skipped + 1))
        else
            gcloud storage cp --billing-project="${GOOGLE_PROJECT}" \
                "${local_path}" "${remote_path}"
            uploaded=$((uploaded + 1))
        fi
    done
done

echo "Uploaded: ${uploaded}, Skipped (already exist): ${skipped}"
