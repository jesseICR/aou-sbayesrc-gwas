#!/bin/bash
# aou_downstream_env.sh - Shared AoU environment defaults for downstream GWAS workflows.
#
# Source this from downstream scripts after SCRIPT_DIR has been set.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "ERROR: source this file from a workflow script; do not execute it directly." >&2
    exit 1
fi

if [[ -z "${SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    export SCRIPT_DIR
fi

if [[ -z "${GOOGLE_PROJECT:-}" ]]; then
    GOOGLE_PROJECT="$(gcloud config get-value project 2>/dev/null || true)"
    if [[ "${GOOGLE_PROJECT}" == "(unset)" ]]; then
        GOOGLE_PROJECT=""
    fi
    export GOOGLE_PROJECT
fi
: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set - are you running inside an AoU Verily Jupyter session?}"

export AOU_DATA_VERSION="${AOU_DATA_VERSION:-v9}"
case "${AOU_DATA_VERSION}" in
    v8|v9) ;;
    *)
        echo "ERROR: AOU_DATA_VERSION must be v8 or v9; got '${AOU_DATA_VERSION}'." >&2
        return 1
        ;;
esac

infer_downstream_cdr_project() {
    if [[ -n "${WORKSPACE_CDR:-}" && "${WORKSPACE_CDR}" == *.* ]]; then
        printf '%s\n' "${WORKSPACE_CDR%%.*}"
    else
        printf '%s\n' "wb-silky-artichoke-2408"
    fi
}

export AOU_CDR_PROJECT="${AOU_CDR_PROJECT:-$(infer_downstream_cdr_project)}"
case "${AOU_DATA_VERSION}" in
    v9)
        export AOU_CDR_DATASET="${AOU_CDR_DATASET:-C2025Q4R6}"
        ;;
    v8)
        if [[ -z "${AOU_CDR_DATASET:-}" && -n "${WORKSPACE_CDR:-}" && "${WORKSPACE_CDR}" == *.* ]]; then
            export AOU_CDR_DATASET="${WORKSPACE_CDR##*.}"
        fi
        export AOU_CDR_DATASET="${AOU_CDR_DATASET:-C2024Q3R9}"
        ;;
esac

export AOU_TARGET_WORKSPACE_CDR="${AOU_CDR_PROJECT}.${AOU_CDR_DATASET}"
if [[ "${AOU_STRICT_WORKSPACE_CDR:-0}" == "1" && -n "${WORKSPACE_CDR:-}" && "${WORKSPACE_CDR}" != "${AOU_TARGET_WORKSPACE_CDR}" ]]; then
    echo "ERROR: WORKSPACE_CDR=${WORKSPACE_CDR}, but AOU_DATA_VERSION=${AOU_DATA_VERSION} expects ${AOU_TARGET_WORKSPACE_CDR}." >&2
    return 1
fi
export WORKSPACE_CDR_ORIGINAL="${WORKSPACE_CDR:-}"
export WORKSPACE_CDR="${AOU_TARGET_WORKSPACE_CDR}"

# Keep off-cycle dataset choices overrideable. ETM task tables are part of the
# v9 main CDR, while the older v8 workflow used the off-cycle ETM dataset.
export WORKSPACE_MHWB_CDR="${WORKSPACE_MHWB_CDR:-${AOU_CDR_PROJECT}.C_V8_R2_offcycle_mhwb}"
if [[ -z "${WORKSPACE_ETM_CDR:-}" ]]; then
    case "${AOU_DATA_VERSION}" in
        v9)
            export WORKSPACE_ETM_CDR="${WORKSPACE_CDR}"
            ;;
        v8)
            export WORKSPACE_ETM_CDR="${AOU_CDR_PROJECT}.C_V8_R2_offcycle_etm"
            ;;
    esac
else
    export WORKSPACE_ETM_CDR
fi

export WORKSPACE_BUCKET_MOUNT="${WORKSPACE_BUCKET_MOUNT:-/home/jupyter/workspace/workspace-bucket}"

is_downstream_workspace_bucket_fuse_mounted() {
    mount | awk -v target="${WORKSPACE_BUCKET_MOUNT}" '
        $2 == "on" && $3 == target && $4 == "type" && $5 == "fuse.gcsfuse" { found = 1 }
        END { exit found ? 0 : 1 }
    '
}

ensure_downstream_workspace_bucket_mount() {
    if is_downstream_workspace_bucket_fuse_mounted && [[ -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
        return 0
    fi
    echo "Workspace bucket is not mounted at ${WORKSPACE_BUCKET_MOUNT}; attempting Workbench mount ..."
    if ! command -v wb >/dev/null 2>&1; then
        echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is not a writable gcsfuse mount, and 'wb' is unavailable." >&2
        return 1
    fi
    if ! wb resource mount --allow-other; then
        echo "ERROR: 'wb resource mount --allow-other' failed." >&2
        return 1
    fi
    if ! is_downstream_workspace_bucket_fuse_mounted || [[ ! -w "${WORKSPACE_BUCKET_MOUNT}" ]]; then
        echo "ERROR: ${WORKSPACE_BUCKET_MOUNT} is still not a writable gcsfuse mount." >&2
        return 1
    fi
}

ensure_downstream_workspace_bucket_mount
WORKSPACE_BUCKET_URI="gs://$(mount | awk -v target="${WORKSPACE_BUCKET_MOUNT}" '$2 == "on" && $3 == target {print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    return 1
fi
export WORKSPACE_BUCKET_URI

export SBAYESRC_OUTPUT_PREFIX="${SBAYESRC_OUTPUT_PREFIX:-sbayesrc_genotypes}"
export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/${SBAYESRC_OUTPUT_PREFIX}"
export DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
export DX_REGENIE_OUTPUT_DIR="${DX_OUTPUT_DIR}/regenie_output"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"
export DX_OUTPUT_URI="${WORKSPACE_BUCKET_URI}/${SBAYESRC_OUTPUT_PREFIX}"
export DX_REGENIE_INPUT_URI="${DX_OUTPUT_URI}/regenie_input"
export DX_REGENIE_OUTPUT_URI="${DX_OUTPUT_URI}/regenie_output"
