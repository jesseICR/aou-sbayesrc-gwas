#!/bin/bash
# Run a continuous REGENIE GWAS for the fold-safe no-teacher
# GradCPT/Flanker factor18 calibrated proxy phenotype.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR

usage() {
    cat <<'EOF'
Usage: bash run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas.sh [OPTIONS]

Creates a fold-safe 3-variable no-teacher calibrated phenotype trained on the
18k GradCPT/Flanker missing-pattern-aware factor target, then runs REGENIE.

Phenotype:
  gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z

Calibration predictors:
  ses_ea_proxy_z
  gradcpt_flanker_finetuned_ea_proxy_z
  gradcpt_flanker_direct_xgb_proxy_z

GWAS covariates:
  sex_c + PC1_AVG ... PC10_AVG

Options:
  --smoke             Run chr22 only and use a chr22_smoke output directory.
  --chroms LIST       Chromosomes for REGENIE Step 2, e.g. 22, 1,2,3, or 1-22.
  --preflight-only    Validate inputs and write metadata; do not submit REGENIE.
  --apply-rint        Apply rank-inverse normal transform (default).
  --no-apply-rint     Disable rank-inverse normal transform.
  --force             Allow rebuilding an existing input package and reusing an
                      existing output directory through the generic runner.
  -h, --help          Show this help.
EOF
}

SMOKE=0
PREFLIGHT_ONLY=0
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --smoke)
            SMOKE=1
            shift
            ;;
        --chroms)
            export REGENIE_CHROMS="${2:?--chroms requires a value}"
            shift 2
            ;;
        --preflight-only)
            PREFLIGHT_ONLY=1
            shift
            ;;
        --apply-rint)
            export REGENIE_APPLY_RINT=1
            shift
            ;;
        --no-apply-rint)
            export REGENIE_APPLY_RINT=0
            shift
            ;;
        --force)
            FORCE=1
            shift
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

: "${GOOGLE_PROJECT:?GOOGLE_PROJECT not set — are you running inside an AoU Verily Jupyter session?}"

export WORKSPACE_BUCKET_MOUNT="/home/jupyter/workspace/workspace-bucket"
WORKSPACE_BUCKET_URI="gs://$(mount | awk '/ on \/home\/jupyter\/workspace\/workspace-bucket /{print $1; exit}')"
if [[ "${WORKSPACE_BUCKET_URI}" == "gs://" ]]; then
    echo "ERROR: could not derive workspace bucket URI from mount table." >&2
    exit 1
fi
export WORKSPACE_BUCKET_URI

export DX_OUTPUT_DIR="${WORKSPACE_BUCKET_MOUNT}/sbayesrc_genotypes"
export DX_GWAS_GENOTYPES_DIR="${DX_OUTPUT_DIR}/gwas_genotypes"
export DX_GWAS_STEP1_BFILE_DIR="${DX_GWAS_GENOTYPES_DIR}/step1_direct"
export DX_GWAS_STEP2_PFILE_DIR="${DX_GWAS_GENOTYPES_DIR}/step2_wgs_pfiles"
export DX_REGENIE_INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input"
export DX_REGENIE_OUTPUT_DIR="${DX_OUTPUT_DIR}/regenie_output"
export DX_LOGS_DIR="${DX_OUTPUT_DIR}/logs"
export DX_GWAS_GENOTYPES_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/gwas_genotypes"
export DX_GWAS_STEP1_BFILE_URI="${DX_GWAS_GENOTYPES_URI}/step1_direct"
export DX_GWAS_STEP2_PFILE_URI="${DX_GWAS_GENOTYPES_URI}/step2_wgs_pfiles"
export DX_REGENIE_INPUT_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/regenie_input"
export DX_REGENIE_OUTPUT_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/regenie_output"

export FACTOR18_GWAS_SOURCE_NAME="${FACTOR18_GWAS_SOURCE_NAME:-gradcpt_flanker_direct_xgb_proxy_ses_ea_proxy_v2_kinholdout}"
export FACTOR18_GWAS_INPUT_NAME="${FACTOR18_GWAS_INPUT_NAME:-gradcpt_flanker_factor18_no_teacher_calibrated_proxy_ses_ea_proxy_v2_kinholdout}"
BASE_OUTPUT_NAME="${FACTOR18_GWAS_OUTPUT_NAME:-gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas}"
if [[ "${SMOKE}" -eq 1 ]]; then
    export REGENIE_CHROMS="22"
    export FACTOR18_GWAS_OUTPUT_NAME="${BASE_OUTPUT_NAME}_chr22_smoke"
else
    export REGENIE_CHROMS="${REGENIE_CHROMS:-1-22}"
    export FACTOR18_GWAS_OUTPUT_NAME="${BASE_OUTPUT_NAME}"
fi

export LOCAL_REGENIE_DIR="${SCRIPT_DIR}/data/regenie"
export REGENIE_APPLY_RINT="${REGENIE_APPLY_RINT:-1}"
export REGENIE_STEP1_BLOCK_SIZE="${REGENIE_STEP1_BLOCK_SIZE:-1000}"
export REGENIE_STEP2_BLOCK_SIZE="${REGENIE_STEP2_BLOCK_SIZE:-200}"
export REGENIE_PHENO_COL="gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z"
export REGENIE_COVAR_COLS="sex_c,PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG,PC5_AVG,PC6_AVG,PC7_AVG,PC8_AVG,PC9_AVG,PC10_AVG"

export DSUB_PROVIDER="${DSUB_PROVIDER:-google-batch}"
export DSUB_REGION="${DSUB_REGION:-us-central1}"
export DSUB_NETWORK="projects/${GOOGLE_PROJECT}/global/networks/network"
export DSUB_SUBNETWORK="projects/${GOOGLE_PROJECT}/regions/${DSUB_REGION}/subnetworks/subnetwork"
export DSUB_IMAGE="${DSUB_IMAGE:-marketplace.gcr.io/google/ubuntu2204}"
DSUB_PET_SA="$(gcloud config get-value account 2>/dev/null || true)"
if [[ -z "${DSUB_PET_SA}" ]]; then
    echo "ERROR: could not determine the pod's pet service account via 'gcloud config get-value account'." >&2
    exit 1
fi
export DSUB_PET_SA
export DSUB_BOOT_DISK_SIZE="${DSUB_BOOT_DISK_SIZE:-50}"
export REGENIE_STEP1_DSUB_MIN_CORES="${REGENIE_STEP1_DSUB_MIN_CORES:-16}"
export REGENIE_STEP1_DSUB_MIN_RAM="${REGENIE_STEP1_DSUB_MIN_RAM:-64}"
export REGENIE_STEP1_DSUB_DISK_SIZE="${REGENIE_STEP1_DSUB_DISK_SIZE:-300}"
export REGENIE_STEP1_DSUB_DISK_TYPE="${REGENIE_STEP1_DSUB_DISK_TYPE:-pd-ssd}"
export REGENIE_STEP2_DSUB_MIN_CORES="${REGENIE_STEP2_DSUB_MIN_CORES:-8}"
export REGENIE_STEP2_DSUB_MIN_RAM="${REGENIE_STEP2_DSUB_MIN_RAM:-32}"
export REGENIE_STEP2_DSUB_DISK_SIZE="${REGENIE_STEP2_DSUB_DISK_SIZE:-300}"
export REGENIE_STEP2_DSUB_DISK_TYPE="${REGENIE_STEP2_DSUB_DISK_TYPE:-pd-ssd}"
export DSUB_BIN_URI="${WORKSPACE_BUCKET_URI}/bin"
export DSUB_REGENIE_BUNDLE_URI="${DSUB_BIN_URI}/regenie_bundle"
export DSUB_LOG_URI="${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/logs/dsub"

if [[ "${PREFLIGHT_ONLY}" -eq 0 ]]; then
    if [[ -n "${REGENIE:-}" && -x "${REGENIE}" ]]; then
        export REGENIE
    elif command -v regenie >/dev/null 2>&1; then
        export REGENIE="$(command -v regenie)"
    elif [[ -x /opt/workbench-tools/binaries/bin/regenie ]]; then
        export REGENIE="/opt/workbench-tools/binaries/bin/regenie"
    else
        echo "ERROR: could not find regenie; set REGENIE to an executable path." >&2
        exit 1
    fi
else
    export REGENIE="${REGENIE:-}"
fi

LOCAL_SOURCE_DIR="${SCRIPT_DIR}/data/regenie/ses_ea_proxy_scrap/${FACTOR18_GWAS_SOURCE_NAME}"
WORKSPACE_SOURCE_DIR="${DX_REGENIE_INPUT_DIR}/${FACTOR18_GWAS_SOURCE_NAME}"
SOURCE_DIR=""
SOURCE_SCORE_TABLE=""
SOURCE_COVAR=""
SOURCE_KEEP=""
INPUT_DIR="${DX_REGENIE_INPUT_DIR}/${FACTOR18_GWAS_INPUT_NAME}"
OUTPUT_DIR="${DX_REGENIE_OUTPUT_DIR}/${FACTOR18_GWAS_OUTPUT_NAME}"
PHEN="${INPUT_DIR}/phen.txt"
COVAR="${INPUT_DIR}/covar.txt"
KEEP="${INPUT_DIR}/training_iids.txt"
PARAMS="${INPUT_DIR}/${FACTOR18_GWAS_INPUT_NAME}.params.tsv"
SUMMARY="${INPUT_DIR}/${FACTOR18_GWAS_INPUT_NAME}.summary.tsv"
DIAG_DIR="${INPUT_DIR}/diagnostics"

mkdir -p "${SCRIPT_DIR}/logs" "${LOCAL_REGENIE_DIR}" "${DX_REGENIE_INPUT_DIR}" "${DX_REGENIE_OUTPUT_DIR}" "${DX_LOGS_DIR}"

LOG_FILE="${SCRIPT_DIR}/logs/run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

resolve_source_dir() {
    local candidate score covar keep
    for candidate in "${LOCAL_SOURCE_DIR}" "${WORKSPACE_SOURCE_DIR}"; do
        score="${candidate}/gradcpt_flanker_direct_xgb_proxy_scores_wide.tsv"
        covar="${candidate}/covar.txt"
        keep="${candidate}/training_iids.txt"
        if [[ -s "${score}" && -s "${covar}" && -s "${keep}" ]]; then
            SOURCE_DIR="${candidate}"
            SOURCE_SCORE_TABLE="${score}"
            SOURCE_COVAR="${covar}"
            SOURCE_KEEP="${keep}"
            return 0
        fi
    done
    echo "ERROR: could not find required direct-XGB source files." >&2
    echo "  Checked local source:     ${LOCAL_SOURCE_DIR}" >&2
    echo "  Checked workspace source: ${WORKSPACE_SOURCE_DIR}" >&2
    echo "  Run bash run_gradcpt_flanker_direct_xgb_proxy.sh --stage-aggregate first." >&2
    exit 1
}

prepare_input_package() {
    for f in "${SOURCE_SCORE_TABLE}" "${SOURCE_COVAR}" "${SOURCE_KEEP}"; do
        if [[ ! -s "${f}" ]]; then
            echo "ERROR: missing required source file ${f}" >&2
            exit 1
        fi
    done

    if [[ -e "${OUTPUT_DIR}/regenie_gwas.summary.tsv" && "${FORCE}" -eq 0 ]]; then
        echo "ERROR: output directory already has completed metadata: ${OUTPUT_DIR}" >&2
        echo "       Use --force only if you intentionally want the generic runner to reuse/skip existing outputs." >&2
        exit 1
    fi

    mkdir -p "${INPUT_DIR}" "${DIAG_DIR}" "${OUTPUT_DIR}"

    python3 - \
        "${SOURCE_SCORE_TABLE}" \
        "${SOURCE_COVAR}" \
        "${SOURCE_KEEP}" \
        "${INPUT_DIR}" \
        "${REGENIE_PHENO_COL}" \
        "${REGENIE_COVAR_COLS}" \
        "${FORCE}" <<'PY'
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew

score_path = Path(sys.argv[1])
covar_path = Path(sys.argv[2])
keep_path = Path(sys.argv[3])
input_dir = Path(sys.argv[4])
pheno_col = sys.argv[5]
covar_cols = sys.argv[6].split(",")
force = bool(int(sys.argv[7]))
diag_dir = input_dir / "diagnostics"
diag_dir.mkdir(parents=True, exist_ok=True)

features = [
    "ses_ea_proxy_z",
    "gradcpt_flanker_finetuned_ea_proxy_z",
    "gradcpt_flanker_direct_xgb_proxy_z",
]
target_col = "gradcpt_flanker_factor_z"
needed_score = [
    "FID", "IID", "role", "fold_id", "final_model_train_allowed", "teacher_z",
    "gradcpt_flanker_mean_z", target_col, *features,
]
df = pd.read_csv(score_path, sep="\t", usecols=needed_score, dtype={"FID": str, "IID": str, "role": str, "fold_id": str})
for col in ["final_model_train_allowed", "teacher_z", "gradcpt_flanker_mean_z", target_col, *features]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

covar = pd.read_csv(covar_path, sep="\t", dtype={"FID": str, "IID": str})
missing_covars = [c for c in ["FID", "IID", *covar_cols] if c not in covar.columns]
if missing_covars:
    raise SystemExit(f"covar source missing required columns: {missing_covars}")
covar_out = covar[["FID", "IID", *covar_cols]].copy()

keep = pd.read_csv(keep_path, sep="\t", header=None, names=["FID", "IID"], dtype=str)
if len(df) != len(covar_out) or len(df) != len(keep):
    raise SystemExit(f"row-count mismatch: score={len(df)} covar={len(covar_out)} keep={len(keep)}")
if not df[["FID", "IID"]].equals(covar_out[["FID", "IID"]]):
    raise SystemExit("score table and covar rows are not in the same FID/IID order")
if not df[["FID", "IID"]].equals(keep[["FID", "IID"]]):
    raise SystemExit("score table and training_iids rows are not in the same FID/IID order")

feature_mat = df[features].to_numpy(float)
target = df[target_col].to_numpy(float)
features_finite = np.all(np.isfinite(feature_mat), axis=1)
labels_finite = features_finite & np.isfinite(target)
role = df["role"].astype(str)
fold = df["fold_id"].astype(str)
oof = role.eq("oof").to_numpy(bool)
applied = role.eq("applied").to_numpy(bool)
final_allowed = df["final_model_train_allowed"].astype(bool).to_numpy(bool)

pred = np.full(len(df), np.nan, dtype=float)
coef_rows = []

def fit_predict(train_mask: np.ndarray, pred_mask: np.ndarray, fit_name: str, predict_group: str) -> None:
    train_n = int(train_mask.sum())
    pred_n = int(pred_mask.sum())
    if train_n < 1000:
        raise SystemExit(f"too few labels for {fit_name}: {train_n}")
    design_train = np.column_stack([np.ones(train_n), feature_mat[train_mask]])
    coef, *_ = np.linalg.lstsq(design_train, target[train_mask], rcond=None)
    if pred_n:
        design_pred = np.column_stack([np.ones(pred_n), feature_mat[pred_mask]])
        pred[pred_mask] = design_pred.dot(coef)
    row = {
        "fit": fit_name,
        "predict_group": predict_group,
        "train_n": train_n,
        "predict_n": pred_n,
        "intercept": float(coef[0]),
    }
    for col, value in zip(features, coef[1:]):
        row[f"coef_{col}"] = float(value)
    coef_rows.append(row)

for k in range(5):
    train_mask = oof & ~fold.eq(str(k)).to_numpy(bool) & labels_finite
    pred_mask = oof & fold.eq(str(k)).to_numpy(bool) & features_finite
    fit_predict(train_mask, pred_mask, f"oof_fold_{k}_train_other_folds", f"oof_fold_{k}")

train_mask = oof & final_allowed & labels_finite
pred_mask = applied & features_finite
fit_predict(train_mask, pred_mask, "applied_model_train_kinholdout_oof", "applied")

if not np.all(np.isfinite(pred[features_finite])):
    raise SystemExit("some finite-feature rows did not receive calibrated predictions")
raw_mean = float(np.mean(pred[features_finite]))
raw_sd = float(np.std(pred[features_finite], ddof=1))
if not math.isfinite(raw_sd) or raw_sd <= 0:
    raise SystemExit("invalid calibrated raw prediction SD")
df[pheno_col] = (pred - raw_mean) / raw_sd
if not np.all(np.isfinite(df[pheno_col].to_numpy(float))):
    raise SystemExit("phenotype is not finite for all rows")

def corr_row(group: str, score: np.ndarray, target_name: str, values: np.ndarray) -> dict[str, object]:
    finite = np.isfinite(score) & np.isfinite(values)
    x = pd.Series(score[finite])
    y = pd.Series(values[finite])
    return {
        "group": group,
        "target": target_name,
        "n": int(finite.sum()),
        "pearson_r": float(x.corr(y, method="pearson")) if int(finite.sum()) > 1 else math.nan,
        "spearman_r": float(x.corr(y, method="spearman")) if int(finite.sum()) > 1 else math.nan,
    }

score = df[pheno_col].to_numpy(float)
both = np.isfinite(df["gradcpt_flanker_mean_z"].to_numpy(float))
either = np.isfinite(df[target_col].to_numpy(float))
corrs = [
    corr_row("full_cohort", score, "teacher_z", df["teacher_z"].to_numpy(float)),
    corr_row("oof", score[oof], "teacher_z", df.loc[oof, "teacher_z"].to_numpy(float)),
    corr_row("applied", score[applied], "teacher_z", df.loc[applied, "teacher_z"].to_numpy(float)),
    corr_row("either_gradcpt_or_flanker", score[either], target_col, df.loc[either, target_col].to_numpy(float)),
    corr_row("both_gradcpt_and_flanker", score[both], target_col, df.loc[both, target_col].to_numpy(float)),
    corr_row("both_gradcpt_and_flanker", score[both], "gradcpt_flanker_mean_z", df.loc[both, "gradcpt_flanker_mean_z"].to_numpy(float)),
]

dist_rows = []
for mask_name, mask in {
    "calibration_target_18k": either,
    "full_cohort": np.ones(len(df), dtype=bool),
}.items():
    for col in [*features, target_col, pheno_col]:
        vals = df.loc[mask, col].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        dist_rows.append({
            "group": mask_name,
            "variable": col,
            "n": int(len(vals)),
            "mean": float(np.mean(vals)) if len(vals) else math.nan,
            "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else math.nan,
            "skew": float(skew(vals, bias=False)) if len(vals) > 2 else math.nan,
            "min": float(np.min(vals)) if len(vals) else math.nan,
            "p50": float(np.quantile(vals, 0.50)) if len(vals) else math.nan,
            "max": float(np.max(vals)) if len(vals) else math.nan,
        })

pd.DataFrame(coef_rows).assign(
    calibrated_raw_mean_full_cohort=raw_mean,
    calibrated_raw_sd_full_cohort=raw_sd,
).to_csv(diag_dir / "factor18_no_teacher_calibration_coefficients.tsv", sep="\t", index=False)
pd.DataFrame(corrs).to_csv(diag_dir / "factor18_no_teacher_calibration_correlations.tsv", sep="\t", index=False)
pd.DataFrame(dist_rows).to_csv(diag_dir / "factor18_no_teacher_calibration_distributions.tsv", sep="\t", index=False)

phen_out = df[["FID", "IID", pheno_col]].copy()
phen_out.to_csv(input_dir / "phen.txt", sep="\t", index=False)
covar_out.to_csv(input_dir / "covar.txt", sep="\t", index=False)
keep.to_csv(input_dir / "training_iids.txt", sep="\t", index=False, header=False)

params = pd.DataFrame([
    ("input_name", input_dir.name),
    ("phenotype_column", pheno_col),
    ("source_score_table", str(score_path)),
    ("calibration_target", target_col),
    ("calibration_predictors", ",".join(features)),
    ("calibration_type", "fold_safe_5_oof_plus_applied_kinholdout"),
    ("teacher_z_in_calibration", "0"),
    ("calibration_target_labels", int(labels_finite.sum())),
    ("full_cohort_samples", len(df)),
    ("calibrated_raw_mean_full_cohort", raw_mean),
    ("calibrated_raw_sd_full_cohort", raw_sd),
    ("covariate_columns", ",".join(covar_cols)),
])
params.to_csv(input_dir / f"{input_dir.name}.params.tsv", sep="\t", index=False, header=["parameter", "value"])

summary = pd.DataFrame([
    ("gwas_samples", len(df)),
    ("phen_rows", len(phen_out)),
    ("covar_rows", len(covar_out)),
    ("training_iids_rows", len(keep)),
    ("pheno_col", pheno_col),
    ("covar_cols", ",".join(covar_cols)),
    ("n_pcs", 10),
    ("calibration_target_labels", int(labels_finite.sum())),
    ("calibration_target_oof_labels", int((labels_finite & oof).sum())),
    ("calibration_target_applied_labels", int((labels_finite & applied).sum())),
    ("final_model_train_allowed_target_labels", int((labels_finite & oof & final_allowed).sum())),
])
for row in corrs:
    key = f"{row['group']}_{row['target']}_pearson"
    summary.loc[len(summary)] = (key, row["pearson_r"])
summary.to_csv(input_dir / f"{input_dir.name}.summary.tsv", sep="\t", index=False, header=["metric", "value"])

print(f"wrote_input_dir\t{input_dir}")
print(f"phenotype_column\t{pheno_col}")
print(f"gwas_samples\t{len(df)}")
print(f"calibration_target_labels\t{int(labels_finite.sum())}")
for row in corrs:
    print(f"correlation\t{row['group']}\t{row['target']}\t{row['n']}\tpearson={row['pearson_r']:.6f}\tspearman={row['spearman_r']:.6f}")
PY
}

validate_input_package() {
    for f in "${PHEN}" "${COVAR}" "${KEEP}" "${PARAMS}" "${SUMMARY}"; do
        if [[ ! -s "${f}" ]]; then
            echo "ERROR: missing generated input file ${f}" >&2
            exit 1
        fi
    done
    python3 - "$PHEN" "$COVAR" "$KEEP" "$REGENIE_PHENO_COL" "$REGENIE_COVAR_COLS" <<'PY'
import csv
import math
import sys

phen_path, covar_path, keep_path, pheno_col, covar_cols = sys.argv[1:6]
covar_cols = covar_cols.split(",")
with open(phen_path, newline="") as handle:
    phen_header = next(csv.reader(handle, delimiter="\t"))
if phen_header != ["FID", "IID", pheno_col]:
    raise SystemExit(f"phen.txt header mismatch: {phen_header}")
with open(covar_path, newline="") as handle:
    covar_header = next(csv.reader(handle, delimiter="\t"))
expected_covar = ["FID", "IID", *covar_cols]
if covar_header != expected_covar:
    raise SystemExit(f"covar.txt header mismatch: {covar_header} != {expected_covar}")
if "yob_c" in covar_header or "yob_c_sex_c_inter" in covar_header:
    raise SystemExit("age/year-of-birth covariate unexpectedly present")

phen_rows = 0
nonfinite = 0
with open(phen_path, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        phen_rows += 1
        try:
            value = float(row[pheno_col])
        except ValueError:
            nonfinite += 1
            continue
        if not math.isfinite(value):
            nonfinite += 1
covar_rows = sum(1 for _ in open(covar_path)) - 1
keep_rows = sum(1 for _ in open(keep_path))
if phen_rows != covar_rows or phen_rows != keep_rows:
    raise SystemExit(f"row-count mismatch: phen={phen_rows} covar={covar_rows} keep={keep_rows}")
if nonfinite:
    raise SystemExit(f"phenotype contains non-finite rows: {nonfinite}")
print(f"validated_samples\t{phen_rows}")
print(f"phenotype_column\t{pheno_col}")
print(f"covariate_columns\t{','.join(covar_cols)}")
PY
}

echo "=================================================================="
echo "AoU factor18 no-teacher GradCPT/Flanker proxy GWAS — $(date)"
echo "=================================================================="
echo "  WORKSPACE_BUCKET_URI = ${WORKSPACE_BUCKET_URI}"
resolve_source_dir
echo "  Source score table   = ${SOURCE_SCORE_TABLE}"
echo "  Input name           = ${FACTOR18_GWAS_INPUT_NAME}"
echo "  Output name          = ${FACTOR18_GWAS_OUTPUT_NAME}"
echo "  Phenotype column     = ${REGENIE_PHENO_COL}"
echo "  Covariates           = ${REGENIE_COVAR_COLS}"
echo "  Chromosomes          = ${REGENIE_CHROMS}"
echo "  Apply RINT           = ${REGENIE_APPLY_RINT}"
echo "  Smoke                = ${SMOKE}"
echo "  Preflight only       = ${PREFLIGHT_ONLY}"
echo "  Force                = ${FORCE}"
echo "  Log file             = ${LOG_FILE}"
echo ""

prepare_input_package
validate_input_package

echo ""
echo "  Lightweight output path will be:"
echo "  ${OUTPUT_DIR}/lightweight/"
echo ""

if [[ "${PREFLIGHT_ONLY}" -eq 1 ]]; then
    echo "Preflight complete; not submitting REGENIE."
    exit 0
fi

bash "${SCRIPT_DIR}/run_continuous_regenie_gwas.sh" \
    "${FACTOR18_GWAS_INPUT_NAME}" \
    "${FACTOR18_GWAS_OUTPUT_NAME}" \
    $( [[ "${REGENIE_APPLY_RINT}" == "1" ]] && printf '%s' "--apply-rint" || printf '%s' "--no-apply-rint" ) \
    --chroms "${REGENIE_CHROMS}" \
    --step1-block-size "${REGENIE_STEP1_BLOCK_SIZE}" \
    --step2-block-size "${REGENIE_STEP2_BLOCK_SIZE}"

echo ""
echo "=== factor18 no-teacher GradCPT/Flanker GWAS command complete ==="
