#!/usr/bin/env bash
# run_pan_aou_gwas.sh - Pan-UKB-style All of Us HapMap3 residualize-first GWAS.
#
# Runs the whole phenotype-wide GWAS from inside an AoU Verily Jupyter terminal:
#   1. resolve the pre-built HapMap3 HQ bfile and the sample-QC / PC / sex inputs;
#   2. build the unrelated-European GWAS keep-list;
#   3. extract survey responses (ds_survey [+ MHWB off-cycle for BHP/EHW]),
#      physical measurements, and ZIP3 SES context to local CSVs via `bq`;
#   4. hand off to pan_aou_gwas.py, which builds residualized phenotypes and runs
#      the covariate-free PLINK2 linear GWAS per phenotype.
#
# See SPECSHEET.md for the method. This script is idempotent: existing extracts
# and GWAS outputs are reused unless --force is given.
#
# Usage:
#   bash run_pan_aou_gwas.sh --setup-only   # extract + build phenotypes, no GWAS
#   bash run_pan_aou_gwas.sh --smoke        # a few phenotypes end-to-end
#   bash run_pan_aou_gwas.sh --gwas-only    # resubmit missing GWAS batches
#   bash run_pan_aou_gwas.sh                 # full run
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Shared AoU environment (WORKSPACE_CDR, workspace bucket, DX_OUTPUT_DIR, ...).
# aou_downstream_env.sh lives one level up in the repo root.
# shellcheck source=/dev/null
source "${REPO_DIR}/aou_downstream_env.sh"

SETUP_ONLY=0
SMOKE=0
FORCE=0
GWAS_ONLY=0
SMOKE_PHENOS="${SMOKE_PHENOS:-height_cm,bmi_kg_m2,ord_educationlevel_highestgrade,ord_overallhealth_generalhealth,bin_overallhealth_generalhealth__excellent}"
for arg in "$@"; do
  case "${arg}" in
    --setup-only) SETUP_ONLY=1 ;;
    --smoke) SMOKE=1 ;;
    --gwas-only) GWAS_ONLY=1 ;;
    --force) FORCE=1 ;;
    *) echo "Unknown arg: ${arg}" >&2; exit 2 ;;
  esac
done
if [[ "${GWAS_ONLY}" == 1 && "${SETUP_ONLY}" == 1 ]]; then
  echo "ERROR: --gwas-only and --setup-only are mutually exclusive" >&2
  exit 2
fi
if [[ "${SMOKE}" == 1 ]]; then
  # The default smoke phenotypes use ds_survey and physical measurements only.
  # Keep smoke focused and fast unless the caller explicitly includes optional
  # phenotypes and opts into these larger extracts.
  export PAN_AOU_SKIP_FITBIT="${PAN_AOU_SKIP_FITBIT:-1}"
  export PAN_AOU_SKIP_MHWB="${PAN_AOU_SKIP_MHWB:-1}"
  if [[ "${SMOKE_PHENOS}" != *zip3_* ]]; then
    export PAN_AOU_SKIP_ZIP3_SES="${PAN_AOU_SKIP_ZIP3_SES:-1}"
  fi
fi

# --- inputs from the main pipeline ----------------------------------------- #
HM3_BFILE="${HM3_BFILE:-${DX_OUTPUT_DIR}/hapmap3_bfile_hq/hapmap3_bfile_hq}"
HM3_BFILE_URI="${HM3_BFILE_URI:-${DX_OUTPUT_URI}/hapmap3_bfile_hq/hapmap3_bfile_hq}"
EUR_IIDS="${DX_OUTPUT_DIR}/europeans/classified_european_iids.txt"
FIT_PCA_IIDS="${DX_OUTPUT_DIR}/pca_eur/fit_pca_iids.txt"
PROJECTED_PCS="${DX_OUTPUT_DIR}/pca_eur/aou_projected.sscore"
SEX_COVAR="${DX_OUTPUT_DIR}/genetic_sex/sex_covar.txt"
SEX_PLOIDY_QC="${DX_OUTPUT_DIR}/genetic_sex/sex_ploidy_qc.tsv"
EXCLUDE_IDENTICAL="${DX_OUTPUT_DIR}/sample_qc/exclude_identical_component_size_ge3_iids.txt"

for f in "${HM3_BFILE}.fam" "${FIT_PCA_IIDS}" "${PROJECTED_PCS}" "${SEX_COVAR}" "${SEX_PLOIDY_QC}"; do
  [[ -s "${f}" ]] || { echo "ERROR: required input missing: ${f}" >&2; exit 1; }
done

# --- output / working dirs ------------------------------------------------- #
# Keep final phenotypes/GWAS outputs in the durable workspace bucket. Keep
# intermediate BigQuery CSV extracts on local ignored disk by default; streaming
# multi-GB query output directly into gcsfuse can stall before flushing data.
PAN_AOU_OUTDIR="${PAN_AOU_OUTDIR:-${DX_OUTPUT_DIR}/pan_aou_gwas}"
PAN_AOU_OUT_URI="${PAN_AOU_OUT_URI:-${DX_OUTPUT_URI}/pan_aou_gwas}"
PAN_AOU_WORKDIR="${PAN_AOU_WORKDIR:-${REPO_DIR}/data/pan_aou_gwas_work}"
WORK="${PAN_AOU_WORKDIR}"
EXTRACT_DIR="${WORK}/extract"
KEEP_DIR="${WORK}/sample_qc"
GWAS_WORKDIR="${WORK}/gwas"
mkdir -p "${EXTRACT_DIR}" "${KEEP_DIR}" "${GWAS_WORKDIR}" "${PAN_AOU_OUTDIR}"
PAN_AOU_EXTRACT_URI="${PAN_AOU_EXTRACT_URI:-${DX_OUTPUT_URI}/pan_aou_gwas/work/extract}"
PAN_AOU_SEX_COVAR="${PAN_AOU_SEX_COVAR:-${KEEP_DIR}/pan_aou_sex_covar.txt}"
PAN_AOU_SEX_COVAR_SUMMARY="${PAN_AOU_SEX_COVAR_SUMMARY:-${KEEP_DIR}/pan_aou_sex_covar.summary.tsv}"
PAN_AOU_SEX_COVAR_AUDIT="${PAN_AOU_SEX_COVAR_AUDIT:-${KEEP_DIR}/pan_aou_sex_covar.imputed_rows.tsv}"
PAN_AOU_PERSON_AGE_REFERENCE_DATE="${PAN_AOU_PERSON_AGE_REFERENCE_DATE:-2026-07-01}"
if [[ ! "${PAN_AOU_PERSON_AGE_REFERENCE_DATE}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
  echo "ERROR: PAN_AOU_PERSON_AGE_REFERENCE_DATE must be YYYY-MM-DD; got ${PAN_AOU_PERSON_AGE_REFERENCE_DATE}" >&2
  exit 2
fi

if [[ "${SETUP_ONLY}" == 1 ]]; then
  PAN_AOU_GWAS_BACKEND="${PAN_AOU_GWAS_BACKEND:-none}"
elif [[ "${GWAS_ONLY}" == 1 ]]; then
  PAN_AOU_GWAS_BACKEND="${PAN_AOU_GWAS_BACKEND:-dsub}"
elif [[ "${SMOKE}" == 1 ]]; then
  PAN_AOU_GWAS_BACKEND="${PAN_AOU_GWAS_BACKEND:-local}"
else
  PAN_AOU_GWAS_BACKEND="${PAN_AOU_GWAS_BACKEND:-dsub}"
fi
case "${PAN_AOU_GWAS_BACKEND}" in
  none|local|dsub) ;;
  *) echo "ERROR: PAN_AOU_GWAS_BACKEND must be none, local, or dsub; got ${PAN_AOU_GWAS_BACKEND}" >&2; exit 2 ;;
esac

choose_bq_tmp_dataset() {
  local requested="${1:-}" candidate
  if [[ -n "${requested}" ]]; then
    printf '%s\n' "${requested}"
    return 0
  fi
  for candidate in sbayesrc_tmp high_quality_cohort dataset_test2; do
    if bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${candidate}" >/dev/null 2>&1; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  candidate=$(bq --project_id="${GOOGLE_PROJECT}" ls --max_results=100 2>/dev/null |
    awk 'NR > 2 && $1 !~ /^-/ {print $1; exit}')
  if [[ -n "${candidate}" ]]; then
    printf '%s\n' "${candidate}"
    return 0
  fi
  echo "ERROR: no writable BigQuery dataset found; set PAN_AOU_BQ_TMP_DATASET or SBAYESRC_BQ_TMP_DATASET." >&2
  return 1
}

PAN_AOU_BQ_TMP_DATASET_RESOLVED=""
resolve_pan_aou_bq_tmp_dataset() {
  if [[ -z "${PAN_AOU_BQ_TMP_DATASET_RESOLVED}" ]]; then
    PAN_AOU_BQ_TMP_DATASET_RESOLVED="$(choose_bq_tmp_dataset "${PAN_AOU_BQ_TMP_DATASET:-${SBAYESRC_BQ_TMP_DATASET:-}}")"
    if ! bq --project_id="${GOOGLE_PROJECT}" show "${GOOGLE_PROJECT}:${PAN_AOU_BQ_TMP_DATASET_RESOLVED}" >/dev/null 2>&1; then
      echo "ERROR: BigQuery dataset ${GOOGLE_PROJECT}:${PAN_AOU_BQ_TMP_DATASET_RESOLVED} does not exist or is not readable." >&2
      return 1
    fi
  fi
  printf '%s\n' "${PAN_AOU_BQ_TMP_DATASET_RESOLVED}"
}

bq_query_to_csv() {
  local label="$1" out_csv="$2" sql="$3"
  local tmp_dataset tmp_table tmp_ref sql_file log_file shard_uri shard_dir tmp_csv
  tmp_dataset="$(resolve_pan_aou_bq_tmp_dataset)"
  tmp_table="pan_aou_${label}_$(date +%Y%m%d_%H%M%S)_$$"
  tmp_ref="${GOOGLE_PROJECT}:${tmp_dataset}.${tmp_table}"
  sql_file="${EXTRACT_DIR}/${label}.sql"
  log_file="${EXTRACT_DIR}/${label}.bq.log"
  shard_uri="${PAN_AOU_EXTRACT_URI}/${label}/${label}-*.csv"
  shard_dir="${EXTRACT_DIR}/${label}_shards"
  tmp_csv="${out_csv}.tmp"

  printf '%s\n' "${sql}" > "${sql_file}"
  rm -rf "${shard_dir}"
  mkdir -p "${shard_dir}"
  rm -f "${tmp_csv}"

  echo "  Querying into temporary BigQuery table ${tmp_ref} ..."
  if ! bq --project_id="${GOOGLE_PROJECT}" query \
      --use_legacy_sql=false \
      --replace \
      --destination_table="${tmp_ref}" \
      < "${sql_file}" > "${log_file}" 2>&1; then
    cat "${log_file}" >&2
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_ref}" >/dev/null 2>&1 || true
    return 1
  fi

  echo "  Exporting ${label} shards to ${PAN_AOU_EXTRACT_URI}/${label}/ ..."
  gcloud storage rm -r "${PAN_AOU_EXTRACT_URI}/${label}" \
    --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true
  if ! bq --project_id="${GOOGLE_PROJECT}" extract \
      --destination_format=CSV \
      --print_header=true \
      "${tmp_ref}" \
      "${shard_uri}" >> "${log_file}" 2>&1; then
    cat "${log_file}" >&2
    bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_ref}" >/dev/null 2>&1 || true
    return 1
  fi
  bq --project_id="${GOOGLE_PROJECT}" rm -f -t "${tmp_ref}" >> "${log_file}" 2>&1 || true

  if ! gcloud storage cp "${shard_uri}" "${shard_dir}/" \
      --billing-project="${GOOGLE_PROJECT}" >> "${log_file}" 2>&1; then
    cat "${log_file}" >&2
    return 1
  fi

  local first=1 shard
  for shard in "${shard_dir}"/*.csv; do
    [[ -s "${shard}" ]] || continue
    if [[ "${first}" == 1 ]]; then
      cat "${shard}" > "${tmp_csv}"
      first=0
    else
      tail -n +2 "${shard}" >> "${tmp_csv}"
    fi
  done
  if [[ ! -s "${tmp_csv}" ]]; then
    echo "ERROR: ${label} export produced no local CSV shards." >&2
    return 1
  fi
  mv "${tmp_csv}" "${out_csv}"
}

configure_pan_aou_dsub() {
  export PLINK2="${PLINK2:-/opt/workbench-tools/binaries/bin/plink2}"
  [[ -x "${PLINK2}" ]] || { echo "ERROR: plink2 not found/executable at ${PLINK2}" >&2; return 1; }
  command -v dsub >/dev/null 2>&1 || { echo "ERROR: dsub not found on PATH" >&2; return 1; }

  export DSUB_PROVIDER="${DSUB_PROVIDER:-google-batch}"
  export DSUB_REGION="${DSUB_REGION:-us-central1}"
  export DSUB_NETWORK="${DSUB_NETWORK:-projects/${GOOGLE_PROJECT}/global/networks/network}"
  export DSUB_SUBNETWORK="${DSUB_SUBNETWORK:-projects/${GOOGLE_PROJECT}/regions/${DSUB_REGION}/subnetworks/subnetwork}"
  export DSUB_IMAGE="${DSUB_IMAGE:-marketplace.gcr.io/google/ubuntu2204}"
  DSUB_PET_SA="${DSUB_PET_SA:-$(gcloud config get-value account 2>/dev/null || true)}"
  [[ -n "${DSUB_PET_SA}" ]] || { echo "ERROR: could not determine pet service account via gcloud config get-value account" >&2; return 1; }
  export DSUB_PET_SA

  export DSUB_BIN_URI="${DSUB_BIN_URI:-${WORKSPACE_BUCKET_URI}/bin}"
  export DSUB_PLINK2_GS="${DSUB_PLINK2_GS:-${DSUB_BIN_URI}/plink2}"
  export DSUB_LOG_URI="${DSUB_LOG_URI:-${DX_OUTPUT_URI}/logs/dsub}"
  export PAN_AOU_GWAS_DSUB_MIN_CORES="${PAN_AOU_GWAS_DSUB_MIN_CORES:-16}"
  export PAN_AOU_GWAS_DSUB_MIN_RAM="${PAN_AOU_GWAS_DSUB_MIN_RAM:-64}"
  export PAN_AOU_GWAS_DSUB_DISK_SIZE="${PAN_AOU_GWAS_DSUB_DISK_SIZE:-220}"
  export PAN_AOU_GWAS_DSUB_DISK_TYPE="${PAN_AOU_GWAS_DSUB_DISK_TYPE:-pd-ssd}"
  export DSUB_BOOT_DISK_SIZE="${DSUB_BOOT_DISK_SIZE:-50}"
}

run_pan_aou_gwas_dsub() {
  local manifest="${PAN_AOU_GWAS_MANIFEST:-${PAN_AOU_OUTDIR}/metadata/phenotype_manifest.tsv}"
  local batch_size="${PAN_AOU_GWAS_BATCH_SIZE:-64}"
  local batch_stage_base="${PAN_AOU_GWAS_BATCH_URI:-${PAN_AOU_OUT_URI}/work/gwas_batches}"
  local batch_stage_uri="${batch_stage_base}/$(date +%Y%m%d_%H%M%S)_$$"
  local plan="${GWAS_WORKDIR}/batch_plan.tsv"
  local summary="${GWAS_WORKDIR}/batch_plan.summary.tsv"
  local pending batches tasks_tsv dsub_out dsub_rc dsub_job_id terminal_count

  [[ -s "${manifest}" ]] || { echo "ERROR: phenotype manifest missing: ${manifest}" >&2; return 1; }
  configure_pan_aou_dsub

  echo "Preparing pending GWAS batches for dsub ..."
  local prep_args=(
    "${SCRIPT_DIR}/scripts/prepare_gwas_batches.py"
    --manifest "${manifest}" \
    --workdir "${GWAS_WORKDIR}" \
    --batch-size "${batch_size}"
  )
  [[ "${FORCE}" == 1 ]] && prep_args+=(--force)
  python3 "${prep_args[@]}"
  pending="$(awk -F'\t' '$1 == "pending_phenotypes" {print $2}' "${summary}")"
  batches="$(awk -F'\t' '$1 == "batches" {print $2}' "${summary}")"
  if [[ "${pending}" == "0" ]]; then
    echo "  All pan-AoU GWAS outputs already exist; no dsub tasks to submit."
    return 0
  fi

  echo "  pending phenotypes=${pending}; batches=${batches}; batch_size=${batch_size}"
  echo "Staging PLINK2 and batch phenotype files ..."
  gcloud storage cp "${PLINK2}" "${DSUB_PLINK2_GS}" --billing-project="${GOOGLE_PROJECT}" >/dev/null
  gcloud storage rm -r "${batch_stage_uri}" --billing-project="${GOOGLE_PROJECT}" >/dev/null 2>&1 || true

  tail -n +2 "${plan}" | while IFS=$'\t' read -r batch_index pheno_tsv keep_tsv manifest_tsv n_phenotypes; do
    [[ -n "${batch_index}" ]] || continue
    gcloud storage cp "${pheno_tsv}" "${keep_tsv}" "${manifest_tsv}" "${batch_stage_uri}/" \
      --billing-project="${GOOGLE_PROJECT}" >/dev/null
  done

  tasks_tsv="${REPO_DIR}/logs/dsub_pan_aou_gwas_$(date +%Y%m%d_%H%M%S).tsv"
  mkdir -p "$(dirname "${tasks_tsv}")"
  {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      '--env BATCH_INDEX' \
      '--input PLINK2' '--input BED' '--input BIM' '--input FAM' \
      '--input PHENO' '--input KEEP' '--input BATCH_MANIFEST' \
      '--output-recursive OUTDIR'
    tail -n +2 "${plan}" | while IFS=$'\t' read -r batch_index pheno_tsv keep_tsv manifest_tsv n_phenotypes; do
      [[ -n "${batch_index}" ]] || continue
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${batch_index}" \
        "${DSUB_PLINK2_GS}" \
        "${HM3_BFILE_URI}.bed" \
        "${HM3_BFILE_URI}.bim" \
        "${HM3_BFILE_URI}.fam" \
        "${batch_stage_uri}/$(basename "${pheno_tsv}")" \
        "${batch_stage_uri}/$(basename "${keep_tsv}")" \
        "${batch_stage_uri}/$(basename "${manifest_tsv}")" \
        "${PAN_AOU_OUT_URI}/gwas/"
    done
  } > "${tasks_tsv}"

  echo "Submitting pan-AoU GWAS dsub array: ${batches} tasks"
  echo "  provider=${DSUB_PROVIDER} region=${DSUB_REGION} cores=${PAN_AOU_GWAS_DSUB_MIN_CORES} ram=${PAN_AOU_GWAS_DSUB_MIN_RAM}G disk=${PAN_AOU_GWAS_DSUB_DISK_SIZE}G"
  echo "  logs=${DSUB_LOG_URI}"
  dsub_out="${tasks_tsv%.tsv}.dsub.out"
  set +e
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
    --name "pan-aou-gwas" \
    --image "${DSUB_IMAGE}" \
    --script "${SCRIPT_DIR}/dsub_pan_aou_gwas_worker.sh" \
    --tasks "${tasks_tsv}" \
    --min-cores "${PAN_AOU_GWAS_DSUB_MIN_CORES}" \
    --min-ram "${PAN_AOU_GWAS_DSUB_MIN_RAM}" \
    --boot-disk-size "${DSUB_BOOT_DISK_SIZE}" \
    --disk-size "${PAN_AOU_GWAS_DSUB_DISK_SIZE}" \
    --disk-type "${PAN_AOU_GWAS_DSUB_DISK_TYPE}" \
    --wait \
    --summary 2>&1 | tee "${dsub_out}"
  dsub_rc=${PIPESTATUS[0]}
  set -e

  dsub_job_id="$(awk '/^Launched job-id:/ {print $NF; exit}' "${dsub_out}")"
  if [[ -n "${dsub_job_id}" ]]; then
    echo "Polling dstat for job ${dsub_job_id} until all ${batches} tasks are terminal ..."
    while true; do
      terminal_count=$(dstat --provider "${DSUB_PROVIDER}" \
        --project "${GOOGLE_PROJECT}" \
        --location "${DSUB_REGION}" \
        --jobs "${dsub_job_id}" \
        --users jupyter \
        --status '*' 2>/dev/null |
        awk 'NR>2 && /SUCCESS|FAILURE|CANCEL/ {c++} END {print c+0}')
      if (( terminal_count >= batches )); then
        echo "  ${terminal_count}/${batches} tasks terminal; verifying outputs"
        break
      fi
      echo "  $(date -u +%H:%M:%SZ) ${terminal_count}/${batches} terminal; waiting 60s ..."
      sleep 60
    done
  fi

  if [[ "${dsub_rc}" -ne 0 ]]; then
    echo "ERROR: pan-AoU GWAS dsub returned ${dsub_rc}; checking for completed shards before failing." >&2
  fi

  python3 "${SCRIPT_DIR}/scripts/prepare_gwas_batches.py" \
    --manifest "${manifest}" \
    --workdir "${GWAS_WORKDIR}" \
    --batch-size "${batch_size}"
  pending="$(awk -F'\t' '$1 == "pending_phenotypes" {print $2}' "${summary}")"
  if [[ "${pending}" != "0" ]]; then
    echo "ERROR: ${pending} pan-AoU GWAS phenotypes are still missing outputs." >&2
    echo "       Re-run the same command to submit only missing batches." >&2
    return 1
  fi
  echo "Pan-AoU GWAS dsub outputs verified complete."
}

# Re-submit missing/stale GWAS batches from an existing phenotype manifest
# without rebuilding BigQuery extracts or phenotype files.
if [[ "${GWAS_ONLY}" == 1 ]]; then
  if [[ "${PAN_AOU_GWAS_BACKEND}" != "dsub" ]]; then
    echo "ERROR: --gwas-only currently requires PAN_AOU_GWAS_BACKEND=dsub" >&2
    exit 2
  fi
  run_pan_aou_gwas_dsub
  echo "Done."
  exit 0
fi

# --- 1. unrelated-EUR keep-list -------------------------------------------- #
# fit_pca_iids is the third-degree-unrelated European PCA-fit set. Intersect
# with the pan-AoU sex covariate and subtract the identical-component exclusion.
#
# The pan-AoU sex covariate starts with the strict main-pipeline sex_covar and
# adds a small pre-specified set of rows from sex_ploidy_qc.tsv:
#   * assigned sex at birth Male + DRAGEN X0/XO -> male;
#   * skipped/prefer-not-to-answer sex at birth + DRAGEN XX/XY -> DRAGEN sex.
if [[ ! -s "${PAN_AOU_SEX_COVAR}" || "${FORCE}" == 1 ||
      "${SEX_COVAR}" -nt "${PAN_AOU_SEX_COVAR}" ||
      "${SEX_PLOIDY_QC}" -nt "${PAN_AOU_SEX_COVAR}" ]]; then
  echo "Building pan-AoU sex covariate ..."
  python3 "${SCRIPT_DIR}/scripts/build_pan_aou_sex_covar.py" \
    --strict-sex-covar "${SEX_COVAR}" \
    --sex-ploidy-qc "${SEX_PLOIDY_QC}" \
    --out "${PAN_AOU_SEX_COVAR}" \
    --summary "${PAN_AOU_SEX_COVAR_SUMMARY}" \
    --audit-out "${PAN_AOU_SEX_COVAR_AUDIT}"
fi

KEEP="${KEEP_DIR}/unrelated_eur.keep"
if [[ ! -s "${KEEP}" || "${FORCE}" == 1 || "${PAN_AOU_SEX_COVAR}" -nt "${KEEP}" ]]; then
  echo "Building unrelated-EUR keep-list ..."
  python3 - "${FIT_PCA_IIDS}" "${PAN_AOU_SEX_COVAR}" "${EXCLUDE_IDENTICAL}" "${KEEP}" <<'PY'
import sys
fit, sexp, excl, out = sys.argv[1:5]
def col(path, idx=-1):
    s=set()
    try:
        with open(path) as f:
            for line in f:
                p=line.split()
                if p: s.add(p[idx])
    except FileNotFoundError:
        pass
    return s
fit_iids=col(fit)
# sex_covar has a header IID/sex_01; take the IID column
sex=set()
with open(sexp) as f:
    header=f.readline().split()
    iidx=header.index("IID") if "IID" in header else 1
    for line in f:
        p=line.split()
        if len(p)>iidx: sex.add(p[iidx])
excl_iids=col(excl)
keep=sorted((fit_iids & sex) - excl_iids, key=lambda x:(0,int(x)) if x.isdigit() else (1,x))
with open(out,"w") as g:
    for iid in keep:
        g.write(f"0\t{iid}\n")
print(f"  keep={len(keep)}  (fit_pca={len(fit_iids)}, sex={len(sex)}, excluded={len(excl_iids)})")
PY
fi

# Person-level age covariate for derived non-survey phenotypes that do not have
# their own measurement/survey date, currently the male-only DRAGEN X0/XO GWAS.
PERSON_AGE_CSV="${EXTRACT_DIR}/person_age.csv"
if [[ ! -s "${PERSON_AGE_CSV}" || "${FORCE}" == 1 ]]; then
  echo "Extracting person age covariate from ${WORKSPACE_CDR}.person ..."
  person_age_sql="
    SELECT
      CAST(p.person_id AS STRING) AS person_id,
      DATE_DIFF(DATE '${PAN_AOU_PERSON_AGE_REFERENCE_DATE}', DATE(p.birth_datetime), DAY)/365.25 AS age_at_reference_date
    FROM \`${WORKSPACE_CDR}.person\` p
    WHERE p.birth_datetime IS NOT NULL
  "
  bq_query_to_csv "person_age" "${PERSON_AGE_CSV}" "${person_age_sql}"
  echo "  wrote ${PERSON_AGE_CSV} ($(wc -l < "${PERSON_AGE_CSV}") lines; reference_date=${PAN_AOU_PERSON_AGE_REFERENCE_DATE})"
fi

# --- 2. extract survey responses ------------------------------------------- #
if [[ "${SMOKE}" == 1 ]]; then
  SURVEY_CSV="${EXTRACT_DIR}/survey_responses.smoke.csv"
else
  SURVEY_CSV="${EXTRACT_DIR}/survey_responses.csv"
fi
if [[ ! -s "${SURVEY_CSV}" || "${FORCE}" == 1 ]]; then
  echo "Extracting survey responses from ${WORKSPACE_CDR}.ds_survey ..."
  SURVEY_SQL_FILTER=""
  if [[ "${SMOKE}" == 1 ]]; then
    SURVEY_SQL_FILTER="AND s.question_concept_id IN (1585940, 1585711)"
  fi
  survey_sql="
    SELECT
      CAST(s.person_id AS STRING)     AS person_id,
      s.survey                        AS survey,
      s.question_concept_id           AS question_concept_id,
      s.question                      AS question,
      s.answer_concept_id             AS answer_concept_id,
      s.answer                        AS answer,
      s.survey_datetime               AS survey_datetime,
      DATE_DIFF(DATE(s.survey_datetime), DATE(p.birth_datetime), DAY)/365.25 AS age_at_survey
    FROM \`${WORKSPACE_CDR}.ds_survey\` s
    JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
    WHERE s.survey_datetime IS NOT NULL
      ${SURVEY_SQL_FILTER}
  "
  bq_query_to_csv "survey_responses" "${SURVEY_CSV}" "${survey_sql}"
  echo "  wrote ${SURVEY_CSV} ($(wc -l < "${SURVEY_CSV}") lines)"
fi

# Behavioral Health / Emotional Health may live in the off-cycle MHWB CDR.
# Extract it too and let the worker union it (deduped by person/question/datetime).
BHP_CSV="${EXTRACT_DIR}/bhp_ehw_responses.csv"
if [[ "${PAN_AOU_SKIP_MHWB:-0}" != 1 && ( ! -s "${BHP_CSV}" || "${FORCE}" == 1 ) ]]; then
  if bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_MHWB_CDR/./:}.survey_conduct" >/dev/null 2>&1; then
    echo "Extracting BHP/EHW responses from ${WORKSPACE_MHWB_CDR} ..."
    bhp_sql="
      SELECT
        CAST(o.person_id AS STRING)                  AS person_id,
        sc.survey_concept_id                         AS survey,
        o.observation_concept_id                     AS question_concept_id,
        qc.concept_name                              AS question,
        o.value_as_concept_id                        AS answer_concept_id,
        COALESCE(a.concept_name, o.value_source_value, o.value_as_string,
                 CAST(o.value_as_number AS STRING))  AS answer,
        o.observation_datetime                       AS survey_datetime,
        DATE_DIFF(DATE(o.observation_datetime), DATE(p.birth_datetime), DAY)/365.25 AS age_at_survey
      FROM \`${WORKSPACE_MHWB_CDR}.observation\` o
      JOIN \`${WORKSPACE_MHWB_CDR}.survey_conduct\` sc
        ON o.questionnaire_response_id = sc.survey_conduct_id
      JOIN \`${WORKSPACE_CDR}.person\` p
        ON p.person_id = o.person_id
      LEFT JOIN \`${WORKSPACE_MHWB_CDR}.concept\` qc ON o.observation_concept_id = qc.concept_id
      LEFT JOIN \`${WORKSPACE_MHWB_CDR}.concept\` a  ON o.value_as_concept_id   = a.concept_id
      WHERE o.observation_datetime IS NOT NULL
    "
    bq_query_to_csv "bhp_ehw_responses" "${BHP_CSV}" "${bhp_sql}" || echo "  WARN: MHWB extract failed; continuing without it."
  else
    echo "  MHWB CDR ${WORKSPACE_MHWB_CDR} not readable; skipping (set PAN_AOU_SKIP_MHWB=1 to silence)."
    : > "${BHP_CSV}"
  fi
fi
if [[ "${PAN_AOU_SKIP_MHWB:-0}" == 1 ]]; then
  BHP_CSV="/dev/null"
fi

# --- 3. extract physical measurements -------------------------------------- #
MEAS_CSV="${EXTRACT_DIR}/measurements.csv"
if [[ ! -s "${MEAS_CSV}" || "${FORCE}" == 1 ]]; then
  echo "Extracting physical measurements from ${WORKSPACE_CDR}.measurement ..."
  measurements_sql="
    SELECT
      CAST(m.person_id AS STRING) AS person_id,
      m.measurement_source_concept_id AS measurement_concept_id,
      c.concept_name              AS measurement_name,
      m.measurement_datetime      AS measurement_datetime,
      m.value_as_number           AS value_as_number,
      m.unit_source_value         AS unit,
      DATE_DIFF(DATE(m.measurement_datetime), DATE(p.birth_datetime), DAY)/365.25 AS age_at_measurement
    FROM \`${WORKSPACE_CDR}.measurement\` m
    JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
    LEFT JOIN \`${WORKSPACE_CDR}.concept\` c ON m.measurement_source_concept_id = c.concept_id
    WHERE m.measurement_source_concept_id IN (
        903133,  -- Height
        903121,  -- Weight
        903124,  -- BMI
        903118,  -- Systolic BP
        903115,  -- Diastolic BP
        903126   -- Heart rate / pulse
      )
      AND m.value_as_number IS NOT NULL
  "
  bq_query_to_csv "measurements" "${MEAS_CSV}" "${measurements_sql}"
  echo "  wrote ${MEAS_CSV} ($(wc -l < "${MEAS_CSV}") lines)"
fi

# --- 3b. extract ZIP3 socioeconomic context -------------------------------- #
ZIP3_SES_CSV="${EXTRACT_DIR}/zip3_ses.csv"
if [[ "${PAN_AOU_SKIP_ZIP3_SES:-0}" != 1 ]]; then
  if [[ ! -s "${ZIP3_SES_CSV}" || "${FORCE}" == 1 ]]; then
    if bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_CDR/./:}.ds_zip_code_socioeconomic" >/dev/null 2>&1; then
      echo "Extracting ZIP3 socioeconomic context from ${WORKSPACE_CDR}.ds_zip_code_socioeconomic ..."
      zip3_ses_sql="
        WITH ranked AS (
          SELECT
            CAST(z.PERSON_ID AS STRING) AS person_id,
            z.OBSERVATION_DATETIME AS observation_datetime,
            DATE_DIFF(DATE(z.OBSERVATION_DATETIME), DATE(p.birth_datetime), DAY)/365.25 AS age_at_observation,
            z.ZIP3_AS_STRING AS zip3_as_string,
            CAST(z.DEPRIVATION_INDEX AS FLOAT64) AS deprivation_index,
            CAST(z.MEDIAN_INCOME AS FLOAT64) AS median_income,
            CAST(z.FRACTION_POVERTY AS FLOAT64) AS fraction_poverty,
            CAST(z.FRACTION_ASSISTED_INCOME AS FLOAT64) AS fraction_assisted_income,
            CAST(z.FRACTION_NO_HEALTH_INS AS FLOAT64) AS fraction_no_health_ins,
            CAST(z.FRACTION_VACANT_HOUSING AS FLOAT64) AS fraction_vacant_housing,
            CAST(z.FRACTION_HIGH_SCHOOL_EDU AS FLOAT64) AS fraction_high_school_edu,
            CAST(z.ACS AS INT64) AS acs,
            ROW_NUMBER() OVER (
              PARTITION BY z.PERSON_ID
              ORDER BY z.OBSERVATION_DATETIME DESC, z.ACS DESC
            ) AS rn
          FROM \`${WORKSPACE_CDR}.ds_zip_code_socioeconomic\` z
          JOIN \`${WORKSPACE_CDR}.person\` p
            ON p.person_id = z.PERSON_ID
          WHERE z.OBSERVATION_DATETIME IS NOT NULL
        )
        SELECT
          person_id,
          observation_datetime,
          age_at_observation,
          zip3_as_string,
          deprivation_index,
          median_income,
          fraction_poverty,
          fraction_assisted_income,
          fraction_no_health_ins,
          fraction_vacant_housing,
          fraction_high_school_edu,
          acs
        FROM ranked
        WHERE rn = 1
        ORDER BY CAST(person_id AS INT64)
      "
      bq_query_to_csv "zip3_ses" "${ZIP3_SES_CSV}" "${zip3_ses_sql}"
      echo "  wrote ${ZIP3_SES_CSV} ($(wc -l < "${ZIP3_SES_CSV}") lines)"
    else
      echo "  WARN: ${WORKSPACE_CDR}.ds_zip_code_socioeconomic not readable; skipping ZIP3 SES phenotypes."
      : > "${ZIP3_SES_CSV}"
    fi
  fi
else
  ZIP3_SES_CSV="/dev/null"
fi

# --- 3c. extract Fitbit activity + sleep (optional) ------------------------ #
FITBIT_ACT_CSV="${EXTRACT_DIR}/fitbit_activity.csv"
FITBIT_SLEEP_CSV="${EXTRACT_DIR}/fitbit_sleep.csv"
if [[ "${PAN_AOU_SKIP_FITBIT:-0}" != 1 ]]; then
  if [[ ( ! -s "${FITBIT_ACT_CSV}" || "${FORCE}" == 1 ) ]] && \
     bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_CDR/./:}.activity_summary" >/dev/null 2>&1; then
    echo "Extracting Fitbit daily activity ..."
    fitbit_activity_sql="
      SELECT
        CAST(a.person_id AS STRING) AS person_id,
        a.steps                      AS steps,
        a.sedentary_minutes          AS sedentary_minutes,
        (a.fairly_active_minutes + a.very_active_minutes) AS active_minutes,
        DATE_DIFF(a.date, DATE(p.birth_datetime), DAY)/365.25 AS age
      FROM \`${WORKSPACE_CDR}.activity_summary\` a
      JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
      WHERE a.steps IS NOT NULL AND a.steps > 0
    "
    bq_query_to_csv "fitbit_activity" "${FITBIT_ACT_CSV}" "${fitbit_activity_sql}" || echo "  WARN: Fitbit activity extract failed."
  fi
  if [[ ( ! -s "${FITBIT_SLEEP_CSV}" || "${FORCE}" == 1 ) ]] && \
     bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_CDR/./:}.sleep_daily_summary" >/dev/null 2>&1; then
    echo "Extracting Fitbit daily sleep ..."
    fitbit_sleep_sql="
      SELECT
        CAST(s.person_id AS STRING) AS person_id,
        s.minute_asleep              AS minute_asleep,
        SAFE_DIVIDE(s.minute_asleep, NULLIF(s.minute_in_bed, 0)) AS sleep_efficiency,
        DATE_DIFF(s.sleep_date, DATE(p.birth_datetime), DAY)/365.25 AS age
      FROM \`${WORKSPACE_CDR}.sleep_daily_summary\` s
      JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
      WHERE LOWER(CAST(s.is_main_sleep AS STRING)) = 'true'
        AND s.minute_asleep IS NOT NULL
    "
    bq_query_to_csv "fitbit_sleep" "${FITBIT_SLEEP_CSV}" "${fitbit_sleep_sql}" || echo "  WARN: Fitbit sleep extract failed."
  fi
  [[ -s "${FITBIT_ACT_CSV}" ]] || echo "  (no Fitbit activity; set PAN_AOU_SKIP_FITBIT=1 to silence)"
  # Chronotype: per-night main-sleep onset clock hour, from sleep_level start times.
  FITBIT_CHRONO_CSV="${EXTRACT_DIR}/fitbit_chronotype.csv"
  if [[ ( ! -s "${FITBIT_CHRONO_CSV}" || "${FORCE}" == 1 ) ]] && \
     bq --project_id="${GOOGLE_PROJECT}" show "${WORKSPACE_CDR/./:}.sleep_level" >/dev/null 2>&1; then
    echo "Extracting Fitbit sleep onset (chronotype) ..."
    fitbit_chronotype_sql="
      WITH onset AS (
        SELECT person_id, sleep_date, MIN(start_datetime) AS sleep_start
        FROM \`${WORKSPACE_CDR}.sleep_level\`
        WHERE LOWER(CAST(is_main_sleep AS STRING)) = 'true'
          AND start_datetime IS NOT NULL
        GROUP BY person_id, sleep_date
      )
      SELECT
        CAST(o.person_id AS STRING) AS person_id,
        EXTRACT(HOUR FROM o.sleep_start) + EXTRACT(MINUTE FROM o.sleep_start)/60.0 AS onset_hour,
        DATE_DIFF(o.sleep_date, DATE(p.birth_datetime), DAY)/365.25 AS age
      FROM onset o JOIN \`${WORKSPACE_CDR}.person\` p USING (person_id)
    "
    bq_query_to_csv "fitbit_chronotype" "${FITBIT_CHRONO_CSV}" "${fitbit_chronotype_sql}" || echo "  WARN: chronotype extract failed (verify sleep_level schema)."
  fi
fi
: "${FITBIT_CHRONO_CSV:=${EXTRACT_DIR}/fitbit_chronotype.csv}"
if [[ "${PAN_AOU_SKIP_FITBIT:-0}" == 1 ]]; then
  FITBIT_ACT_CSV="/dev/null"
  FITBIT_SLEEP_CSV="/dev/null"
  FITBIT_CHRONO_CSV="/dev/null"
fi

# --- 4. build phenotypes + run GWAS ---------------------------------------- #
PY_ARGS=(
  --bfile "${HM3_BFILE}"
  --keep "${KEEP}"
  --sex "${PAN_AOU_SEX_COVAR}"
  --pcs "${PROJECTED_PCS}"
  --sex-ploidy-qc "${SEX_PLOIDY_QC}"
  --person-age-csv "${PERSON_AGE_CSV}"
  --survey-csv "${SURVEY_CSV}"
  --bhp-csv "${BHP_CSV}"
  --measurements-csv "${MEAS_CSV}"
  --zip3-ses-csv "${ZIP3_SES_CSV}"
  --fitbit-activity-csv "${FITBIT_ACT_CSV}"
  --fitbit-sleep-csv "${FITBIT_SLEEP_CSV}"
  --question-manifest "${SCRIPT_DIR}/metadata/survey_question_manifest.tsv"
  --aou-question-concepts "${REPO_DIR}/data/aou_metadata/aou_ds_survey_question_concepts.tsv"
  --ea-proxy-feature-manifest "${SCRIPT_DIR}/metadata/ea_proxy_feature_sources.tsv"
  --ordinal-manifest "${SCRIPT_DIR}/metadata/ordinal_mapping_manifest.tsv"
  --item-inventory "${SCRIPT_DIR}/metadata/survey_item_inventory.tsv"
  --state-clusters "${SCRIPT_DIR}/metadata/state_clusters.tsv"
  --fitbit-chronotype-csv "${FITBIT_CHRONO_CSV}"
  --pfhh-allowlist "${SCRIPT_DIR}/metadata/pfhh_self_allowlist.tsv"
  --composite-manifest "${SCRIPT_DIR}/metadata/composite_items_manifest.tsv"
  --external-scores "${SCRIPT_DIR}/metadata/external_scores.tsv"
  --sex-specific-items "${SCRIPT_DIR}/metadata/sex_specific_items.tsv"
  --outdir "${PAN_AOU_OUTDIR}"
  --gwas-workdir "${GWAS_WORKDIR}"
  --gwas-batch-size "${PAN_AOU_GWAS_BATCH_SIZE:-64}"
)
# Where the ea_proxy ETM/proxy score files live (registry paths expand these).
export PAN_AOU_ETM_COG_DIR="${PAN_AOU_ETM_COG_DIR:-${REPO_DIR}/data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors}"
export PAN_AOU_SES_EA_DIR="${PAN_AOU_SES_EA_DIR:-${DX_REGENIE_INPUT_DIR}/ses_ea_proxy_v2_kinholdout}"
export PAN_AOU_FINETUNED_DIR="${PAN_AOU_FINETUNED_DIR:-${DX_REGENIE_INPUT_DIR}/gradcpt_flanker_finetuned_ea_proxy_ses_ea_proxy_v2_kinholdout}"
export PAN_AOU_DIRECT_XGB_DIR="${PAN_AOU_DIRECT_XGB_DIR:-${DX_REGENIE_INPUT_DIR}/gradcpt_flanker_direct_xgb_proxy_ses_ea_proxy_v2_kinholdout}"
# Back-compat for older registry paths.
export PAN_AOU_COG_DIR="${PAN_AOU_COG_DIR:-${PAN_AOU_ETM_COG_DIR}}"
# Optional residence-state CSV (person_id,state,age) from ZIP3 geography; else
# the worker uses the survey work-address state.
[[ -n "${PAN_AOU_STATE_CSV:-}" && -s "${PAN_AOU_STATE_CSV}" ]] && PY_ARGS+=(--state-csv "${PAN_AOU_STATE_CSV}")
[[ "${FORCE}" == 1 ]] && PY_ARGS+=(--force)
if [[ "${SETUP_ONLY}" == 1 || "${PAN_AOU_GWAS_BACKEND}" == "none" || "${PAN_AOU_GWAS_BACKEND}" == "dsub" ]]; then
  PY_ARGS+=(--skip-gwas)
fi
[[ "${SMOKE}" == 1 ]] && PY_ARGS+=(--phenotypes "${SMOKE_PHENOS}")

echo "Running pan_aou_gwas.py (GWAS backend: ${PAN_AOU_GWAS_BACKEND}) ..."
python3 "${SCRIPT_DIR}/scripts/pan_aou_gwas.py" "${PY_ARGS[@]}"
if [[ "${PAN_AOU_GWAS_BACKEND}" == "dsub" ]]; then
  run_pan_aou_gwas_dsub
fi
echo "Done."
