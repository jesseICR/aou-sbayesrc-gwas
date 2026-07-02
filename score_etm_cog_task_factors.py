#!/usr/bin/env python3
"""Build ETM cognitive task factor scores for the SES-EA proxy cohort.

This is a phenotype-scoring/diagnostic helper. It does not run GWAS. The sample
universe is the already-built SES-EA proxy phenotype cohort in all_scores.tsv,
including both OOF and applied proxy scores.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import FactorAnalysis, PCA
from sklearn.exceptions import ConvergenceWarning


RNG_SEED = 2026
REDUNDANCY_R = 0.95
WEAK_LOADING = 0.20
WINSOR_LO = 0.005
WINSOR_HI = 0.995
LOGIT_EPS = 0.001
RECOMMENDED_SCORE_NAMES = ("dd_patience", "gradcpt_perf", "flanker_efficiency", "emorecog_perf")


@dataclass(frozen=True)
class IndicatorSpec:
    name: str
    source: str
    transform: str
    priority: int
    role: str = "primary"


@dataclass(frozen=True)
class ScoreConfig:
    score_name: str
    task: str
    score_type: str
    indicators: tuple[IndicatorSpec, ...]
    simple_compare: tuple[str, ...] = ()
    notes: str = ""


DD_PRIMARY = ScoreConfig(
    score_name="dd_patience_factor",
    task="delaydiscounting",
    score_type="fa_primary",
    indicators=(
        IndicatorSpec("dd_two_weeks_lnk_rev", "dd_two_weeks_lnk", "neg_identity", 1),
        IndicatorSpec("dd_one_month_lnk_rev", "dd_one_month_lnk", "neg_identity", 2),
        IndicatorSpec("dd_one_year_lnk_rev", "dd_one_year_lnk", "neg_identity", 3),
        IndicatorSpec("dd_ten_years_lnk_rev", "dd_ten_years_lnk", "neg_identity", 4),
    ),
    simple_compare=("dd_simple_lnk", "dd_simple_score"),
    notes="Delay Discounting factor from four delay-specific reversed log-k fields.",
)

GRADCPT_PRIMARY = ScoreConfig(
    score_name="gradcpt_perf_factor",
    task="gradcpt",
    score_type="fa_primary",
    indicators=(
        IndicatorSpec("gradcpt_dprime", "gradcpt_dprime", "identity", 1),
        IndicatorSpec("gradcpt_cv_rtc_neglog", "gradcpt_cv_rtc", "neg_log", 2),
        IndicatorSpec("gradcpt_median_rtc_neglog", "gradcpt_median_rtc", "neg_log", 3),
    ),
    simple_compare=("gradcpt_simple_dprime", "gradcpt_simple_score"),
    notes="GradCPT primary factor: dprime plus RT consistency/speed, no component accuracy rates.",
)

GRADCPT_COMPONENT = ScoreConfig(
    score_name="gradcpt_component_factor",
    task="gradcpt",
    score_type="fa_sensitivity",
    indicators=(
        IndicatorSpec("gradcpt_go_accuracy_logit", "gradcpt_go_accuracy", "logit", 1),
        IndicatorSpec("gradcpt_nogo_accuracy_logit", "gradcpt_nogo_accuracy", "logit", 2),
        IndicatorSpec("gradcpt_cv_rtc_neglog", "gradcpt_cv_rtc", "neg_log", 3),
        IndicatorSpec("gradcpt_median_rtc_neglog", "gradcpt_median_rtc", "neg_log", 4),
    ),
    simple_compare=("gradcpt_simple_dprime", "gradcpt_simple_score"),
    notes="GradCPT sensitivity factor using go/no-go component accuracies instead of dprime.",
)

FLANKER_PRIMARY = ScoreConfig(
    score_name="flanker_perf_factor",
    task="flanker",
    score_type="fa_primary_candidate",
    indicators=(
        IndicatorSpec("flanker_rcs_incongruent_log", "flanker_rcs_incongruent", "log_plus_eps", 1),
        IndicatorSpec("flanker_rcs_congruent_log", "flanker_rcs_congruent", "log_plus_eps", 2),
        IndicatorSpec("flanker_accuracy_interference_rev", "flanker_accuracy_interference", "neg_identity", 3),
        IndicatorSpec("flanker_median_rt_interference_rev", "flanker_median_rt_interference", "neg_identity", 4),
    ),
    simple_compare=("flanker_simple_score", "flanker_simple_rcs_interference"),
    notes="Flanker one-factor candidate blending efficiency and interference.",
)

FLANKER_EFFICIENCY = ScoreConfig(
    score_name="flanker_efficiency_unit_mean",
    task="flanker",
    score_type="split_unit_mean",
    indicators=(
        IndicatorSpec("flanker_rcs_incongruent_log", "flanker_rcs_incongruent", "log_plus_eps", 1),
        IndicatorSpec("flanker_rcs_congruent_log", "flanker_rcs_congruent", "log_plus_eps", 2),
    ),
    simple_compare=("flanker_simple_score",),
    notes="Predeclared Flanker split score for condition-level speed-accuracy efficiency.",
)

FLANKER_INTERFERENCE = ScoreConfig(
    score_name="flanker_interference_unit_mean",
    task="flanker",
    score_type="split_unit_mean",
    indicators=(
        IndicatorSpec("flanker_accuracy_interference_rev", "flanker_accuracy_interference", "neg_identity", 1),
        IndicatorSpec("flanker_median_rt_interference_rev", "flanker_median_rt_interference", "neg_identity", 2),
    ),
    simple_compare=("flanker_simple_rcs_interference",),
    notes="Predeclared Flanker split score for lower interference cost.",
)

EMORECOG_EFFICIENCY = ScoreConfig(
    score_name="emorecog_efficiency_factor",
    task="emorecog",
    score_type="fa_primary",
    indicators=(
        IndicatorSpec("emorecog_happy_rcs_log", "emorecog_happy_rcs", "log_plus_eps", 1),
        IndicatorSpec("emorecog_angry_rcs_log", "emorecog_angry_rcs", "log_plus_eps", 2),
        IndicatorSpec("emorecog_fearful_rcs_log", "emorecog_fearful_rcs", "log_plus_eps", 3),
        IndicatorSpec("emorecog_sad_rcs_log", "emorecog_sad_rcs", "log_plus_eps", 4),
    ),
    simple_compare=("emorecog_simple_accuracy", "emorecog_simple_score", "emorecog_accuracy_factor"),
    notes="Emotional Recognition per-emotion trial-derived rate-correct efficiency factor.",
)

EMORECOG_SUMMARY_EFFICIENCY = ScoreConfig(
    score_name="emorecog_summary_efficiency_factor",
    task="emorecog",
    score_type="fa_sensitivity",
    indicators=(
        IndicatorSpec("emorecog_happy_summary_eff_log", "emorecog_happy_summary_eff", "log_plus_eps", 1),
        IndicatorSpec("emorecog_angry_summary_eff_log", "emorecog_angry_summary_eff", "log_plus_eps", 2),
        IndicatorSpec("emorecog_fearful_summary_eff_log", "emorecog_fearful_summary_eff", "log_plus_eps", 3),
        IndicatorSpec("emorecog_sad_summary_eff_log", "emorecog_sad_summary_eff", "log_plus_eps", 4),
    ),
    simple_compare=("emorecog_simple_accuracy", "emorecog_simple_score"),
    notes="Emotional Recognition summary-field per-emotion accuracy / median correct RT efficiency factor.",
)

EMORECOG_ACCURACY = ScoreConfig(
    score_name="emorecog_accuracy_factor",
    task="emorecog",
    score_type="fa_sensitivity",
    indicators=(
        IndicatorSpec("emorecog_happy_accuracy_logit", "emorecog_happy_accuracy", "logit", 1),
        IndicatorSpec("emorecog_angry_accuracy_logit", "emorecog_angry_accuracy", "logit", 2),
        IndicatorSpec("emorecog_fearful_accuracy_logit", "emorecog_fearful_accuracy", "logit", 3),
        IndicatorSpec("emorecog_sad_accuracy_logit", "emorecog_sad_accuracy", "logit", 4),
    ),
    simple_compare=("emorecog_simple_accuracy", "emorecog_simple_score"),
    notes="Emotional Recognition accuracy-only per-emotion sensitivity factor.",
)

EMORECOG_EFFICIENCY_UNIT_MEAN = ScoreConfig(
    score_name="emorecog_efficiency_unit_mean",
    task="emorecog",
    score_type="unit_mean_fallback",
    indicators=EMORECOG_EFFICIENCY.indicators,
    simple_compare=("emorecog_simple_accuracy", "emorecog_efficiency_factor"),
    notes="Emotional Recognition fallback unit mean of z-scored per-emotion efficiency indicators.",
)

EMORECOG_SPEED = ScoreConfig(
    score_name="emorecog_speed_diagnostic",
    task="emorecog",
    score_type="diagnostic",
    indicators=(
        IndicatorSpec("emorecog_median_rtc_neglog", "emorecog_median_rtc", "neg_log", 1),
        IndicatorSpec("emorecog_mean_rtc_neglog", "emorecog_mean_rtc", "neg_log", 2),
        IndicatorSpec("emorecog_cv_rtc_neglog", "emorecog_cv_rtc", "neg_log", 3),
    ),
    simple_compare=("emorecog_simple_median_rtc",),
    notes="Emotional Recognition speed-only diagnostic, not a recommended primary construct.",
)

EMORECOG_SCORE_RT = ScoreConfig(
    score_name="emorecog_score_rt_factor",
    task="emorecog",
    score_type="pca_primary",
    indicators=(
        IndicatorSpec("emorecog_score", "emorecog_score", "identity", 1),
        IndicatorSpec("emorecog_cv_rtc_neglog", "emorecog_cv_rtc", "neg_log", 2),
        IndicatorSpec("emorecog_median_rtc_neglog", "emorecog_median_rtc", "neg_log", 3),
    ),
    simple_compare=("emorecog_simple_score", "emorecog_simple_accuracy", "emorecog_simple_median_rtc"),
    notes="Emotional Recognition GradCPT-analog PC1: score plus RT consistency/speed.",
)

SIMPLE_SCORES = (
    ScoreConfig(
        score_name="dd_simple_lnk",
        task="delaydiscounting",
        score_type="official_simple",
        indicators=(IndicatorSpec("dd_lnk_rev", "dd_lnk", "neg_identity", 1),),
        notes="Official/simple Delay Discounting comparison: reversed lnk.",
    ),
    ScoreConfig(
        score_name="dd_simple_score",
        task="delaydiscounting",
        score_type="simple_validation",
        indicators=(IndicatorSpec("dd_score_log1p", "dd_score", "log1p", 1),),
        notes="Participant-facing Delay Discounting halving-time score, log1p-transformed.",
    ),
    ScoreConfig(
        score_name="gradcpt_simple_dprime",
        task="gradcpt",
        score_type="official_simple",
        indicators=(IndicatorSpec("gradcpt_dprime", "gradcpt_dprime", "identity", 1),),
        notes="Official/simple GradCPT comparison: dprime.",
    ),
    ScoreConfig(
        score_name="gradcpt_simple_score",
        task="gradcpt",
        score_type="simple_validation",
        indicators=(IndicatorSpec("gradcpt_score_scaled", "gradcpt_score", "scale_100", 1),),
        notes="Participant-facing GradCPT no-go percent-correct score.",
    ),
    ScoreConfig(
        score_name="flanker_simple_score",
        task="flanker",
        score_type="official_simple",
        indicators=(IndicatorSpec("flanker_score", "flanker_score", "identity", 1),),
        notes="Participant-facing Flanker overall rate-correct score.",
    ),
    ScoreConfig(
        score_name="flanker_simple_rcs_interference",
        task="flanker",
        score_type="simple_validation",
        indicators=(IndicatorSpec("flanker_rcs_interference_rev", "flanker_rcs_interference", "neg_identity", 1),),
        notes="Reversed Flanker RCS interference comparison.",
    ),
    ScoreConfig(
        score_name="flanker_accuracy_interference_single",
        task="flanker",
        score_type="split_component",
        indicators=(IndicatorSpec("flanker_accuracy_interference_rev", "flanker_accuracy_interference", "neg_identity", 1),),
        notes="Single reversed Flanker accuracy-interference component for unstable interference composites.",
    ),
    ScoreConfig(
        score_name="flanker_rt_interference_single",
        task="flanker",
        score_type="split_component",
        indicators=(IndicatorSpec("flanker_median_rt_interference_rev", "flanker_median_rt_interference", "neg_identity", 1),),
        notes="Single reversed Flanker RT-interference component for unstable interference composites.",
    ),
    ScoreConfig(
        score_name="emorecog_simple_accuracy",
        task="emorecog",
        score_type="official_simple",
        indicators=(IndicatorSpec("emorecog_accuracy_logit", "emorecog_accuracy", "logit", 1),),
        notes="Overall Emotional Recognition accuracy, clipped-logit transformed.",
    ),
    ScoreConfig(
        score_name="emorecog_simple_score",
        task="emorecog",
        score_type="simple_validation",
        indicators=(IndicatorSpec("emorecog_score", "emorecog_score", "identity", 1),),
        notes="Participant-facing Emotional Recognition number-correct score.",
    ),
    ScoreConfig(
        score_name="emorecog_simple_median_rtc",
        task="emorecog",
        score_type="simple_validation",
        indicators=(IndicatorSpec("emorecog_median_rtc_neglog", "emorecog_median_rtc", "neg_log", 1),),
        notes="Overall Emotional Recognition median correct RT speed diagnostic.",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--etm-dataset", required=True, help="ETM BigQuery dataset as project.dataset.")
    parser.add_argument("--workspace-cdr", required=True, help="Main CDR BigQuery dataset as project.dataset; used for person.birth_datetime.")
    parser.add_argument("--bq-temp-dataset", required=True, help="Existing writable BigQuery temp dataset name in GOOGLE_PROJECT.")
    parser.add_argument("--ses-ea-dir", type=Path, required=True, help="SES-EA proxy regenie_input directory.")
    parser.add_argument("--work-dir", type=Path, required=True, help="Local on-platform scratch output directory.")
    parser.add_argument("--workspace-scrap-dir", type=Path, required=True, help="Workspace bucket scrap dir for aggregate outputs.")
    parser.add_argument("--reuse-extracts", action="store_true", help="Reuse existing ETM extract CSV; fail if it is missing.")
    parser.add_argument("--force", action="store_true", help="Force re-extraction and overwrite existing output files.")
    parser.add_argument("--stage-aggregate", action="store_true", help="Stage aggregate diagnostics to workspace bucket scrap directory.")
    parser.add_argument("--seed", type=int, default=RNG_SEED)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"ERROR: {message}")


def run_command(cmd: list[str], *, stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        subprocess.run(cmd, check=True)
        return
    with stdout_path.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, check=True, stdout=out)


def sql_ref(dataset: str, table: str) -> str:
    return f"`{dataset}.{table}`"


def cli_ref(dataset: str, table: str) -> str:
    project, ds = dataset.split(".", 1)
    return f"{project}:{ds}.{table}"


def zscore_array(x: np.ndarray) -> tuple[np.ndarray, float, float]:
    x = np.asarray(x, dtype=float)
    mean = float(np.nanmean(x))
    sd = float(np.nanstd(x, ddof=0))
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(x, np.nan, dtype=float), mean, sd
    return (x - mean) / sd, mean, sd


def pearson_spearman(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> tuple[float, float, int]:
    xs = pd.Series(np.asarray(x, dtype=float))
    ys = pd.Series(np.asarray(y, dtype=float))
    ok = np.isfinite(xs) & np.isfinite(ys)
    n = int(ok.sum())
    if n < 3:
        return np.nan, np.nan, n
    return (
        float(stats.pearsonr(xs[ok], ys[ok]).statistic),
        float(stats.spearmanr(xs[ok], ys[ok]).statistic),
        n,
    )


def load_proxy_cohort(ses_ea_dir: Path) -> pd.DataFrame:
    all_scores = ses_ea_dir / "all_scores.tsv"
    base_covar = ses_ea_dir / "base_covar.txt"
    require(all_scores.exists() and all_scores.stat().st_size > 0, f"missing {all_scores}")
    require(base_covar.exists() and base_covar.stat().st_size > 0, f"missing {base_covar}")
    scores = pd.read_csv(all_scores, sep="\t", dtype={"IID": str, "FID": str})
    covar = pd.read_csv(base_covar, sep="\t", dtype={"IID": str, "FID": str}, usecols=["IID", "sex_c"])
    required = {"FID", "IID", "role", "fold_id", "ea_years", "teacher_z", "ses_ea_proxy_z"}
    missing = required - set(scores.columns)
    require(not missing, f"{all_scores} missing columns: {sorted(missing)}")
    cohort = scores.merge(covar, on="IID", how="left", validate="one_to_one")
    cohort["IID"] = cohort["IID"].astype(str)
    for col in ["ea_years", "teacher_z", "ses_ea_proxy_z", "sex_c"]:
        cohort[col] = pd.to_numeric(cohort[col], errors="coerce")
    require(cohort["IID"].is_unique, "all_scores.tsv has duplicate IIDs")
    return cohort


def extract_has_emorecog(path: Path) -> tuple[bool, bool]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split(",")
        has_columns = "emorecog_score" in header and "emorecog_happy_rcs" in header
        if not has_columns or "task" not in header:
            return has_columns, False
        task_idx = header.index("task")
        for line in f:
            fields = line.rstrip("\n").split(",")
            if len(fields) > task_idx and fields[task_idx] == "emorecog":
                return True, True
    return has_columns, False


def existing_extract_matches_args(args: argparse.Namespace) -> bool:
    params_path = args.work_dir / "etm_cog_task_factors_params.tsv"
    if not params_path.exists() or params_path.stat().st_size == 0:
        return False
    params = pd.read_csv(params_path, sep="\t", dtype=str)
    if not {"parameter", "value"}.issubset(params.columns):
        return False
    if "score_name" in params.columns:
        params = params[params["score_name"].isna()]
    observed = dict(zip(params["parameter"], params["value"]))
    return (
        observed.get("etm_dataset") == args.etm_dataset
        and observed.get("workspace_cdr") == args.workspace_cdr
    )


def ensure_bq_extract(args: argparse.Namespace, cohort: pd.DataFrame) -> Path:
    extract_path = args.work_dir / "etm_cog_task_factor_valid_sittings.csv"
    if extract_path.exists() and extract_path.stat().st_size > 0 and not args.force:
        require(
            existing_extract_matches_args(args),
            f"{extract_path} exists but does not match --etm-dataset/--workspace-cdr; rerun with --force",
        )
        has_columns, has_rows = extract_has_emorecog(extract_path)
        require(
            has_columns and has_rows,
            f"{extract_path} predates Emotional Recognition extraction or has no emorecog rows; rerun with --force",
        )
        print(f"Reusing ETM extract: {extract_path}", flush=True)
        return extract_path
    if args.reuse_extracts:
        require(extract_path.exists() and extract_path.stat().st_size > 0, f"--reuse-extracts requested but {extract_path} is missing")
        require(
            existing_extract_matches_args(args),
            f"--reuse-extracts requested but {extract_path} does not match --etm-dataset/--workspace-cdr",
        )
        has_columns, has_rows = extract_has_emorecog(extract_path)
        require(
            has_columns and has_rows,
            f"--reuse-extracts requested but {extract_path} lacks emorecog columns/rows; rerun without --reuse-extracts and with --force",
        )

    google_project = os.environ.get("GOOGLE_PROJECT")
    require(bool(google_project), "GOOGLE_PROJECT is not set")
    require("." in args.etm_dataset, "--etm-dataset must be project.dataset")
    require("." in args.workspace_cdr, "--workspace-cdr must be project.dataset")

    for dataset, table in [
        (args.etm_dataset, "delaydiscounting"),
        (args.etm_dataset, "gradcpt"),
        (args.etm_dataset, "flanker"),
        (args.etm_dataset, "emorecog"),
        (args.workspace_cdr, "person"),
    ]:
        run_command(["bq", "show", cli_ref(dataset, table)])

    args.work_dir.mkdir(parents=True, exist_ok=True)
    temp_table = f"etm_factor_iids_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
    temp_csv = args.work_dir / f"{temp_table}.csv"
    cohort[["IID"]].assign(IID=lambda d: d["IID"].astype(str)).to_csv(temp_csv, index=False)
    temp_ref = f"{google_project}:{args.bq_temp_dataset}.{temp_table}"
    temp_sql = f"`{google_project}.{args.bq_temp_dataset}.{temp_table}`"

    print(f"Loading proxy cohort IIDs to temporary BigQuery table {temp_ref}", flush=True)
    run_command(
        [
            "bq",
            "--project_id",
            google_project,
            "load",
            "--replace",
            "--source_format=CSV",
            "--skip_leading_rows=1",
            temp_ref,
            str(temp_csv),
            "IID:INTEGER",
        ]
    )

    dd = sql_ref(args.etm_dataset, "delaydiscounting")
    grad = sql_ref(args.etm_dataset, "gradcpt")
    flank = sql_ref(args.etm_dataset, "flanker")
    emo = sql_ref(args.etm_dataset, "emorecog")
    person = sql_ref(args.workspace_cdr, "person")
    query = f"""
WITH cohort AS (
  SELECT DISTINCT IID
  FROM {temp_sql}
),
person_birth AS (
  SELECT person_id, birth_datetime
  FROM {person}
),
delaydiscounting AS (
  SELECT
    'delaydiscounting' AS task,
    d.person_id AS IID,
    d.sitting_id,
    d.metadata.test_start_date_time AS test_start_date_time,
    d.metadata.test_end_date_time AS test_end_date_time,
    DATE_DIFF(DATE(d.metadata.test_start_date_time), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_test,
    CAST(d.metadata.response_device AS STRING) AS response_device,
    CAST(d.metadata.touch AS BOOL) AS touch,
    CAST(d.metadata.test_language AS STRING) AS test_language,
    CAST(d.metadata.test_version AS STRING) AS test_version,
    CAST(d.metadata.aou_version AS STRING) AS aou_version,
    CAST(d.metadata.operating_system AS STRING) AS operating_system,
    CAST(d.metadata.user_agent AS STRING) AS user_agent,
    CAST(d.metadata.test_restarted AS BOOL) AS test_restarted,
    CAST(d.outcomes.any_timeouts AS INT64) AS any_timeouts,
    CAST(d.outcomes.flag_median_rt AS INT64) AS dd_flag_median_rt,
    CAST(d.outcomes.flag_catch_trials AS INT64) AS dd_flag_catch_trials,
    CAST(NULL AS INT64) AS gradcpt_flag_trial_flags,
    CAST(NULL AS INT64) AS gradcpt_flag_non_response,
    CAST(NULL AS INT64) AS gradcpt_flag_omission_error_rate,
    CAST(NULL AS INT64) AS flanker_flag_accuracy,
    CAST(NULL AS INT64) AS flanker_flag_trial_flags,
    d.outcomes.score AS dd_score,
    d.outcomes.catch_score AS dd_catch_score,
    d.outcomes.lnk AS dd_lnk,
    d.outcomes.two_weeks_lnk AS dd_two_weeks_lnk,
    d.outcomes.one_month_lnk AS dd_one_month_lnk,
    d.outcomes.one_year_lnk AS dd_one_year_lnk,
    d.outcomes.ten_years_lnk AS dd_ten_years_lnk,
    d.outcomes.mean_rt AS dd_mean_rt,
    d.outcomes.median_rt AS dd_median_rt,
    d.outcomes.sd_rt AS dd_sd_rt,
    CAST(NULL AS FLOAT64) AS gradcpt_dprime,
    CAST(NULL AS FLOAT64) AS gradcpt_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_go_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_nogo_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_score,
    CAST(NULL AS FLOAT64) AS gradcpt_crit,
    CAST(NULL AS FLOAT64) AS gradcpt_mean_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_median_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_sd_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_cv_rtc,
    CAST(NULL AS FLOAT64) AS flanker_score,
    CAST(NULL AS FLOAT64) AS flanker_accuracy,
    CAST(NULL AS FLOAT64) AS flanker_mean_rtc,
    CAST(NULL AS FLOAT64) AS flanker_median_rtc,
    CAST(NULL AS FLOAT64) AS flanker_sd_rtc,
    CAST(NULL AS FLOAT64) AS flanker_rcs_congruent,
    CAST(NULL AS FLOAT64) AS flanker_rcs_incongruent,
    CAST(NULL AS FLOAT64) AS flanker_accuracy_interference,
    CAST(NULL AS FLOAT64) AS flanker_median_rt_interference,
    CAST(NULL AS FLOAT64) AS flanker_rcs_interference,
    CAST(NULL AS INT64) AS emorecog_flag_median_rtc,
    CAST(NULL AS INT64) AS emorecog_flag_same_response,
    CAST(NULL AS INT64) AS emorecog_flag_trial_flags,
    CAST(NULL AS FLOAT64) AS emorecog_score,
    CAST(NULL AS FLOAT64) AS emorecog_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_mean_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_sd_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_cv_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_happy_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_happy_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_happy_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_happy_rcs,
    CAST(NULL AS INT64) AS emorecog_happy_trials,
    CAST(NULL AS INT64) AS emorecog_happy_correct,
    CAST(NULL AS FLOAT64) AS emorecog_happy_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_angry_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_angry_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_angry_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_angry_rcs,
    CAST(NULL AS INT64) AS emorecog_angry_trials,
    CAST(NULL AS INT64) AS emorecog_angry_correct,
    CAST(NULL AS FLOAT64) AS emorecog_angry_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_rcs,
    CAST(NULL AS INT64) AS emorecog_fearful_trials,
    CAST(NULL AS INT64) AS emorecog_fearful_correct,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_sad_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_sad_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_sad_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_sad_rcs,
    CAST(NULL AS INT64) AS emorecog_sad_trials,
    CAST(NULL AS INT64) AS emorecog_sad_correct,
    CAST(NULL AS FLOAT64) AS emorecog_sad_rt_seconds
  FROM {dd} AS d
  JOIN cohort AS c ON c.IID = d.person_id
  LEFT JOIN person_birth AS p ON p.person_id = d.person_id
  WHERE COALESCE(d.outcomes.flag_median_rt, 0) = 0
    AND COALESCE(d.outcomes.flag_catch_trials, 0) = 0
    AND COALESCE(d.metadata.test_restarted, FALSE) = FALSE
),
gradcpt AS (
  SELECT
    'gradcpt' AS task,
    g.person_id AS IID,
    g.sitting_id,
    g.metadata.test_start_date_time AS test_start_date_time,
    g.metadata.test_end_date_time AS test_end_date_time,
    DATE_DIFF(DATE(g.metadata.test_start_date_time), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_test,
    CAST(g.metadata.response_device AS STRING) AS response_device,
    CAST(g.metadata.touch AS BOOL) AS touch,
    CAST(g.metadata.test_language AS STRING) AS test_language,
    CAST(g.metadata.test_version AS STRING) AS test_version,
    CAST(g.metadata.aou_version AS STRING) AS aou_version,
    CAST(g.metadata.operating_system AS STRING) AS operating_system,
    CAST(g.metadata.user_agent AS STRING) AS user_agent,
    CAST(g.metadata.test_restarted AS BOOL) AS test_restarted,
    CAST(NULL AS INT64) AS any_timeouts,
    CAST(NULL AS INT64) AS dd_flag_median_rt,
    CAST(NULL AS INT64) AS dd_flag_catch_trials,
    CAST(g.outcomes.flag_trial_flags AS INT64) AS gradcpt_flag_trial_flags,
    CAST(g.outcomes.flag_non_response AS INT64) AS gradcpt_flag_non_response,
    CAST(g.outcomes.flag_omission_error_rate AS INT64) AS gradcpt_flag_omission_error_rate,
    CAST(NULL AS INT64) AS flanker_flag_accuracy,
    CAST(NULL AS INT64) AS flanker_flag_trial_flags,
    CAST(NULL AS FLOAT64) AS dd_score,
    CAST(NULL AS FLOAT64) AS dd_catch_score,
    CAST(NULL AS FLOAT64) AS dd_lnk,
    CAST(NULL AS FLOAT64) AS dd_two_weeks_lnk,
    CAST(NULL AS FLOAT64) AS dd_one_month_lnk,
    CAST(NULL AS FLOAT64) AS dd_one_year_lnk,
    CAST(NULL AS FLOAT64) AS dd_ten_years_lnk,
    CAST(NULL AS FLOAT64) AS dd_mean_rt,
    CAST(NULL AS FLOAT64) AS dd_median_rt,
    CAST(NULL AS FLOAT64) AS dd_sd_rt,
    g.outcomes.dprime AS gradcpt_dprime,
    g.outcomes.accuracy AS gradcpt_accuracy,
    g.outcomes.go_accuracy AS gradcpt_go_accuracy,
    g.outcomes.nogo_accuracy AS gradcpt_nogo_accuracy,
    g.outcomes.score AS gradcpt_score,
    g.outcomes.crit AS gradcpt_crit,
    g.outcomes.mean_rtc AS gradcpt_mean_rtc,
    g.outcomes.median_rtc AS gradcpt_median_rtc,
    g.outcomes.sd_rtc AS gradcpt_sd_rtc,
    g.outcomes.cv_rtc AS gradcpt_cv_rtc,
    CAST(NULL AS FLOAT64) AS flanker_score,
    CAST(NULL AS FLOAT64) AS flanker_accuracy,
    CAST(NULL AS FLOAT64) AS flanker_mean_rtc,
    CAST(NULL AS FLOAT64) AS flanker_median_rtc,
    CAST(NULL AS FLOAT64) AS flanker_sd_rtc,
    CAST(NULL AS FLOAT64) AS flanker_rcs_congruent,
    CAST(NULL AS FLOAT64) AS flanker_rcs_incongruent,
    CAST(NULL AS FLOAT64) AS flanker_accuracy_interference,
    CAST(NULL AS FLOAT64) AS flanker_median_rt_interference,
    CAST(NULL AS FLOAT64) AS flanker_rcs_interference,
    CAST(NULL AS INT64) AS emorecog_flag_median_rtc,
    CAST(NULL AS INT64) AS emorecog_flag_same_response,
    CAST(NULL AS INT64) AS emorecog_flag_trial_flags,
    CAST(NULL AS FLOAT64) AS emorecog_score,
    CAST(NULL AS FLOAT64) AS emorecog_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_mean_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_sd_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_cv_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_happy_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_happy_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_happy_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_happy_rcs,
    CAST(NULL AS INT64) AS emorecog_happy_trials,
    CAST(NULL AS INT64) AS emorecog_happy_correct,
    CAST(NULL AS FLOAT64) AS emorecog_happy_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_angry_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_angry_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_angry_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_angry_rcs,
    CAST(NULL AS INT64) AS emorecog_angry_trials,
    CAST(NULL AS INT64) AS emorecog_angry_correct,
    CAST(NULL AS FLOAT64) AS emorecog_angry_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_rcs,
    CAST(NULL AS INT64) AS emorecog_fearful_trials,
    CAST(NULL AS INT64) AS emorecog_fearful_correct,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_sad_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_sad_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_sad_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_sad_rcs,
    CAST(NULL AS INT64) AS emorecog_sad_trials,
    CAST(NULL AS INT64) AS emorecog_sad_correct,
    CAST(NULL AS FLOAT64) AS emorecog_sad_rt_seconds
  FROM {grad} AS g
  JOIN cohort AS c ON c.IID = g.person_id
  LEFT JOIN person_birth AS p ON p.person_id = g.person_id
  WHERE COALESCE(g.outcomes.flag_trial_flags, 0) = 0
    AND COALESCE(g.outcomes.flag_non_response, 0) = 0
    AND COALESCE(g.outcomes.flag_omission_error_rate, 0) = 0
    AND COALESCE(g.metadata.test_restarted, FALSE) = FALSE
),
flanker AS (
  SELECT
    'flanker' AS task,
    f.person_id AS IID,
    f.sitting_id,
    f.metadata.test_start_date_time AS test_start_date_time,
    f.metadata.test_end_date_time AS test_end_date_time,
    DATE_DIFF(DATE(f.metadata.test_start_date_time), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_test,
    CAST(f.metadata.response_device AS STRING) AS response_device,
    CAST(f.metadata.touch AS BOOL) AS touch,
    CAST(f.metadata.test_language AS STRING) AS test_language,
    CAST(f.metadata.test_version AS STRING) AS test_version,
    CAST(f.metadata.aou_version AS STRING) AS aou_version,
    CAST(f.metadata.operating_system AS STRING) AS operating_system,
    CAST(f.metadata.user_agent AS STRING) AS user_agent,
    CAST(f.metadata.test_restarted AS BOOL) AS test_restarted,
    CAST(f.outcomes.any_timeouts AS INT64) AS any_timeouts,
    CAST(NULL AS INT64) AS dd_flag_median_rt,
    CAST(NULL AS INT64) AS dd_flag_catch_trials,
    CAST(NULL AS INT64) AS gradcpt_flag_trial_flags,
    CAST(NULL AS INT64) AS gradcpt_flag_non_response,
    CAST(NULL AS INT64) AS gradcpt_flag_omission_error_rate,
    CAST(f.outcomes.flag_accuracy AS INT64) AS flanker_flag_accuracy,
    CAST(f.outcomes.flag_trial_flags AS INT64) AS flanker_flag_trial_flags,
    CAST(NULL AS FLOAT64) AS dd_score,
    CAST(NULL AS FLOAT64) AS dd_catch_score,
    CAST(NULL AS FLOAT64) AS dd_lnk,
    CAST(NULL AS FLOAT64) AS dd_two_weeks_lnk,
    CAST(NULL AS FLOAT64) AS dd_one_month_lnk,
    CAST(NULL AS FLOAT64) AS dd_one_year_lnk,
    CAST(NULL AS FLOAT64) AS dd_ten_years_lnk,
    CAST(NULL AS FLOAT64) AS dd_mean_rt,
    CAST(NULL AS FLOAT64) AS dd_median_rt,
    CAST(NULL AS FLOAT64) AS dd_sd_rt,
    CAST(NULL AS FLOAT64) AS gradcpt_dprime,
    CAST(NULL AS FLOAT64) AS gradcpt_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_go_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_nogo_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_score,
    CAST(NULL AS FLOAT64) AS gradcpt_crit,
    CAST(NULL AS FLOAT64) AS gradcpt_mean_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_median_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_sd_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_cv_rtc,
    f.outcomes.score AS flanker_score,
    f.outcomes.accuracy AS flanker_accuracy,
    f.outcomes.mean_rtc AS flanker_mean_rtc,
    f.outcomes.median_rtc AS flanker_median_rtc,
    f.outcomes.sd_rtc AS flanker_sd_rtc,
    f.outcomes.rcs_congruent AS flanker_rcs_congruent,
    f.outcomes.rcs_incongruent AS flanker_rcs_incongruent,
    f.outcomes.accuracy_interference AS flanker_accuracy_interference,
    f.outcomes.median_rt_interference AS flanker_median_rt_interference,
    f.outcomes.rcs_interference AS flanker_rcs_interference,
    CAST(NULL AS INT64) AS emorecog_flag_median_rtc,
    CAST(NULL AS INT64) AS emorecog_flag_same_response,
    CAST(NULL AS INT64) AS emorecog_flag_trial_flags,
    CAST(NULL AS FLOAT64) AS emorecog_score,
    CAST(NULL AS FLOAT64) AS emorecog_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_mean_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_sd_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_cv_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_happy_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_happy_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_happy_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_happy_rcs,
    CAST(NULL AS INT64) AS emorecog_happy_trials,
    CAST(NULL AS INT64) AS emorecog_happy_correct,
    CAST(NULL AS FLOAT64) AS emorecog_happy_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_angry_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_angry_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_angry_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_angry_rcs,
    CAST(NULL AS INT64) AS emorecog_angry_trials,
    CAST(NULL AS INT64) AS emorecog_angry_correct,
    CAST(NULL AS FLOAT64) AS emorecog_angry_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_rcs,
    CAST(NULL AS INT64) AS emorecog_fearful_trials,
    CAST(NULL AS INT64) AS emorecog_fearful_correct,
    CAST(NULL AS FLOAT64) AS emorecog_fearful_rt_seconds,
    CAST(NULL AS FLOAT64) AS emorecog_sad_accuracy,
    CAST(NULL AS FLOAT64) AS emorecog_sad_median_rtc,
    CAST(NULL AS FLOAT64) AS emorecog_sad_summary_eff,
    CAST(NULL AS FLOAT64) AS emorecog_sad_rcs,
    CAST(NULL AS INT64) AS emorecog_sad_trials,
    CAST(NULL AS INT64) AS emorecog_sad_correct,
    CAST(NULL AS FLOAT64) AS emorecog_sad_rt_seconds
  FROM {flank} AS f
  JOIN cohort AS c ON c.IID = f.person_id
  LEFT JOIN person_birth AS p ON p.person_id = f.person_id
  WHERE COALESCE(f.outcomes.flag_accuracy, 0) = 0
    AND COALESCE(f.outcomes.flag_trial_flags, 0) = 0
    AND COALESCE(f.metadata.test_restarted, FALSE) = FALSE
),
emorecog_base AS (
  SELECT
    e.person_id,
    e.sitting_id,
    e.metadata.test_start_date_time AS test_start_date_time,
    e.metadata.test_end_date_time AS test_end_date_time,
    DATE_DIFF(DATE(e.metadata.test_start_date_time), DATE(p.birth_datetime), DAY) / 365.25 AS age_at_test,
    CAST(e.metadata.response_device AS STRING) AS response_device,
    CAST(e.metadata.touch AS BOOL) AS touch,
    CAST(e.metadata.test_language AS STRING) AS test_language,
    CAST(e.metadata.test_version AS STRING) AS test_version,
    CAST(e.metadata.aou_version AS STRING) AS aou_version,
    CAST(e.metadata.operating_system AS STRING) AS operating_system,
    CAST(e.metadata.user_agent AS STRING) AS user_agent,
    CAST(e.metadata.test_restarted AS BOOL) AS test_restarted,
    CAST(e.outcomes.any_timeouts AS INT64) AS any_timeouts,
    CAST(e.outcomes.flag_median_rtc AS INT64) AS emorecog_flag_median_rtc,
    CAST(e.outcomes.flag_same_response AS INT64) AS emorecog_flag_same_response,
    CAST(e.outcomes.flag_trial_flags AS INT64) AS emorecog_flag_trial_flags,
    CAST(e.outcomes.score AS FLOAT64) AS emorecog_score,
    CAST(e.outcomes.accuracy AS FLOAT64) AS emorecog_accuracy,
    CAST(e.outcomes.mean_rtc AS FLOAT64) AS emorecog_mean_rtc,
    CAST(e.outcomes.median_rtc AS FLOAT64) AS emorecog_median_rtc,
    CAST(e.outcomes.sd_rtc AS FLOAT64) AS emorecog_sd_rtc,
    SAFE_DIVIDE(CAST(e.outcomes.sd_rtc AS FLOAT64), CAST(e.outcomes.mean_rtc AS FLOAT64)) AS emorecog_cv_rtc,
    CAST(e.outcomes.happy_accuracy AS FLOAT64) AS emorecog_happy_accuracy,
    CAST(e.outcomes.happy_median_rtc AS FLOAT64) AS emorecog_happy_median_rtc,
    SAFE_DIVIDE(CAST(e.outcomes.happy_accuracy AS FLOAT64), CAST(e.outcomes.happy_median_rtc AS FLOAT64) / 1000.0) AS emorecog_happy_summary_eff,
    CAST(e.outcomes.angry_accuracy AS FLOAT64) AS emorecog_angry_accuracy,
    CAST(e.outcomes.angry_median_rtc AS FLOAT64) AS emorecog_angry_median_rtc,
    SAFE_DIVIDE(CAST(e.outcomes.angry_accuracy AS FLOAT64), CAST(e.outcomes.angry_median_rtc AS FLOAT64) / 1000.0) AS emorecog_angry_summary_eff,
    CAST(e.outcomes.fearful_accuracy AS FLOAT64) AS emorecog_fearful_accuracy,
    CAST(e.outcomes.fearful_median_rtc AS FLOAT64) AS emorecog_fearful_median_rtc,
    SAFE_DIVIDE(CAST(e.outcomes.fearful_accuracy AS FLOAT64), CAST(e.outcomes.fearful_median_rtc AS FLOAT64) / 1000.0) AS emorecog_fearful_summary_eff,
    CAST(e.outcomes.sad_accuracy AS FLOAT64) AS emorecog_sad_accuracy,
    CAST(e.outcomes.sad_median_rtc AS FLOAT64) AS emorecog_sad_median_rtc,
    SAFE_DIVIDE(CAST(e.outcomes.sad_accuracy AS FLOAT64), CAST(e.outcomes.sad_median_rtc AS FLOAT64) / 1000.0) AS emorecog_sad_summary_eff
  FROM {emo} AS e
  JOIN cohort AS c ON c.IID = e.person_id
  LEFT JOIN person_birth AS p ON p.person_id = e.person_id
  WHERE COALESCE(e.outcomes.flag_median_rtc, 0) = 0
    AND COALESCE(e.outcomes.flag_same_response, 0) = 0
    AND COALESCE(e.outcomes.flag_trial_flags, 0) = 0
    AND COALESCE(e.metadata.test_restarted, FALSE) = FALSE
),
emorecog_trial_rows AS (
  SELECT
    e.person_id,
    e.sitting_id,
    CASE
      WHEN LOWER(CAST(td.emotion AS STRING)) IN ('h', 'happy') THEN 'happy'
      WHEN LOWER(CAST(td.emotion AS STRING)) IN ('a', 'angry') THEN 'angry'
      WHEN LOWER(CAST(td.emotion AS STRING)) IN ('f', 'fearful') THEN 'fearful'
      WHEN LOWER(CAST(td.emotion AS STRING)) IN ('s', 'sad') THEN 'sad'
      ELSE NULL
    END AS emotion,
    SAFE_CAST(td.correct AS INT64) AS correct,
    SAFE_CAST(td.reaction_time AS FLOAT64) AS reaction_time
  FROM {emo} AS e
  JOIN emorecog_base AS b
    ON b.person_id = e.person_id
   AND b.sitting_id = e.sitting_id
  LEFT JOIN UNNEST(e.trial_data) AS td
),
emorecog_trial_agg AS (
  SELECT
    person_id,
    sitting_id,
    COUNTIF(emotion = 'happy') AS emorecog_happy_trials,
    SUM(IF(emotion = 'happy' AND correct = 1, 1, 0)) AS emorecog_happy_correct,
    SUM(IF(emotion = 'happy' AND reaction_time > 0, reaction_time, 0)) / 1000.0 AS emorecog_happy_rt_seconds,
    SAFE_DIVIDE(SUM(IF(emotion = 'happy' AND correct = 1, 1, 0)), SUM(IF(emotion = 'happy' AND reaction_time > 0, reaction_time, 0)) / 1000.0) AS emorecog_happy_rcs,
    COUNTIF(emotion = 'angry') AS emorecog_angry_trials,
    SUM(IF(emotion = 'angry' AND correct = 1, 1, 0)) AS emorecog_angry_correct,
    SUM(IF(emotion = 'angry' AND reaction_time > 0, reaction_time, 0)) / 1000.0 AS emorecog_angry_rt_seconds,
    SAFE_DIVIDE(SUM(IF(emotion = 'angry' AND correct = 1, 1, 0)), SUM(IF(emotion = 'angry' AND reaction_time > 0, reaction_time, 0)) / 1000.0) AS emorecog_angry_rcs,
    COUNTIF(emotion = 'fearful') AS emorecog_fearful_trials,
    SUM(IF(emotion = 'fearful' AND correct = 1, 1, 0)) AS emorecog_fearful_correct,
    SUM(IF(emotion = 'fearful' AND reaction_time > 0, reaction_time, 0)) / 1000.0 AS emorecog_fearful_rt_seconds,
    SAFE_DIVIDE(SUM(IF(emotion = 'fearful' AND correct = 1, 1, 0)), SUM(IF(emotion = 'fearful' AND reaction_time > 0, reaction_time, 0)) / 1000.0) AS emorecog_fearful_rcs,
    COUNTIF(emotion = 'sad') AS emorecog_sad_trials,
    SUM(IF(emotion = 'sad' AND correct = 1, 1, 0)) AS emorecog_sad_correct,
    SUM(IF(emotion = 'sad' AND reaction_time > 0, reaction_time, 0)) / 1000.0 AS emorecog_sad_rt_seconds,
    SAFE_DIVIDE(SUM(IF(emotion = 'sad' AND correct = 1, 1, 0)), SUM(IF(emotion = 'sad' AND reaction_time > 0, reaction_time, 0)) / 1000.0) AS emorecog_sad_rcs
  FROM emorecog_trial_rows
  GROUP BY person_id, sitting_id
),
emorecog AS (
  SELECT
    'emorecog' AS task,
    b.person_id AS IID,
    b.sitting_id,
    b.test_start_date_time,
    b.test_end_date_time,
    b.age_at_test,
    b.response_device,
    b.touch,
    b.test_language,
    b.test_version,
    b.aou_version,
    b.operating_system,
    b.user_agent,
    b.test_restarted,
    b.any_timeouts,
    CAST(NULL AS INT64) AS dd_flag_median_rt,
    CAST(NULL AS INT64) AS dd_flag_catch_trials,
    CAST(NULL AS INT64) AS gradcpt_flag_trial_flags,
    CAST(NULL AS INT64) AS gradcpt_flag_non_response,
    CAST(NULL AS INT64) AS gradcpt_flag_omission_error_rate,
    CAST(NULL AS INT64) AS flanker_flag_accuracy,
    CAST(NULL AS INT64) AS flanker_flag_trial_flags,
    CAST(NULL AS FLOAT64) AS dd_score,
    CAST(NULL AS FLOAT64) AS dd_catch_score,
    CAST(NULL AS FLOAT64) AS dd_lnk,
    CAST(NULL AS FLOAT64) AS dd_two_weeks_lnk,
    CAST(NULL AS FLOAT64) AS dd_one_month_lnk,
    CAST(NULL AS FLOAT64) AS dd_one_year_lnk,
    CAST(NULL AS FLOAT64) AS dd_ten_years_lnk,
    CAST(NULL AS FLOAT64) AS dd_mean_rt,
    CAST(NULL AS FLOAT64) AS dd_median_rt,
    CAST(NULL AS FLOAT64) AS dd_sd_rt,
    CAST(NULL AS FLOAT64) AS gradcpt_dprime,
    CAST(NULL AS FLOAT64) AS gradcpt_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_go_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_nogo_accuracy,
    CAST(NULL AS FLOAT64) AS gradcpt_score,
    CAST(NULL AS FLOAT64) AS gradcpt_crit,
    CAST(NULL AS FLOAT64) AS gradcpt_mean_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_median_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_sd_rtc,
    CAST(NULL AS FLOAT64) AS gradcpt_cv_rtc,
    CAST(NULL AS FLOAT64) AS flanker_score,
    CAST(NULL AS FLOAT64) AS flanker_accuracy,
    CAST(NULL AS FLOAT64) AS flanker_mean_rtc,
    CAST(NULL AS FLOAT64) AS flanker_median_rtc,
    CAST(NULL AS FLOAT64) AS flanker_sd_rtc,
    CAST(NULL AS FLOAT64) AS flanker_rcs_congruent,
    CAST(NULL AS FLOAT64) AS flanker_rcs_incongruent,
    CAST(NULL AS FLOAT64) AS flanker_accuracy_interference,
    CAST(NULL AS FLOAT64) AS flanker_median_rt_interference,
    CAST(NULL AS FLOAT64) AS flanker_rcs_interference,
    b.emorecog_flag_median_rtc,
    b.emorecog_flag_same_response,
    b.emorecog_flag_trial_flags,
    b.emorecog_score,
    b.emorecog_accuracy,
    b.emorecog_mean_rtc,
    b.emorecog_median_rtc,
    b.emorecog_sd_rtc,
    b.emorecog_cv_rtc,
    b.emorecog_happy_accuracy,
    b.emorecog_happy_median_rtc,
    b.emorecog_happy_summary_eff,
    t.emorecog_happy_rcs,
    t.emorecog_happy_trials,
    t.emorecog_happy_correct,
    t.emorecog_happy_rt_seconds,
    b.emorecog_angry_accuracy,
    b.emorecog_angry_median_rtc,
    b.emorecog_angry_summary_eff,
    t.emorecog_angry_rcs,
    t.emorecog_angry_trials,
    t.emorecog_angry_correct,
    t.emorecog_angry_rt_seconds,
    b.emorecog_fearful_accuracy,
    b.emorecog_fearful_median_rtc,
    b.emorecog_fearful_summary_eff,
    t.emorecog_fearful_rcs,
    t.emorecog_fearful_trials,
    t.emorecog_fearful_correct,
    t.emorecog_fearful_rt_seconds,
    b.emorecog_sad_accuracy,
    b.emorecog_sad_median_rtc,
    b.emorecog_sad_summary_eff,
    t.emorecog_sad_rcs,
    t.emorecog_sad_trials,
    t.emorecog_sad_correct,
    t.emorecog_sad_rt_seconds
  FROM emorecog_base AS b
  LEFT JOIN emorecog_trial_agg AS t
    ON t.person_id = b.person_id
   AND t.sitting_id = b.sitting_id
)
SELECT *
FROM delaydiscounting
UNION ALL
SELECT * FROM gradcpt
UNION ALL
SELECT * FROM flanker
UNION ALL
SELECT * FROM emorecog
ORDER BY task, IID, test_start_date_time, sitting_id
"""

    print(f"Querying valid ETM sittings from {args.etm_dataset}", flush=True)
    try:
        run_command(
            [
                "bq",
                "--project_id",
                google_project,
                "query",
                "--use_legacy_sql=false",
                "--format=csv",
                "--max_rows=1000000",
                query,
            ],
            stdout_path=extract_path,
        )
    finally:
        subprocess.run(["bq", "--project_id", google_project, "rm", "-f", temp_ref], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        temp_csv.unlink(missing_ok=True)

    return extract_path


def select_first_valid_sittings(valid: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    valid = valid.copy()
    valid["IID"] = valid["IID"].astype(str)
    valid["test_start_date_time"] = pd.to_datetime(valid["test_start_date_time"], errors="coerce", utc=True)
    valid["sitting_id_sort"] = pd.to_numeric(valid["sitting_id"], errors="coerce")
    repeat = (
        valid.groupby(["task", "IID"], as_index=False)
        .agg(n_valid_sittings=("sitting_id", "nunique"))
        .groupby("task", as_index=False)
        .agg(
            people_with_valid_sitting=("IID", "size"),
            people_with_2plus_valid_sittings=("n_valid_sittings", lambda x: int((x >= 2).sum())),
            mean_valid_sittings=("n_valid_sittings", "mean"),
            max_valid_sittings=("n_valid_sittings", "max"),
        )
    )
    repeat["frac_with_2plus_valid_sittings"] = repeat["people_with_2plus_valid_sittings"] / repeat["people_with_valid_sitting"]
    first = (
        valid.sort_values(["task", "IID", "test_start_date_time", "sitting_id_sort"], kind="mergesort")
        .drop_duplicates(["task", "IID"], keep="first")
        .drop(columns=["sitting_id_sort"])
    )
    return first, repeat


def transform_raw(values: pd.Series, spec: IndicatorSpec, eps: float | None = None) -> tuple[pd.Series, pd.Series]:
    raw = pd.to_numeric(values, errors="coerce").astype(float)
    invalid = pd.Series(False, index=raw.index)
    if spec.transform == "identity":
        out = raw.copy()
    elif spec.transform == "neg_identity":
        out = -raw
    elif spec.transform == "neg_log":
        invalid = raw <= 0
        out = pd.Series(np.nan, index=raw.index, dtype=float)
        out.loc[~invalid] = -np.log(raw.loc[~invalid])
    elif spec.transform == "log_plus_eps":
        require(eps is not None and eps > 0, f"positive eps required for {spec.name}")
        invalid = raw + eps <= 0
        out = pd.Series(np.nan, index=raw.index, dtype=float)
        out.loc[~invalid] = np.log(raw.loc[~invalid] + eps)
    elif spec.transform == "logit":
        clipped = raw.clip(LOGIT_EPS, 1.0 - LOGIT_EPS)
        invalid = (raw < 0) | (raw > 1)
        out = np.log(clipped / (1.0 - clipped))
    elif spec.transform == "log1p":
        invalid = raw < 0
        out = pd.Series(np.nan, index=raw.index, dtype=float)
        out.loc[~invalid] = np.log1p(raw.loc[~invalid])
    elif spec.transform == "scale_100":
        out = raw / 100.0
    else:
        raise ValueError(f"unknown transform {spec.transform}")
    return pd.Series(out, index=raw.index, dtype=float), invalid | raw.isna()


def prepare_indicators(
    df: pd.DataFrame,
    config: ScoreConfig,
    *,
    suffix: str = "",
) -> tuple[pd.DataFrame, list[dict[str, object]], list[dict[str, object]]]:
    params: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    indicators: dict[str, pd.Series] = {}

    eps = None
    if any(spec.transform == "log_plus_eps" for spec in config.indicators):
        vals = []
        for spec in config.indicators:
            if spec.transform == "log_plus_eps":
                x = pd.to_numeric(df[spec.source], errors="coerce")
                vals.append(x[x > 0])
        positive = pd.concat(vals) if vals else pd.Series(dtype=float)
        eps = float(0.5 * positive.min()) if len(positive) else np.nan
        require(np.isfinite(eps) and eps > 0, f"cannot compute positive RCS eps for {config.score_name}")
        params.append({"score_name": config.score_name, "parameter": f"rcs_eps{suffix}", "value": eps})

    for spec in config.indicators:
        transformed, invalid = transform_raw(df[spec.source], spec, eps)
        finite = transformed[np.isfinite(transformed)]
        if len(finite) == 0:
            lo = hi = np.nan
            winsor = transformed.copy()
        else:
            lo = float(np.nanquantile(finite, WINSOR_LO))
            hi = float(np.nanquantile(finite, WINSOR_HI))
            winsor = transformed.clip(lo, hi)
        z, mean, sd = zscore_array(winsor.to_numpy(dtype=float))
        indicators[spec.name] = pd.Series(z, index=df.index, dtype=float)
        params.extend(
            [
                {"score_name": config.score_name, "indicator": spec.name, "parameter": "source", "value": spec.source},
                {"score_name": config.score_name, "indicator": spec.name, "parameter": "transform", "value": spec.transform},
                {"score_name": config.score_name, "indicator": spec.name, "parameter": "winsor_lo", "value": lo},
                {"score_name": config.score_name, "indicator": spec.name, "parameter": "winsor_hi", "value": hi},
                {"score_name": config.score_name, "indicator": spec.name, "parameter": "mean_after_transform_winsor", "value": mean},
                {"score_name": config.score_name, "indicator": spec.name, "parameter": "sd_after_transform_winsor", "value": sd},
            ]
        )
        missing_rows.append(
            {
                "score_name": config.score_name,
                "task": config.task,
                "indicator": spec.name,
                "source": spec.source,
                "transform": spec.transform,
                "n_task_rows": len(df),
                "n_raw_missing_or_invalid": int(invalid.sum()),
                "n_finite_transformed": int(np.isfinite(transformed).sum()),
                "n_finite_z": int(np.isfinite(z).sum()),
            }
        )

    return pd.DataFrame(indicators, index=df.index), params, missing_rows


def drop_redundant_indicators(
    z: pd.DataFrame,
    config: ScoreConfig,
) -> tuple[list[str], pd.DataFrame, list[dict[str, object]]]:
    corr = z.corr(method="pearson", min_periods=20)
    corr_rows = []
    for a in corr.index:
        for b in corr.columns:
            corr_rows.append(
                {
                    "score_name": config.score_name,
                    "task": config.task,
                    "indicator_a": a,
                    "indicator_b": b,
                    "pearson": corr.loc[a, b],
                }
            )
    priority = {spec.name: spec.priority for spec in config.indicators}
    selected = set(z.columns)
    decisions = []
    pairs = []
    for i, a in enumerate(z.columns):
        for b in z.columns[i + 1 :]:
            r = corr.loc[a, b]
            if np.isfinite(r) and abs(r) > REDUNDANCY_R:
                pairs.append((abs(r), a, b, r))
    for _, a, b, r in sorted(pairs, reverse=True):
        if a not in selected or b not in selected:
            continue
        drop = b if priority[a] <= priority[b] else a
        keep = a if drop == b else b
        selected.remove(drop)
        decisions.append(
            {
                "score_name": config.score_name,
                "task": config.task,
                "kept_indicator": keep,
                "dropped_indicator": drop,
                "pearson": r,
                "reason": f"abs_r_gt_{REDUNDANCY_R}_fixed_priority",
            }
        )
    return [spec.name for spec in sorted(config.indicators, key=lambda s: s.priority) if spec.name in selected], pd.DataFrame(corr_rows), decisions


def fit_pca(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    pca = PCA(n_components=1, random_state=RNG_SEED)
    score = pca.fit_transform(X)[:, 0]
    loading = pca.components_[0].copy()
    orientation = 1.0 if np.nansum(loading) >= 0 else -1.0
    return score * orientation, loading * orientation, float(pca.explained_variance_ratio_[0])


def fit_fa(X: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, str, float, float]:
    try:
        fa = FactorAnalysis(n_components=1, random_state=seed, max_iter=2000, tol=1e-4)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            score = fa.fit_transform(X)[:, 0]
        if any(issubclass(w.category, ConvergenceWarning) for w in caught):
            raise RuntimeError("factor_analysis_nonconverged")
        loading = fa.components_[0].copy()
        orientation = 1.0 if np.nansum(loading) >= 0 else -1.0
        score *= orientation
        loading *= orientation
        observed = np.corrcoef(X, rowvar=False)
        model_cov = fa.get_covariance()
        # Inputs are standardized, so covariance residuals are interpretable as correlation residuals.
        resid = observed - model_cov
        offdiag = resid[~np.eye(resid.shape[0], dtype=bool)]
        max_abs_resid = float(np.nanmax(np.abs(offdiag))) if offdiag.size else np.nan
        mean_abs_resid = float(np.nanmean(np.abs(offdiag))) if offdiag.size else np.nan
        return score, loading, "factor_analysis", max_abs_resid, mean_abs_resid
    except Exception:
        score, loading, _ = fit_pca(X)
        return score, loading, "pca_fallback", np.nan, np.nan


def cronbach_alpha(X: np.ndarray) -> float:
    if X.shape[1] < 2:
        return np.nan
    item_vars = np.var(X, axis=0, ddof=1)
    total_var = np.var(X.sum(axis=1), ddof=1)
    if not np.isfinite(total_var) or total_var == 0:
        return np.nan
    k = X.shape[1]
    return float(k / (k - 1) * (1 - item_vars.sum() / total_var))


def design_age_sex(df: pd.DataFrame) -> np.ndarray:
    age = pd.to_numeric(df["age_at_test"], errors="coerce").to_numpy(dtype=float)
    sex = pd.to_numeric(df["sex_c"], errors="coerce").to_numpy(dtype=float)
    return np.column_stack([np.ones(len(df)), sex, age, age**2])


def residualize_z(raw: np.ndarray, meta: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, float, float]:
    x = design_age_sex(meta)
    beta = np.linalg.lstsq(x, raw, rcond=None)[0]
    resid = raw - x @ beta
    z, mean, sd = zscore_array(resid)
    return z, beta, mean, sd


def residualize_indicator_matrix(z: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    out = {}
    x = design_age_sex(meta)
    for col in z.columns:
        y = z[col].to_numpy(dtype=float)
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
        resid = y - x @ beta
        out[col] = zscore_array(resid)[0]
    return pd.DataFrame(out, index=z.index)


def make_score_frame(
    meta: pd.DataFrame,
    raw_score: np.ndarray,
    *,
    score_name: str,
    task: str,
    score_type: str,
    status: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    z_unadjusted, raw_mean, raw_sd = zscore_array(raw_score)
    z_age_sex, beta, resid_mean, resid_sd = residualize_z(raw_score, meta)
    out = meta[
        [
            "FID",
            "IID",
            "role",
            "fold_id",
            "ea_years",
            "teacher_z",
            "ses_ea_proxy_z",
            "sex_c",
            "task",
            "sitting_id",
            "test_start_date_time",
            "age_at_test",
            "response_device",
            "touch",
            "test_language",
            "test_version",
            "aou_version",
            "any_timeouts",
        ]
    ].copy()
    out["score_name"] = score_name
    out["score_type"] = score_type
    out["score_status"] = status
    out["score_raw"] = raw_score
    out["score_z_unadjusted"] = z_unadjusted
    out["score_z_age_sex"] = z_age_sex
    model = {
        "score_name": score_name,
        "task": task,
        "n": len(out),
        "status": status,
        "raw_mean": raw_mean,
        "raw_sd": raw_sd,
        "resid_intercept": float(beta[0]),
        "resid_sex_c": float(beta[1]),
        "resid_age_at_test": float(beta[2]),
        "resid_age_at_test_sq": float(beta[3]),
        "resid_mean": resid_mean,
        "resid_sd": resid_sd,
    }
    return out, model


def compute_factor_score(
    task_df: pd.DataFrame,
    config: ScoreConfig,
    *,
    seed: int,
    allow_grad_median_drop: bool = False,
    allow_emorecog_happy_drop: bool = False,
    dominance_threshold: float | None = None,
    force_pca: bool = False,
) -> dict[str, object]:
    z, params, missing_rows = prepare_indicators(task_df, config)
    selected, corr_df, redundancy = drop_redundant_indicators(z, config)
    selected_note = "initial"

    if allow_emorecog_happy_drop and "emorecog_happy_rcs_log" in selected:
        complete_with_happy = z[selected].apply(np.isfinite).all(axis=1)
        happy = z["emorecog_happy_rcs_log"].to_numpy(dtype=float)
        happy_finite = happy[np.isfinite(happy)]
        if int(complete_with_happy.sum()) < 20 or len(happy_finite) < 20 or float(np.nanstd(happy_finite, ddof=0)) == 0.0:
            selected = [x for x in selected if x != "emorecog_happy_rcs_log"]
            redundancy.append(
                {
                    "score_name": config.score_name,
                    "task": config.task,
                    "kept_indicator": ",".join(selected),
                    "dropped_indicator": "emorecog_happy_rcs_log",
                    "pearson": np.nan,
                    "reason": "happy_near_zero_or_missing_variance",
                }
            )
            selected_note = "emorecog_happy_dropped_before_fit"

    def fit_with_selected(selected_cols: list[str]) -> dict[str, object]:
        complete = z[selected_cols].apply(np.isfinite).all(axis=1) & np.isfinite(task_df["age_at_test"]) & np.isfinite(task_df["sex_c"])
        meta = task_df.loc[complete].copy()
        X = z.loc[complete, selected_cols].to_numpy(dtype=float)
        require(len(meta) >= 20, f"too few complete rows for {config.score_name}")
        pca_score, pca_loading, pca_var = fit_pca(X)
        if force_pca:
            fa_score, fa_loading, method, max_resid, mean_resid = pca_score, pca_loading, "pca_primary", np.nan, np.nan
        else:
            fa_score, fa_loading, method, max_resid, mean_resid = fit_fa(X, seed)
        unit_mean = X.mean(axis=1)
        determ_pearson, _, _ = pearson_spearman(fa_score, unit_mean)
        loadings = [
            {
                "score_name": config.score_name,
                "task": config.task,
                "method": method,
                "indicator": ind,
                "loading": float(load),
                "selected_for_primary": 1,
            }
            for ind, load in zip(selected_cols, fa_loading)
        ]
        pca_loadings = [
            {
                "score_name": config.score_name,
                "task": config.task,
                "method": "pca_diagnostic",
                "indicator": ind,
                "loading": float(load),
                "pc1_explained_variance_ratio": pca_var,
            }
            for ind, load in zip(selected_cols, pca_loading)
        ]
        weak = [ind for ind, load in zip(selected_cols, fa_loading) if abs(load) < WEAK_LOADING]
        wrong = [ind for ind, load in zip(selected_cols, fa_loading) if load < 0]
        return {
            "meta": meta,
            "raw_score": fa_score,
            "X": X,
            "z": z,
            "selected": selected_cols,
            "method": method,
            "loadings": loadings,
            "pca_loadings": pca_loadings,
            "weak": weak,
            "wrong": wrong,
            "complete_mask": complete,
            "diagnostics": {
                "score_name": config.score_name,
                "task": config.task,
                "n_complete": int(complete.sum()),
                "n_indicators": len(selected_cols),
                "method": method,
                "cronbach_alpha": cronbach_alpha(X),
                "factor_unit_mean_pearson": determ_pearson,
                "pc1_explained_variance_ratio": pca_var,
                "fa_max_abs_resid_corr_offdiag": max_resid,
                "fa_mean_abs_resid_corr_offdiag": mean_resid,
                "weak_loading_indicators": ",".join(weak),
                "wrong_sign_indicators": ",".join(wrong),
                "selected_indicators": ",".join(selected_cols),
            },
        }

    fit = fit_with_selected(selected)

    if allow_emorecog_happy_drop and "emorecog_happy_rcs_log" in fit["selected"]:
        load_map = {row["indicator"]: row["loading"] for row in fit["loadings"]}
        if abs(load_map.get("emorecog_happy_rcs_log", 0.0)) < WEAK_LOADING and len(fit["selected"]) > 3:
            selected = [x for x in fit["selected"] if x != "emorecog_happy_rcs_log"]
            redundancy.append(
                {
                    "score_name": config.score_name,
                    "task": config.task,
                    "kept_indicator": ",".join(selected),
                    "dropped_indicator": "emorecog_happy_rcs_log",
                    "pearson": np.nan,
                    "reason": "happy_loading_below_0.20",
                }
            )
            selected_note = "emorecog_happy_dropped_weak_loading"
            fit = fit_with_selected(selected)

    if allow_grad_median_drop and "gradcpt_median_rtc_neglog" in fit["selected"]:
        load_map = {row["indicator"]: row["loading"] for row in fit["loadings"]}
        if load_map.get("gradcpt_median_rtc_neglog", 0.0) < 0 and load_map.get("gradcpt_dprime", 0.0) > 0 and load_map.get("gradcpt_cv_rtc_neglog", 0.0) > 0:
            selected = [x for x in fit["selected"] if x != "gradcpt_median_rtc_neglog"]
            redundancy.append(
                {
                    "score_name": config.score_name,
                    "task": config.task,
                    "kept_indicator": ",".join(selected),
                    "dropped_indicator": "gradcpt_median_rtc_neglog",
                    "pearson": np.nan,
                    "reason": "median_rtc_loading_opposite_dprime_and_cv_rtc",
                }
            )
            selected_note = "grad_median_dropped_opposite_loading"
            fit = fit_with_selected(selected)

    status = "accepted"
    if fit["wrong"] or any(ind in fit["weak"] for ind in fit["selected"]):
        status = "flagged_loading_rule"
    dominance_ratio = np.nan
    if dominance_threshold is not None and fit["loadings"]:
        abs_load = np.asarray([abs(row["loading"]) for row in fit["loadings"]], dtype=float)
        denom = float(abs_load.sum())
        dominance_ratio = float(abs_load.max() / denom) if denom > 0 else np.nan
        if np.isfinite(dominance_ratio) and dominance_ratio > dominance_threshold:
            status = "flagged_dominant_indicator"
    fit["diagnostics"]["status"] = status
    fit["diagnostics"]["selected_note"] = selected_note
    fit["diagnostics"]["dominance_ratio"] = dominance_ratio
    score_df, age_model = make_score_frame(
        fit["meta"],
        fit["raw_score"],
        score_name=config.score_name,
        task=config.task,
        score_type=config.score_type,
        status=status,
    )
    score_df["n_indicators_used"] = len(fit["selected"])
    score_df["selected_indicators"] = ",".join(fit["selected"])
    score_df["score_method"] = fit["method"]
    fit["score_df"] = score_df
    fit["age_model"] = age_model
    fit["params"] = params
    fit["missing_rows"] = missing_rows
    fit["corr_df"] = corr_df
    fit["redundancy"] = redundancy
    return fit


def compute_unit_mean_score(task_df: pd.DataFrame, config: ScoreConfig) -> dict[str, object]:
    z, params, missing_rows = prepare_indicators(task_df, config)
    selected, corr_df, redundancy = drop_redundant_indicators(z, config)
    complete = z[selected].apply(np.isfinite).all(axis=1) & np.isfinite(task_df["age_at_test"]) & np.isfinite(task_df["sex_c"])
    meta = task_df.loc[complete].copy()
    X = z.loc[complete, selected].to_numpy(dtype=float)
    raw = X.mean(axis=1)
    status = "accepted"
    if config.score_name == "flanker_interference_unit_mean" and X.shape[1] == 2:
        r = float(np.corrcoef(X, rowvar=False)[0, 1])
        if not np.isfinite(r) or r < 0:
            status = "unstable_interference_indicators_opposite"
    score_df, age_model = make_score_frame(
        meta,
        raw,
        score_name=config.score_name,
        task=config.task,
        score_type=config.score_type,
        status=status,
    )
    score_df["n_indicators_used"] = len(selected)
    score_df["selected_indicators"] = ",".join(selected)
    score_df["score_method"] = "unit_mean"
    pca_score, pca_loading, pca_var = fit_pca(X)
    return {
        "score_df": score_df,
        "age_model": age_model,
        "params": params,
        "missing_rows": missing_rows,
        "corr_df": corr_df,
        "redundancy": redundancy,
        "loadings": [
            {
                "score_name": config.score_name,
                "task": config.task,
                "method": "unit_mean",
                "indicator": ind,
                "loading": 1.0 / math.sqrt(len(selected)),
                "selected_for_primary": 1,
            }
            for ind in selected
        ],
        "pca_loadings": [
            {
                "score_name": config.score_name,
                "task": config.task,
                "method": "pca_diagnostic",
                "indicator": ind,
                "loading": float(load),
                "pc1_explained_variance_ratio": pca_var,
            }
            for ind, load in zip(selected, pca_loading)
        ],
        "diagnostics": {
            "score_name": config.score_name,
            "task": config.task,
            "n_complete": int(complete.sum()),
            "n_indicators": len(selected),
            "method": "unit_mean",
            "cronbach_alpha": cronbach_alpha(X),
            "factor_unit_mean_pearson": 1.0,
            "pc1_explained_variance_ratio": pca_var,
            "fa_max_abs_resid_corr_offdiag": np.nan,
            "fa_mean_abs_resid_corr_offdiag": np.nan,
            "weak_loading_indicators": "",
            "wrong_sign_indicators": "",
            "selected_indicators": ",".join(selected),
            "status": status,
            "selected_note": "unit_mean_predeclared_split",
        },
    }


def compute_simple_score(task_df: pd.DataFrame, config: ScoreConfig) -> dict[str, object]:
    z, params, missing_rows = prepare_indicators(task_df, config)
    indicator = config.indicators[0].name
    complete = np.isfinite(z[indicator]) & np.isfinite(task_df["age_at_test"]) & np.isfinite(task_df["sex_c"])
    meta = task_df.loc[complete].copy()
    raw = z.loc[complete, indicator].to_numpy(dtype=float)
    score_df, age_model = make_score_frame(
        meta,
        raw,
        score_name=config.score_name,
        task=config.task,
        score_type=config.score_type,
        status="accepted",
    )
    score_df["n_indicators_used"] = 1
    score_df["selected_indicators"] = indicator
    score_df["score_method"] = "single_indicator"
    return {
        "score_df": score_df,
        "age_model": age_model,
        "params": params,
        "missing_rows": missing_rows,
        "corr_df": pd.DataFrame(
            [
                {
                    "score_name": config.score_name,
                    "task": config.task,
                    "indicator_a": indicator,
                    "indicator_b": indicator,
                    "pearson": 1.0,
                }
            ]
        ),
        "redundancy": [],
        "loadings": [
            {
                "score_name": config.score_name,
                "task": config.task,
                "method": "single_indicator",
                "indicator": indicator,
                "loading": 1.0,
                "selected_for_primary": 1,
            }
        ],
        "pca_loadings": [],
        "diagnostics": {
            "score_name": config.score_name,
            "task": config.task,
            "n_complete": int(complete.sum()),
            "n_indicators": 1,
            "method": "single_indicator",
            "cronbach_alpha": np.nan,
            "factor_unit_mean_pearson": 1.0,
            "pc1_explained_variance_ratio": np.nan,
            "fa_max_abs_resid_corr_offdiag": np.nan,
            "fa_mean_abs_resid_corr_offdiag": np.nan,
            "weak_loading_indicators": "",
            "wrong_sign_indicators": "",
            "selected_indicators": indicator,
            "status": "accepted",
            "selected_note": "single_indicator_simple",
        },
    }


def age_residualized_indicator_sensitivity(task_df: pd.DataFrame, config: ScoreConfig, primary_score: pd.DataFrame, seed: int) -> dict[str, object]:
    z, _, _ = prepare_indicators(task_df, config, suffix="_age_resid_sensitivity")
    selected, _, _ = drop_redundant_indicators(z, config)
    complete = z[selected].apply(np.isfinite).all(axis=1) & np.isfinite(task_df["age_at_test"]) & np.isfinite(task_df["sex_c"])
    meta = task_df.loc[complete].copy()
    if len(meta) < 20:
        return {"score_name": config.score_name, "task": config.task, "n": len(meta), "pearson_with_primary": np.nan, "note": "too_few_rows"}
    z_resid = residualize_indicator_matrix(z.loc[complete, selected], meta)
    if config.score_name == "emorecog_score_rt_factor":
        score, loading, _ = fit_pca(z_resid.to_numpy(dtype=float))
        method = "pca_primary"
    else:
        score, loading, method, _, _ = fit_fa(z_resid.to_numpy(dtype=float), seed)
    sensitivity_df, _ = make_score_frame(
        meta,
        score,
        score_name=f"{config.score_name}_age_resid_indicator_sensitivity",
        task=config.task,
        score_type="diagnostic",
        status="diagnostic",
    )
    merged = primary_score[["IID", "score_z_age_sex"]].merge(
        sensitivity_df[["IID", "score_z_age_sex"]],
        on="IID",
        suffixes=("_primary", "_age_resid_indicator"),
    )
    pearson, spearman, n = pearson_spearman(merged["score_z_age_sex_primary"], merged["score_z_age_sex_age_resid_indicator"])
    return {
        "score_name": config.score_name,
        "task": config.task,
        "n": n,
        "method": method,
        "pearson_with_primary": pearson,
        "spearman_with_primary": spearman,
        "loadings": ",".join(f"{ind}:{val:.6g}" for ind, val in zip(selected, loading)),
        "note": "factor_fit_on_age_sex_residualized_indicators",
    }


def score_correlations(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for score_name, sub in scores.groupby("score_name"):
        sub = sub.copy()
        sub["proxy_plus_teacher_z"] = zscore_array(sub["ses_ea_proxy_z"].to_numpy(dtype=float))[0] + zscore_array(sub["teacher_z"].to_numpy(dtype=float))[0]
        groups = [("combined", sub)]
        for role, role_df in sub.groupby("role"):
            groups.append((f"role_{role}", role_df))
        for group_name, g in groups:
            for target in ["ses_ea_proxy_z", "teacher_z", "ea_years", "proxy_plus_teacher_z"]:
                pearson, spearman, n = pearson_spearman(g["score_z_age_sex"], g[target])
                rows.append(
                    {
                        "score_name": score_name,
                        "task": str(g["task"].iloc[0]) if len(g) else "",
                        "group": group_name,
                        "target": target,
                        "n": n,
                        "pearson": pearson,
                        "spearman": spearman,
                    }
                )
    return pd.DataFrame(rows)


def simple_score_correlations(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    comparisons = {
        "dd_patience_factor": ["dd_simple_lnk", "dd_simple_score"],
        "gradcpt_perf_factor": ["gradcpt_simple_dprime", "gradcpt_simple_score"],
        "gradcpt_component_factor": ["gradcpt_simple_dprime", "gradcpt_simple_score"],
        "flanker_perf_factor": ["flanker_simple_score", "flanker_simple_rcs_interference"],
        "flanker_efficiency_unit_mean": ["flanker_simple_score"],
        "flanker_interference_unit_mean": ["flanker_simple_rcs_interference"],
        "emorecog_efficiency_factor": ["emorecog_simple_accuracy", "emorecog_simple_score", "emorecog_accuracy_factor"],
        "emorecog_summary_efficiency_factor": ["emorecog_simple_accuracy", "emorecog_simple_score"],
        "emorecog_accuracy_factor": ["emorecog_simple_accuracy", "emorecog_simple_score"],
        "emorecog_efficiency_unit_mean": ["emorecog_simple_accuracy", "emorecog_efficiency_factor"],
        "emorecog_speed_diagnostic": ["emorecog_simple_median_rtc"],
        "emorecog_score_rt_factor": ["emorecog_simple_score", "emorecog_simple_accuracy", "emorecog_simple_median_rtc", "emorecog_efficiency_factor"],
    }
    slim = scores[["IID", "score_name", "score_z_age_sex"]].copy()
    for score_name, simple_names in comparisons.items():
        base = slim.loc[slim["score_name"] == score_name, ["IID", "score_z_age_sex"]].rename(columns={"score_z_age_sex": "score"})
        for simple_name in simple_names:
            simp = slim.loc[slim["score_name"] == simple_name, ["IID", "score_z_age_sex"]].rename(columns={"score_z_age_sex": "simple"})
            merged = base.merge(simp, on="IID")
            pearson, spearman, n = pearson_spearman(merged["score"], merged["simple"])
            rows.append({"score_name": score_name, "simple_score_name": simple_name, "n": n, "pearson": pearson, "spearman": spearman})
    return pd.DataFrame(rows)


def choose_recommended_scores(scores: pd.DataFrame, simple_corr: pd.DataFrame, diagnostics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    recommended = []

    def add_alias(source: str, alias: str, reason: str) -> None:
        src = scores.loc[scores["score_name"] == source].copy()
        if src.empty:
            rows.append({"recommended_score": alias, "source_score": source, "reason": "source_missing"})
            return
        src["score_source"] = source
        src["score_name"] = alias
        src["score_type"] = "recommended"
        src["score_status"] = reason
        recommended.append(src)
        rows.append({"recommended_score": alias, "source_score": source, "reason": reason, "n": len(src)})

    def corr_ge(source: str, simple: str) -> bool:
        hit = simple_corr.loc[(simple_corr["score_name"] == source) & (simple_corr["simple_score_name"] == simple)]
        return bool(len(hit) and pd.notna(hit["pearson"].iloc[0]) and abs(hit["pearson"].iloc[0]) >= 0.95)

    # Final primary set after diagnostics:
    # - Delay Discounting uses the official aggregate -lnk. It summarizes all
    #   four delay-specific log-k fields and outperformed the four-delay factor
    #   against the proxy and teacher diagnostics.
    # - GradCPT uses the accepted PC1 fallback from dprime + RT consistency/speed.
    # - Flanker uses the predeclared efficiency score; interference remains a
    #   diagnostic because the split interference indicators were unstable.
    add_alias("dd_simple_lnk", "dd_patience", "official_lnk_primary")
    add_alias("gradcpt_simple_dprime" if corr_ge("gradcpt_perf_factor", "gradcpt_simple_dprime") else "gradcpt_perf_factor", "gradcpt_perf", "simple_ge_0.95_with_factor" if corr_ge("gradcpt_perf_factor", "gradcpt_simple_dprime") else "factor_primary")

    flanker_diag = diagnostics.loc[diagnostics["score_name"] == "flanker_perf_factor"]
    flanker_ok = bool(len(flanker_diag) and flanker_diag["status"].iloc[0] == "accepted")
    if flanker_ok:
        add_alias("flanker_simple_score" if corr_ge("flanker_perf_factor", "flanker_simple_score") else "flanker_perf_factor", "flanker_perf", "simple_ge_0.95_with_factor" if corr_ge("flanker_perf_factor", "flanker_simple_score") else "factor_primary")
    else:
        add_alias("flanker_simple_score" if corr_ge("flanker_efficiency_unit_mean", "flanker_simple_score") else "flanker_efficiency_unit_mean", "flanker_efficiency", "simple_ge_0.95_with_efficiency" if corr_ge("flanker_efficiency_unit_mean", "flanker_simple_score") else "split_efficiency")

    emorecog_diag = diagnostics.loc[diagnostics["score_name"] == "emorecog_score_rt_factor"]
    emorecog_has_score_rt = bool(len(emorecog_diag))
    if emorecog_has_score_rt:
        add_alias(
            "emorecog_score_rt_factor",
            "emorecog_perf",
            "score_cv_median_pc1_primary",
        )
    else:
        unit_diag = diagnostics.loc[diagnostics["score_name"] == "emorecog_efficiency_unit_mean"]
        unit_ok = bool(len(unit_diag) and unit_diag["status"].iloc[0] == "accepted")
        add_alias(
            "emorecog_efficiency_unit_mean" if unit_ok else "emorecog_simple_accuracy",
            "emorecog_perf",
            "efficiency_unit_mean_fallback" if unit_ok else "simple_accuracy_fallback",
        )

    rec_df = pd.concat(recommended, ignore_index=True) if recommended else pd.DataFrame()
    return rec_df, pd.DataFrame(rows)


def recommended_wide_table(scores: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    base_cols = ["FID", "IID", "role", "fold_id", "ea_years", "teacher_z", "ses_ea_proxy_z", "sex_c"]
    base = cohort[base_cols].drop_duplicates("IID").copy()
    for score_name in RECOMMENDED_SCORE_NAMES:
        keep_cols = ["IID", "score_z_age_sex", "sitting_id", "age_at_test"]
        if score_name == "emorecog_perf":
            keep_cols.extend(
                [
                    "score_raw",
                    "score_source",
                    "n_indicators_used",
                    "selected_indicators",
                    "score_method",
                    "score_status",
                    "any_timeouts",
                ]
            )
        keep_cols = [col for col in keep_cols if col in scores.columns]
        sub = scores.loc[
            (scores["score_name"] == score_name) & (scores["score_type"] == "recommended"),
            keep_cols,
        ].copy()
        rename = {
            "score_z_age_sex": f"{score_name}_z_age_sex",
            "sitting_id": f"{score_name}_sitting_id",
            "age_at_test": f"{score_name}_age_at_test",
        }
        if score_name == "emorecog_perf":
            rename.update(
                {
                    "score_raw": "emorecog_perf_raw",
                    "score_source": "emorecog_score_source",
                    "n_indicators_used": "emorecog_n_indicators_used",
                    "selected_indicators": "emorecog_selected_indicators",
                    "score_method": "emorecog_score_method",
                    "score_status": "emorecog_score_status",
                    "any_timeouts": "emorecog_any_timeouts",
                    "sitting_id": "emorecog_valid_sitting_id",
                    "age_at_test": "emorecog_age_at_test",
                }
            )
        sub = sub.rename(columns=rename)
        base = base.merge(sub, on="IID", how="left", validate="one_to_one")
    return base


def admin_sensitivity(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for score_name, sub in scores.groupby("score_name"):
        if score_name not in set(RECOMMENDED_SCORE_NAMES):
            continue
        sub = sub.copy()
        y = sub["score_raw"].to_numpy(dtype=float)
        base_z, _, _, _ = residualize_z(y, sub)
        cats = sub[["response_device", "touch", "test_language", "test_version"]].copy()
        for col in cats.columns:
            cats[col] = cats[col].astype("string").fillna("missing")
            counts = cats[col].value_counts()
            cats[col] = cats[col].where(cats[col].map(counts) >= 20, "other_small")
        dummies = pd.get_dummies(cats, drop_first=True, dtype=float)
        age_sex = pd.DataFrame(
            {
                "intercept": 1.0,
                "sex_c": pd.to_numeric(sub["sex_c"], errors="coerce"),
                "age": pd.to_numeric(sub["age_at_test"], errors="coerce"),
            },
            index=sub.index,
        )
        age_sex["age_sq"] = age_sex["age"] ** 2
        design = pd.concat([age_sex, dummies], axis=1)
        ok = design.apply(np.isfinite).all(axis=1) & np.isfinite(y)
        if ok.sum() < 20:
            rows.append({"score_name": score_name, "n": int(ok.sum()), "pearson_age_sex_vs_admin": np.nan, "note": "too_few_rows"})
            continue
        beta = np.linalg.lstsq(design.loc[ok].to_numpy(dtype=float), y[ok], rcond=None)[0]
        resid = y[ok] - design.loc[ok].to_numpy(dtype=float) @ beta
        admin_z = zscore_array(resid)[0]
        pearson, spearman, n = pearson_spearman(pd.Series(base_z, index=sub.index).loc[ok], admin_z)
        rows.append(
            {
                "score_name": score_name,
                "task": str(sub["task"].iloc[0]),
                "n": n,
                "admin_covariate_count": dummies.shape[1],
                "pearson_age_sex_vs_admin": pearson,
                "spearman_age_sex_vs_admin": spearman,
            }
        )
    return pd.DataFrame(rows)


def recommended_cross_task_correlations(wide: pd.DataFrame, work_dir: Path) -> pd.DataFrame:
    score_cols = {
        "dd_patience": "dd_patience_z_age_sex",
        "gradcpt_perf": "gradcpt_perf_z_age_sex",
        "flanker_efficiency": "flanker_efficiency_z_age_sex",
        "emorecog_perf": "emorecog_perf_z_age_sex",
    }
    etmg = work_dir / "etm_general_factor" / "etm_general_factor_scores_wide.tsv"
    if etmg.exists() and etmg.stat().st_size > 0:
        g = pd.read_csv(etmg, sep="\t", dtype={"IID": str}, usecols=["IID", "etm_g_z"])
        wide = wide.merge(g, on="IID", how="left", validate="one_to_one")
        score_cols["etm_g_three_domain"] = "etm_g_z"
    rows = []
    names = [name for name, col in score_cols.items() if col in wide.columns]
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            lcol = score_cols[left]
            rcol = score_cols[right]
            pearson, spearman, n = pearson_spearman(wide[lcol], wide[rcol])
            rows.append(
                {
                    "score_a": left,
                    "score_b": right,
                    "column_a": lcol,
                    "column_b": rcol,
                    "n": n,
                    "pearson": pearson,
                    "spearman": spearman,
                }
            )
    return pd.DataFrame(rows)


def emorecog_qc_counts(valid: pd.DataFrame, first: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, df in [("all_valid_sittings", valid.loc[valid["task"] == "emorecog"]), ("first_valid_sittings", first.loc[first["task"] == "emorecog"])]:
        row = {
            "subset": label,
            "n_rows": int(len(df)),
            "n_iids": int(df["IID"].nunique()) if "IID" in df.columns else 0,
            "n_any_timeouts": int(pd.to_numeric(df.get("any_timeouts", pd.Series(dtype=float)), errors="coerce").fillna(0).astype(float).gt(0).sum()),
        }
        for flag in ["emorecog_flag_median_rtc", "emorecog_flag_same_response", "emorecog_flag_trial_flags"]:
            if flag in df.columns:
                row[f"{flag}_sum"] = int(pd.to_numeric(df[flag], errors="coerce").fillna(0).sum())
        for emotion in ["happy", "angry", "fearful", "sad"]:
            rcs = f"emorecog_{emotion}_rcs"
            trials = f"emorecog_{emotion}_trials"
            if rcs in df.columns:
                row[f"{emotion}_rcs_nonmissing"] = int(pd.to_numeric(df[rcs], errors="coerce").notna().sum())
            if trials in df.columns:
                row[f"{emotion}_trial_count_sum"] = float(pd.to_numeric(df[trials], errors="coerce").sum())
        rows.append(row)
    return pd.DataFrame(rows)


def emorecog_timeout_sensitivity(task_df: pd.DataFrame, recommended_scores: pd.DataFrame, seed: int) -> pd.DataFrame:
    primary = recommended_scores.loc[recommended_scores["score_name"] == "emorecog_perf"].copy()
    if task_df.empty or primary.empty or "any_timeouts" not in task_df.columns:
        return pd.DataFrame([{"n": 0, "note": "missing_emorecog_task_or_primary"}])
    no_timeout = task_df.loc[pd.to_numeric(task_df["any_timeouts"], errors="coerce").fillna(0).eq(0)].copy()
    if len(no_timeout) < 20:
        return pd.DataFrame([{"n": int(len(no_timeout)), "note": "too_few_no_timeout_rows"}])
    try:
        result = compute_factor_score(
            no_timeout,
            EMORECOG_SCORE_RT,
            seed=seed,
            force_pca=True,
        )
    except SystemExit as exc:
        return pd.DataFrame([{"n": int(len(no_timeout)), "note": f"timeout_exclusion_fit_failed:{exc}"}])
    sens = result["score_df"]
    merged = primary[["IID", "score_z_age_sex"]].merge(
        sens[["IID", "score_z_age_sex"]],
        on="IID",
        suffixes=("_primary", "_no_timeout_sensitivity"),
    )
    pearson, spearman, n = pearson_spearman(merged["score_z_age_sex_primary"], merged["score_z_age_sex_no_timeout_sensitivity"])
    return pd.DataFrame(
        [
            {
                "n_all_primary": int(len(primary)),
                "n_no_timeout_fit": int(len(no_timeout)),
                "n_overlap": n,
                "pearson_primary_vs_no_timeout_sensitivity": pearson,
                "spearman_primary_vs_no_timeout_sensitivity": spearman,
                "no_timeout_method": result["diagnostics"].get("method"),
                "no_timeout_status": result["diagnostics"].get("status"),
                "no_timeout_selected_indicators": result["diagnostics"].get("selected_indicators"),
            }
        ]
    )


def write_outputs(args: argparse.Namespace, outputs: dict[str, pd.DataFrame]) -> None:
    args.work_dir.mkdir(parents=True, exist_ok=True)
    for name, df in outputs.items():
        path = args.work_dir / f"etm_cog_task_factors_{name}.tsv"
        df.to_csv(path, sep="\t", index=False)
        print(f"Wrote {path}", flush=True)

    if args.stage_aggregate:
        aggregate = {
            "summary",
            "correlations",
            "simple_score_correlations",
            "factor_diagnostics",
            "loadings",
            "pca_loadings",
            "indicator_correlations",
            "redundancy_decisions",
            "indicator_missingness",
            "age_sex_models",
            "repeat_sittings",
            "recommended_sources",
            "admin_sensitivity",
            "age_residualized_indicator_sensitivity",
            "params",
            "recommended_cross_task_correlations",
            "emorecog_indicator_correlations",
            "emorecog_fa_loadings",
            "emorecog_pca_loadings",
            "emorecog_factor_vs_simple_correlations",
            "emorecog_qc_counts",
            "emorecog_redundancy_decisions",
            "emorecog_age_sex_diagnostics",
            "emorecog_device_language_version_diagnostics",
            "emorecog_timeout_sensitivity",
        }
        args.workspace_scrap_dir.mkdir(parents=True, exist_ok=True)
        for name in aggregate:
            src = args.work_dir / f"etm_cog_task_factors_{name}.tsv"
            if src.exists():
                dest = args.workspace_scrap_dir / src.name
                shutil.copyfile(src, dest)
                print(f"Staged aggregate {dest}", flush=True)


def summarize_scores(scores: pd.DataFrame) -> pd.DataFrame:
    return (
        scores.groupby(["score_name", "score_type", "score_status", "task"], as_index=False)
        .agg(
            n=("IID", "size"),
            unique_iids=("IID", "nunique"),
            age_mean=("age_at_test", "mean"),
            age_min=("age_at_test", "min"),
            age_max=("age_at_test", "max"),
            z_mean=("score_z_age_sex", "mean"),
            z_sd=("score_z_age_sex", lambda x: float(np.nanstd(x, ddof=0))),
        )
        .sort_values(["task", "score_name"])
    )


def main() -> None:
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    print("=== ETM cognitive task factor scoring ===", flush=True)
    print(f"ETM dataset: {args.etm_dataset}", flush=True)
    print(f"Workspace CDR: {args.workspace_cdr}", flush=True)
    print(f"Temp BQ dataset: {args.bq_temp_dataset}", flush=True)
    print(f"SES-EA proxy dir: {args.ses_ea_dir}", flush=True)

    cohort = load_proxy_cohort(args.ses_ea_dir)
    print(f"Proxy cohort rows: {len(cohort)}", flush=True)
    print(cohort["role"].value_counts(dropna=False).to_string(), flush=True)
    bad_sex = cohort.loc[~cohort["sex_c"].isin([-0.5, 0.5])]
    require(bad_sex.empty, "sex_c is not the expected confirmed genetic-sex centered binary covariate")

    extract_path = ensure_bq_extract(args, cohort)
    valid = pd.read_csv(extract_path, dtype={"IID": str})
    require(not valid.empty, "ETM extract returned no rows")
    first, repeat = select_first_valid_sittings(valid)
    first = first.merge(cohort, on="IID", how="inner", validate="many_to_one")
    first["age_at_test"] = pd.to_numeric(first["age_at_test"], errors="coerce")
    missing_age = int(first["age_at_test"].isna().sum())
    print(f"First valid task rows in proxy cohort: {len(first)}", flush=True)
    print(f"Task rows missing age_at_test from birth_datetime/test_start_date_time: {missing_age}", flush=True)

    all_results: list[dict[str, object]] = []
    scores = []
    loadings = []
    pca_loadings = []
    corr_mats = []
    redundancy = []
    missingness = []
    params = [
        {"parameter": "etm_dataset", "value": args.etm_dataset},
        {"parameter": "workspace_cdr", "value": args.workspace_cdr},
        {"parameter": "cohort_source", "value": str(args.ses_ea_dir / "all_scores.tsv")},
        {"parameter": "sex_source", "value": str(args.ses_ea_dir / "base_covar.txt")},
        {"parameter": "age_source", "value": "DATE_DIFF(test_start_date_time, person.birth_datetime) / 365.25"},
        {"parameter": "valid_sitting_sort", "value": "test_start_date_time,sitting_id"},
        {"parameter": "winsor_quantiles", "value": f"{WINSOR_LO},{WINSOR_HI}"},
        {"parameter": "redundancy_abs_r_threshold", "value": REDUNDANCY_R},
        {"parameter": "weak_loading_threshold", "value": WEAK_LOADING},
    ]
    age_models = []
    age_resid_sens = []

    task_map = {task: df.copy() for task, df in first.groupby("task")}

    # Primary/sensitivity factor scores.
    for config in [DD_PRIMARY, GRADCPT_PRIMARY, GRADCPT_COMPONENT, FLANKER_PRIMARY, EMORECOG_SCORE_RT, EMORECOG_EFFICIENCY, EMORECOG_SUMMARY_EFFICIENCY, EMORECOG_ACCURACY, EMORECOG_SPEED]:
        task_df = task_map.get(config.task, pd.DataFrame())
        require(not task_df.empty, f"no first-valid rows for task {config.task}")
        result = compute_factor_score(
            task_df,
            config,
            seed=args.seed,
            allow_grad_median_drop=(config.score_name == "gradcpt_perf_factor"),
            allow_emorecog_happy_drop=(config.score_name == "emorecog_efficiency_factor"),
            dominance_threshold=0.60 if config.score_name == "emorecog_efficiency_factor" else None,
            force_pca=(config.score_name == "emorecog_score_rt_factor"),
        )
        all_results.append(result)
        scores.append(result["score_df"])
        loadings.extend(result["loadings"])
        pca_loadings.extend(result["pca_loadings"])
        corr_mats.append(result["corr_df"])
        redundancy.extend(result["redundancy"])
        missingness.extend(result["missing_rows"])
        params.extend(result["params"])
        age_models.append(result["age_model"])
        age_resid_sens.append(age_residualized_indicator_sensitivity(task_df, config, result["score_df"], args.seed))

    # Flanker split scores are always computed; recommendation chooses them only if one-factor is rejected.
    for config in [FLANKER_EFFICIENCY, FLANKER_INTERFERENCE, EMORECOG_EFFICIENCY_UNIT_MEAN]:
        task_df = task_map.get(config.task, pd.DataFrame())
        require(not task_df.empty, f"no first-valid rows for task {config.task}")
        result = compute_unit_mean_score(task_df, config)
        all_results.append(result)
        scores.append(result["score_df"])
        loadings.extend(result["loadings"])
        pca_loadings.extend(result["pca_loadings"])
        corr_mats.append(result["corr_df"])
        redundancy.extend(result["redundancy"])
        missingness.extend(result["missing_rows"])
        params.extend(result["params"])
        age_models.append(result["age_model"])

    # Simple official/validation scores.
    for config in SIMPLE_SCORES:
        task_df = task_map.get(config.task, pd.DataFrame())
        result = compute_simple_score(task_df, config)
        all_results.append(result)
        scores.append(result["score_df"])
        loadings.extend(result["loadings"])
        corr_mats.append(result["corr_df"])
        missingness.extend(result["missing_rows"])
        params.extend(result["params"])
        age_models.append(result["age_model"])

    factor_diag = pd.DataFrame([r["diagnostics"] for r in all_results])
    all_scores = pd.concat(scores, ignore_index=True)
    simple_corr = simple_score_correlations(all_scores)
    recommended_scores, recommended_sources = choose_recommended_scores(all_scores, simple_corr, factor_diag)
    if not recommended_scores.empty:
        all_scores = pd.concat([all_scores, recommended_scores], ignore_index=True)
    recommended_wide = recommended_wide_table(all_scores, cohort)
    loadings_df = pd.DataFrame(loadings)
    pca_loadings_df = pd.DataFrame(pca_loadings)
    indicator_corr_df = pd.concat(corr_mats, ignore_index=True) if corr_mats else pd.DataFrame()
    redundancy_df = pd.DataFrame(
        redundancy,
        columns=["score_name", "task", "kept_indicator", "dropped_indicator", "pearson", "reason"],
    )
    indicator_missingness_df = pd.DataFrame(missingness)
    age_sex_models_df = pd.DataFrame(age_models)
    admin_sensitivity_df = admin_sensitivity(all_scores)
    age_resid_sens_df = pd.DataFrame(age_resid_sens)

    def emorecog_filter(df: pd.DataFrame, *cols: str) -> pd.DataFrame:
        if df.empty:
            return df
        mask = pd.Series(False, index=df.index)
        for col in cols:
            if col in df.columns:
                mask |= df[col].astype(str).str.contains("emorecog", na=False)
        return df.loc[mask].copy()

    outputs = {
        "scores": all_scores,
        "recommended_wide": recommended_wide,
        "summary": summarize_scores(all_scores),
        "correlations": score_correlations(all_scores),
        "simple_score_correlations": simple_corr,
        "factor_diagnostics": factor_diag,
        "loadings": loadings_df,
        "pca_loadings": pca_loadings_df,
        "indicator_correlations": indicator_corr_df,
        "redundancy_decisions": redundancy_df,
        "indicator_missingness": indicator_missingness_df,
        "age_sex_models": age_sex_models_df,
        "repeat_sittings": repeat,
        "recommended_sources": recommended_sources,
        "admin_sensitivity": admin_sensitivity_df,
        "age_residualized_indicator_sensitivity": age_resid_sens_df,
        "params": pd.DataFrame(params),
        "recommended_cross_task_correlations": recommended_cross_task_correlations(recommended_wide, args.work_dir),
        "emorecog_indicator_correlations": emorecog_filter(indicator_corr_df, "score_name", "task"),
        "emorecog_fa_loadings": emorecog_filter(loadings_df, "score_name", "task"),
        "emorecog_pca_loadings": emorecog_filter(pca_loadings_df, "score_name", "task"),
        "emorecog_factor_vs_simple_correlations": emorecog_filter(simple_corr, "score_name", "simple_score_name"),
        "emorecog_qc_counts": emorecog_qc_counts(valid, first),
        "emorecog_redundancy_decisions": emorecog_filter(redundancy_df, "score_name", "task"),
        "emorecog_age_sex_diagnostics": emorecog_filter(pd.concat([age_sex_models_df, age_resid_sens_df], ignore_index=True, sort=False), "score_name", "task"),
        "emorecog_device_language_version_diagnostics": emorecog_filter(admin_sensitivity_df, "score_name", "task"),
        "emorecog_timeout_sensitivity": emorecog_timeout_sensitivity(task_map.get("emorecog", pd.DataFrame()), recommended_scores, args.seed),
    }
    write_outputs(args, outputs)

    print("\nRecommended score sources:", flush=True)
    print(outputs["recommended_sources"].to_string(index=False), flush=True)
    print("\nRecommended score correlations:", flush=True)
    rec_corr = outputs["correlations"].loc[outputs["correlations"]["score_name"].isin(set(RECOMMENDED_SCORE_NAMES))]
    print(rec_corr.to_string(index=False), flush=True)
    print("\n=== ETM cognitive task factor scoring complete ===", flush=True)


if __name__ == "__main__":
    main()
