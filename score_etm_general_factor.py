#!/usr/bin/env python3
"""Build ETM general cognitive/performance factor scores.

This command consumes already-built task-specific ETM scores. It does not query
ETM tables, run GWAS, or choose the phenotype based on SES-EA/teacher targets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shutil
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy
from sklearn import __version__ as sklearn_version
from sklearn.decomposition import FactorAnalysis, PCA
from sklearn.exceptions import ConvergenceWarning


THREE_TASK_LABELS = ("DD", "GRADCPT", "FLANKER")
FOUR_TASK_LABELS = ("DD", "GRADCPT", "FLANKER", "EMORECOG")
TASK_LABELS = THREE_TASK_LABELS
DEFAULT_TASK_COLUMNS = {
    "DD": "dd_patience_z_age_sex",
    "GRADCPT": "gradcpt_perf_z_age_sex",
    "FLANKER": "flanker_efficiency_z_age_sex",
    "EMORECOG": "emorecog_perf_z_age_sex",
}
ALTERNATE_FLANKER_COLUMN = "flanker_perf_z_age_sex"
MIN_GROUP_N = 20
HASH_MAX_BYTES = 100 * 1024 * 1024
UNIQUENESS_FLOOR = 1e-8


@dataclass(frozen=True)
class FactorScoreResult:
    labels: tuple[str, ...]
    columns: tuple[str, ...]
    means: np.ndarray
    sds: np.ndarray
    complete_case_mask: np.ndarray
    standardized: np.ndarray
    loadings: np.ndarray
    uniquenesses: np.ndarray
    uniquenesses_for_scoring: np.ndarray
    floored_uniqueness: np.ndarray
    orientation_sign: float
    weights_by_pattern: pd.DataFrame
    g_hat: np.ndarray
    g_z_candidate: np.ndarray
    g_hat_mean_complete_case: float
    g_hat_sd_complete_case: float
    observed_corr_complete_case: np.ndarray
    model_corr: np.ndarray
    residual_corr: np.ndarray
    pca_loadings: np.ndarray
    pca_explained_variance_ratio: float
    fa_n_iter: int | None
    fa_loglike_last: float | None
    warning_messages: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-score-file", type=Path, required=True)
    parser.add_argument("--all-scores-file", type=Path, required=True)
    parser.add_argument("--base-covar-file", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workspace-scrap-dir", type=Path, default=None)
    parser.add_argument("--stage-aggregate", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--flanker-input",
        choices=("auto", "flanker_efficiency_z_age_sex", "flanker_perf_z_age_sex"),
        default="auto",
    )
    parser.add_argument("--min-complete-case-n", type=int, default=500)
    parser.add_argument("--force-three-domain-g", action="store_true")
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--uniqueness-floor", type=float, default=UNIQUENESS_FLOOR)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"ERROR: {message}")


def read_header(path: Path) -> list[str]:
    require(path.exists() and path.stat().st_size > 0, f"missing or empty input file: {path}")
    with path.open() as handle:
        return handle.readline().rstrip("\n").split("\t")


def sha256_file(path: Path) -> str | None:
    size = path.stat().st_size
    if size > HASH_MAX_BYTES:
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_fingerprint(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": sha256_file(path),
    }


def resolve_flanker_column(task_header: Iterable[str], requested: str) -> tuple[str, str]:
    cols = set(task_header)
    if requested == "auto":
        if DEFAULT_TASK_COLUMNS["FLANKER"] in cols:
            return DEFAULT_TASK_COLUMNS["FLANKER"], "auto:flanker_efficiency_z_age_sex"
        if ALTERNATE_FLANKER_COLUMN in cols:
            return ALTERNATE_FLANKER_COLUMN, "auto:flanker_perf_z_age_sex"
        raise SystemExit("ERROR: no supported Flanker input column found in task-score file")
    require(requested in cols, f"requested --flanker-input column is absent: {requested}")
    return requested, f"explicit:{requested}"


def build_run_signature(args: argparse.Namespace, task_columns: dict[str, str], flanker_source: str) -> tuple[str, pd.DataFrame]:
    base_covar_file = args.base_covar_file or (args.all_scores_file.parent / "base_covar.txt")
    task_order = [label for label in FOUR_TASK_LABELS if label in task_columns]
    payload = {
        "task_score_file": file_fingerprint(args.task_score_file),
        "all_scores_file": file_fingerprint(args.all_scores_file),
        "base_covar_file": file_fingerprint(base_covar_file),
        "task_order": task_order,
        "task_columns": task_columns,
        "flanker_input_source": flanker_source,
        "four_domain_enabled": "EMORECOG" in task_columns,
        "min_complete_case_n": args.min_complete_case_n,
        "force_three_domain_g": bool(args.force_three_domain_g),
        "random_state": args.random_state,
        "uniqueness_floor": args.uniqueness_floor,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    signature = hashlib.sha256(blob.encode("utf-8")).hexdigest()
    rows = [{"parameter": "run_signature", "value": signature}]
    for key, value in payload.items():
        rows.append({"parameter": key, "value": json.dumps(value, sort_keys=True)})
    return signature, pd.DataFrame(rows)


def saved_run_signature(output_dir: Path) -> str | None:
    params = output_dir / "diagnostics" / "etm_g_reproducibility_params.tsv"
    if not params.exists():
        return None
    try:
        df = pd.read_csv(params, sep="\t", dtype=str)
    except Exception:
        return None
    if {"parameter", "value"} - set(df.columns):
        return None
    hit = df.loc[df["parameter"] == "run_signature", "value"]
    if hit.empty:
        return None
    return str(hit.iloc[0])


def check_idempotency(args: argparse.Namespace, signature: str) -> None:
    wide = args.output_dir / "etm_general_factor_scores_wide.tsv"
    params = args.output_dir / "diagnostics" / "etm_g_reproducibility_params.tsv"
    if args.force:
        return
    if not wide.exists() and not params.exists():
        return
    saved = saved_run_signature(args.output_dir)
    if saved == signature and wide.exists():
        print(f"Existing ETM-g outputs match inputs/parameters; skipping: {args.output_dir}", flush=True)
        raise SystemExit(0)
    raise SystemExit(
        "ERROR: ETM-g outputs already exist but do not match the current inputs/parameters. "
        "Rerun with --force to overwrite."
    )


def load_inputs(args: argparse.Namespace, task_columns: dict[str, str]) -> pd.DataFrame:
    base_covar_file = args.base_covar_file or (args.all_scores_file.parent / "base_covar.txt")
    require(args.all_scores_file.exists() and args.all_scores_file.stat().st_size > 0, f"missing {args.all_scores_file}")
    require(base_covar_file.exists() and base_covar_file.stat().st_size > 0, f"missing {base_covar_file}")
    require(args.task_score_file.exists() and args.task_score_file.stat().st_size > 0, f"missing {args.task_score_file}")

    all_scores = pd.read_csv(args.all_scores_file, sep="\t", dtype={"FID": str, "IID": str, "role": str, "fold_id": str})
    required_scores = {"FID", "IID", "role", "fold_id", "ea_years", "teacher_z", "ses_ea_proxy_z"}
    missing_scores = required_scores - set(all_scores.columns)
    require(not missing_scores, f"{args.all_scores_file} missing columns: {sorted(missing_scores)}")
    require(all_scores["IID"].is_unique, "all_scores.tsv has duplicate IIDs")

    base_header = set(read_header(base_covar_file))
    base_cols = ["IID"]
    for col in ("sex_c", "yob_c", "fractional_yob", "yob", "birth_year"):
        if col in base_header:
            base_cols.append(col)
    require("sex_c" in base_cols, f"{base_covar_file} missing required sex_c column")
    covar = pd.read_csv(base_covar_file, sep="\t", dtype={"IID": str}, usecols=base_cols)
    require(covar["IID"].is_unique, "base_covar file has duplicate IIDs")

    task_usecols = ["IID"] + list(task_columns.values())
    task_header = set(read_header(args.task_score_file))
    if "person_id" in task_header:
        task_usecols.append("person_id")
    missing_task = set(task_usecols) - task_header
    require(not missing_task, f"{args.task_score_file} missing columns: {sorted(missing_task)}")
    tasks = pd.read_csv(args.task_score_file, sep="\t", dtype={"IID": str}, usecols=task_usecols)
    require(tasks["IID"].is_unique, "task-score file has duplicate IIDs")

    df = all_scores.merge(covar, on="IID", how="left", validate="one_to_one")
    df = df.merge(tasks, on="IID", how="left", validate="one_to_one")
    if "person_id" not in df.columns:
        df["person_id"] = df["IID"]

    for col in ["ea_years", "teacher_z", "ses_ea_proxy_z", "sex_c", "yob_c", "fractional_yob", "yob", "birth_year"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in task_columns.values():
        df[col] = pd.to_numeric(df[col], errors="coerce")

    yob_col = choose_yob_column(df)
    require(yob_col is not None, "no YOB-like column available after merging all_scores/base_covar/task scores")
    return df


def choose_yob_column(df: pd.DataFrame) -> str | None:
    for col in ("fractional_yob", "yob", "birth_year", "yob_c"):
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col
    return None


def finite_mask(*arrays: Iterable[float]) -> np.ndarray:
    masks = [np.isfinite(np.asarray(a, dtype=float)) for a in arrays]
    out = masks[0].copy()
    for mask in masks[1:]:
        out &= mask
    return out


def pearson_corr(x: Iterable[float], y: Iterable[float]) -> tuple[float, int]:
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    ok = finite_mask(xs, ys)
    n = int(ok.sum())
    if n < 3:
        return math.nan, n
    xs = xs[ok]
    ys = ys[ok]
    if np.nanstd(xs, ddof=1) <= 0 or np.nanstd(ys, ddof=1) <= 0:
        return math.nan, n
    return float(np.corrcoef(xs, ys)[0, 1]), n


def pearson_spearman(x: Iterable[float], y: Iterable[float]) -> tuple[float, float, int]:
    pearson, n = pearson_corr(x, y)
    xs = pd.Series(np.asarray(x, dtype=float))
    ys = pd.Series(np.asarray(y, dtype=float))
    ok = finite_mask(xs.to_numpy(), ys.to_numpy())
    if int(ok.sum()) < 3:
        return pearson, math.nan, n
    rx = xs.loc[ok].rank(method="average").to_numpy(dtype=float)
    ry = ys.loc[ok].rank(method="average").to_numpy(dtype=float)
    spearman, _ = pearson_corr(rx, ry)
    return pearson, spearman, n


def pattern_name(labels: tuple[str, ...], observed: np.ndarray) -> str:
    if not observed.any():
        return "NONE"
    if observed.all():
        return f"ALL{len(labels)}"
    if observed.sum() == 1:
        return f"{labels[int(np.flatnonzero(observed)[0])]}_ONLY"
    return "_".join(label for label, keep in zip(labels, observed) if keep)


def all_patterns(n_labels: int) -> list[np.ndarray]:
    return [np.asarray([(mask >> i) & 1 == 1 for i in range(n_labels)], dtype=bool) for mask in range(2**n_labels)]


def standardize_inputs(df: pd.DataFrame, labels: tuple[str, ...], columns: tuple[str, ...], min_n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    raw = df.loc[:, columns].to_numpy(dtype=float)
    cc = np.isfinite(raw).all(axis=1)
    n_cc = int(cc.sum())
    require(n_cc >= min_n, f"complete-case reference N={n_cc} is below minimum {min_n} for {labels}")
    means = np.nanmean(raw[cc], axis=0)
    sds = np.nanstd(raw[cc], axis=0, ddof=1)
    bad = (~np.isfinite(sds)) | (sds <= 0)
    require(not bool(bad.any()), f"invalid complete-case task SDs for {np.asarray(columns)[bad].tolist()}")
    standardized = (raw - means) / sds
    return standardized, cc, means, sds


def corr_long(labels: tuple[str, ...], values: np.ndarray, mask: np.ndarray, table_name: str) -> pd.DataFrame:
    rows = []
    for i, left in enumerate(labels):
        for j, right in enumerate(labels):
            if j <= i:
                continue
            ok = mask & np.isfinite(values[:, i]) & np.isfinite(values[:, j])
            pearson, spearman, n = pearson_spearman(values[:, i][ok], values[:, j][ok])
            rows.append(
                {
                    "table": table_name,
                    "task_x": left,
                    "task_y": right,
                    "n": n,
                    "pearson": pearson,
                    "spearman": spearman,
                }
            )
    return pd.DataFrame(rows)


def matrix_long(labels: tuple[str, ...], matrix: np.ndarray, value_name: str) -> pd.DataFrame:
    rows = []
    for i, row_label in enumerate(labels):
        for j, col_label in enumerate(labels):
            rows.append({"task_x": row_label, "task_y": col_label, value_name: float(matrix[i, j])})
    return pd.DataFrame(rows)


def safe_corrcoef(x: np.ndarray) -> np.ndarray:
    if x.shape[1] == 1:
        return np.ones((1, 1), dtype=float)
    return np.corrcoef(x, rowvar=False)


def model_implied_corr(loadings: np.ndarray, uniquenesses: np.ndarray) -> np.ndarray:
    cov = np.outer(loadings, loadings) + np.diag(uniquenesses)
    diag = np.sqrt(np.diag(cov))
    return cov / np.outer(diag, diag)


def fit_factor_and_score(
    df: pd.DataFrame,
    labels: tuple[str, ...],
    columns: tuple[str, ...],
    *,
    min_n: int,
    random_state: int,
    uniqueness_floor: float,
) -> FactorScoreResult:
    x_all, cc, means, sds = standardize_inputs(df, labels, columns, min_n)
    x_cc = x_all[cc]

    warning_messages: list[str] = []
    fa = FactorAnalysis(n_components=1, random_state=random_state)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        fa.fit(x_cc)
    for item in caught:
        warning_messages.append(str(item.message))

    loadings = fa.components_.T[:, 0].astype(float)
    uniquenesses = fa.noise_variance_.astype(float)
    require(np.isfinite(loadings).all(), f"non-finite FA loadings for {labels}")
    require(np.isfinite(uniquenesses).all(), f"non-finite FA uniquenesses for {labels}")

    orientation_sign = 1.0
    if float(loadings.sum()) < 0:
        loadings = -loadings
        orientation_sign = -1.0

    floored = uniquenesses < uniqueness_floor
    uniquenesses_for_scoring = np.maximum(uniquenesses, uniqueness_floor)
    if bool(floored.any()):
        warning_messages.append(
            "Uniqueness floor applied for scoring: "
            + ",".join(label for label, hit in zip(labels, floored) if hit)
        )

    pca = PCA(n_components=1, random_state=random_state)
    pca.fit(x_cc)
    pca_loadings = pca.components_[0].astype(float)
    if float(pca_loadings.sum()) < 0:
        pca_loadings = -pca_loadings

    weights_by_pattern, pattern_to_weights = scoring_weights(labels, loadings, uniquenesses_for_scoring, x_all)
    g_hat = np.full(x_all.shape[0], np.nan, dtype=float)
    task_patterns = observed_patterns(labels, x_all)
    for pattern, weights in pattern_to_weights.items():
        if pattern == "NONE":
            continue
        row_mask = task_patterns == pattern
        if not bool(row_mask.any()):
            continue
        g_hat[row_mask] = np.nansum(x_all[row_mask] * weights.reshape(1, -1), axis=1)

    mu_g = float(np.nanmean(g_hat[cc]))
    sd_g = float(np.nanstd(g_hat[cc], ddof=1))
    require(np.isfinite(sd_g) and sd_g > 0, f"invalid complete-case g_hat SD for {labels}")
    g_z = (g_hat - mu_g) / sd_g

    obs_corr = safe_corrcoef(x_cc)
    mod_corr = model_implied_corr(loadings, uniquenesses)
    residual_corr = obs_corr - mod_corr
    fa_loglike_last = None
    if getattr(fa, "loglike_", None):
        fa_loglike_last = float(fa.loglike_[-1])

    return FactorScoreResult(
        labels=labels,
        columns=columns,
        means=means,
        sds=sds,
        complete_case_mask=cc,
        standardized=x_all,
        loadings=loadings,
        uniquenesses=uniquenesses,
        uniquenesses_for_scoring=uniquenesses_for_scoring,
        floored_uniqueness=floored,
        orientation_sign=orientation_sign,
        weights_by_pattern=weights_by_pattern,
        g_hat=g_hat,
        g_z_candidate=g_z,
        g_hat_mean_complete_case=mu_g,
        g_hat_sd_complete_case=sd_g,
        observed_corr_complete_case=obs_corr,
        model_corr=mod_corr,
        residual_corr=residual_corr,
        pca_loadings=pca_loadings,
        pca_explained_variance_ratio=float(pca.explained_variance_ratio_[0]),
        fa_n_iter=int(getattr(fa, "n_iter_", -1)) if getattr(fa, "n_iter_", None) is not None else None,
        fa_loglike_last=fa_loglike_last,
        warning_messages=tuple(warning_messages),
    )


def observed_patterns(labels: tuple[str, ...], x_all: np.ndarray) -> np.ndarray:
    return np.asarray([pattern_name(labels, np.isfinite(row)) for row in x_all], dtype=object)


def scoring_weights(
    labels: tuple[str, ...], loadings: np.ndarray, uniquenesses_for_scoring: np.ndarray, x_all: np.ndarray
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    task_patterns = observed_patterns(labels, x_all)
    rows = []
    lookup: dict[str, np.ndarray] = {}
    for observed in all_patterns(len(labels)):
        pattern = pattern_name(labels, observed)
        weights_full = np.zeros(len(labels), dtype=float)
        used_pinv = 0
        if observed.any():
            lam = loadings[observed]
            psi = uniquenesses_for_scoring[observed]
            sigma = np.outer(lam, lam) + np.diag(psi)
            try:
                weights_obs = np.linalg.solve(sigma, lam)
            except np.linalg.LinAlgError:
                weights_obs = np.linalg.pinv(sigma) @ lam
                used_pinv = 1
            weights_full[observed] = weights_obs
        lookup[pattern] = weights_full
        row = {
            "task_pattern": pattern,
            "n_tasks_observed": int(observed.sum()),
            "n_samples": int((task_patterns == pattern).sum()),
            "used_pinv": used_pinv,
        }
        for label, weight in zip(labels, weights_full):
            row[f"weight_{label}"] = float(weight)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_tasks_observed", "task_pattern"]), lookup


def accepted_domain_factor(
    result: FactorScoreResult,
    complete_case_corr: pd.DataFrame,
    *,
    min_loadings_ge20: int,
) -> tuple[bool, str]:
    positive_loadings = bool((result.loadings > 0).all())
    enough_loadings = int((result.loadings >= 0.20).sum()) >= min_loadings_ge20
    min_pair_corr = float(complete_case_corr["pearson"].min()) if not complete_case_corr.empty else math.nan
    severe_contradiction = np.isfinite(min_pair_corr) and min_pair_corr < -0.05
    accepted = positive_loadings and enough_loadings and not severe_contradiction
    reasons = []
    if not positive_loadings:
        reasons.append("nonpositive_loading")
    if not enough_loadings:
        reasons.append(f"fewer_than_{min_loadings_ge20}_loadings_ge_0.20")
    if severe_contradiction:
        reasons.append("complete_case_pair_corr_lt_-0.05")
    if not reasons:
        reasons.append("accepted")
    return accepted, ",".join(reasons)


def accepted_three_domain(result: FactorScoreResult, complete_case_corr: pd.DataFrame) -> tuple[bool, str]:
    return accepted_domain_factor(result, complete_case_corr, min_loadings_ge20=2)


def accepted_four_domain(result: FactorScoreResult, complete_case_corr: pd.DataFrame) -> tuple[bool, str]:
    return accepted_domain_factor(result, complete_case_corr, min_loadings_ge20=3)


def add_observed_task_columns(df: pd.DataFrame, x_all: np.ndarray, labels: tuple[str, ...]) -> pd.DataFrame:
    out = df.copy()
    obs = np.isfinite(x_all)
    for idx, label in enumerate(labels):
        out[f"has_{label.lower()}"] = obs[:, idx].astype(int)
    out["n_tasks_observed"] = obs.sum(axis=1).astype(int)
    out["task_pattern"] = observed_patterns(labels, x_all)
    return out


def quantile_or_nan(values: pd.Series, q: float) -> float:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    if vals.empty:
        return math.nan
    return float(vals.quantile(q))


def score_distributions(scored: pd.DataFrame, score_col: str) -> pd.DataFrame:
    group_defs = [("overall", None), ("n_tasks_observed", "n_tasks_observed"), ("task_pattern", "task_pattern")]
    rows = []
    for group_type, col in group_defs:
        if col is None:
            groups = [("all", scored)]
        else:
            groups = list(scored.groupby(col, dropna=False))
        for group_value, sub in groups:
            vals = pd.to_numeric(sub[score_col], errors="coerce").dropna()
            rows.append(
                {
                    "group_type": group_type,
                    "group_value": str(group_value),
                    "n": int(vals.shape[0]),
                    "mean": float(vals.mean()) if not vals.empty else math.nan,
                    "sd": float(vals.std(ddof=1)) if vals.shape[0] >= 2 else math.nan,
                    "min": float(vals.min()) if not vals.empty else math.nan,
                    "p01": quantile_or_nan(vals, 0.01),
                    "p05": quantile_or_nan(vals, 0.05),
                    "p50": quantile_or_nan(vals, 0.50),
                    "p95": quantile_or_nan(vals, 0.95),
                    "p99": quantile_or_nan(vals, 0.99),
                    "max": float(vals.max()) if not vals.empty else math.nan,
                }
            )
    return pd.DataFrame(rows)


def validation_groups(df: pd.DataFrame, score_col: str | None = None) -> list[tuple[str, str, pd.DataFrame]]:
    if score_col is not None and score_col in df.columns:
        df = df.loc[pd.to_numeric(df[score_col], errors="coerce").notna()].copy()
    groups: list[tuple[str, str, pd.DataFrame]] = [("overall", "all", df)]
    for col in ("role", "n_tasks_observed", "task_pattern"):
        if col not in df.columns:
            continue
        for value, sub in df.groupby(col, dropna=False):
            if len(sub) >= MIN_GROUP_N:
                groups.append((col, str(value), sub))
    return groups


def regression_fit(y: pd.Series, design: pd.DataFrame) -> dict[str, object]:
    yv = pd.to_numeric(y, errors="coerce").to_numpy(dtype=float)
    x = design.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(yv) & np.isfinite(x).all(axis=1)
    n = int(ok.sum())
    p = x.shape[1]
    if n <= p:
        return {"n": n, "r2": math.nan, "rank": math.nan, "coefficients": {}}
    beta, _, rank, _ = np.linalg.lstsq(x[ok], yv[ok], rcond=None)
    pred = x[ok] @ beta
    ss_res = float(np.sum((yv[ok] - pred) ** 2))
    ss_tot = float(np.sum((yv[ok] - np.mean(yv[ok])) ** 2))
    r2 = math.nan if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    return {
        "n": n,
        "r2": float(r2),
        "rank": int(rank),
        "coefficients": {name: float(val) for name, val in zip(design.columns, beta)},
    }


def age_sex_validation(scored: pd.DataFrame, score_col: str, yob_col: str) -> pd.DataFrame:
    rows = []
    for group_type, group_value, sub in validation_groups(scored, score_col):
        y = pd.to_numeric(sub[score_col], errors="coerce")
        sex = pd.to_numeric(sub["sex_c"], errors="coerce") if "sex_c" in sub.columns else pd.Series(np.nan, index=sub.index)
        yob = pd.to_numeric(sub[yob_col], errors="coerce")
        p_sex, s_sex, n_sex = pearson_spearman(y, sex)
        p_yob, s_yob, n_yob = pearson_spearman(y, yob)
        design_linear = pd.DataFrame({"intercept": 1.0, "sex_c": sex, "yob_like": yob}, index=sub.index)
        design_quad = design_linear.copy()
        design_quad["yob_like_sq"] = yob**2
        fit_linear = regression_fit(y, design_linear)
        fit_quad = regression_fit(y, design_quad)
        row = {
            "group_type": group_type,
            "group_value": group_value,
            "score_column": score_col,
            "yob_column": yob_col,
            "n": int(y.notna().sum()),
            "pearson_sex_c": p_sex,
            "spearman_sex_c": s_sex,
            "n_sex_c": n_sex,
            "pearson_yob_like": p_yob,
            "spearman_yob_like": s_yob,
            "n_yob_like": n_yob,
            "linear_n": fit_linear["n"],
            "linear_r2": fit_linear["r2"],
            "linear_rank": fit_linear["rank"],
            "quad_n": fit_quad["n"],
            "quad_r2": fit_quad["r2"],
            "quad_rank": fit_quad["rank"],
        }
        for name, value in fit_linear["coefficients"].items():
            row[f"linear_beta_{name}"] = value
        for name, value in fit_quad["coefficients"].items():
            row[f"quad_beta_{name}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


def external_validation(scored: pd.DataFrame, score_col: str) -> pd.DataFrame:
    rows = []
    for group_type, group_value, sub in validation_groups(scored, score_col):
        for target in ("ses_ea_proxy_z", "teacher_z", "ea_years"):
            if target not in sub.columns:
                continue
            pearson, spearman, n = pearson_spearman(sub[score_col], sub[target])
            rows.append(
                {
                    "group_type": group_type,
                    "group_value": group_value,
                    "score_column": score_col,
                    "target": target,
                    "n": n,
                    "pearson": pearson,
                    "spearman": spearman,
                }
            )
    return pd.DataFrame(rows)


def z_by_reference(values: np.ndarray, ref_mask: np.ndarray) -> np.ndarray:
    ref = values[ref_mask & np.isfinite(values)]
    if ref.shape[0] < 3:
        return np.full(values.shape[0], np.nan)
    sd = float(np.nanstd(ref, ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        return np.full(values.shape[0], np.nan)
    return (values - float(np.nanmean(ref))) / sd


def comparison_scores(df: pd.DataFrame, main: FactorScoreResult, attention_z: np.ndarray, primary_score_col: str) -> pd.DataFrame:
    x = main.standardized
    cc = main.complete_case_mask
    obs = np.isfinite(x)
    pca_complete = np.full(x.shape[0], np.nan, dtype=float)
    pca_complete[cc] = x[cc] @ main.pca_loadings
    available_mean = np.full(x.shape[0], np.nan, dtype=float)
    has_any = obs.any(axis=1)
    available_mean[has_any] = np.nanmean(x[has_any], axis=1)
    pca_weighted = np.full(x.shape[0], np.nan, dtype=float)
    for idx, row_obs in enumerate(obs):
        if not row_obs.any():
            continue
        denom = float(main.pca_loadings[row_obs].sum())
        if denom == 0 or not np.isfinite(denom):
            continue
        pca_weighted[idx] = float(np.nansum(x[idx, row_obs] * main.pca_loadings[row_obs]) / denom)

    candidates = {
        "complete_case_pca_pc1": z_by_reference(pca_complete, cc),
        "available_task_unweighted_mean": z_by_reference(available_mean, cc),
        "available_task_pca_weighted_mean_div_available_weights": z_by_reference(pca_weighted, cc),
        "gradcpt_flanker_attention_exec": attention_z,
    }

    rows = []
    primary = pd.to_numeric(df[primary_score_col], errors="coerce")
    for name, values in candidates.items():
        row = {"comparison_score": name}
        row["n"] = int(np.isfinite(values).sum())
        row["mean"] = float(np.nanmean(values)) if np.isfinite(values).any() else math.nan
        row["sd"] = float(np.nanstd(values, ddof=1)) if np.isfinite(values).sum() >= 2 else math.nan
        p, s, n = pearson_spearman(values, primary)
        row.update({"target": primary_score_col, "target_n": n, "pearson": p, "spearman": s})
        rows.append(row.copy())
        for target in ("ses_ea_proxy_z", "teacher_z", "ea_years"):
            p, s, n = pearson_spearman(values, df[target])
            row2 = {
                "comparison_score": name,
                "n": int(np.isfinite(values).sum()),
                "mean": row["mean"],
                "sd": row["sd"],
                "target": target,
                "target_n": n,
                "pearson": p,
                "spearman": s,
            }
            rows.append(row2)
    return pd.DataFrame(rows)


def main_score_tables(
    df: pd.DataFrame,
    task_columns: dict[str, str],
    flanker_source: str,
    main: FactorScoreResult,
    attention: FactorScoreResult,
    accepted: bool,
    acceptance_reason: str,
    force_three_domain_g: bool,
    four: FactorScoreResult | None = None,
    accepted_four: bool | None = None,
    acceptance_reason_four: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scored = add_observed_task_columns(df, main.standardized, main.labels)
    scored["etm_g_hat"] = main.g_hat
    scored["accepted_three_domain_g"] = int(accepted)
    scored["etm_g_z"] = main.g_z_candidate if accepted else np.nan
    if force_three_domain_g and not accepted:
        scored["etm_g_z_forced"] = main.g_z_candidate
    scored["flanker_input_source"] = flanker_source
    scored["etm_attention_exec_z"] = attention.g_z_candidate
    scored["etm_g_acceptance_reason"] = acceptance_reason

    if four is not None:
        four_obs = np.isfinite(four.standardized)
        scored["etm_g4_hat"] = four.g_hat
        scored["accepted_four_domain_g"] = int(bool(accepted_four))
        scored["etm_g4_z"] = four.g_z_candidate if accepted_four else np.nan
        scored["n_tasks_observed_four_domain"] = four_obs.sum(axis=1).astype(int)
        scored["task_pattern_four_domain"] = observed_patterns(four.labels, four.standardized)
        for idx, label in enumerate(four.labels):
            scored[f"has_four_{label.lower()}"] = four_obs[:, idx].astype(int)
        scored["etm_g4_acceptance_reason"] = acceptance_reason_four or ""

    base_cols = [
        "FID",
        "IID",
        "person_id",
        "role",
        "fold_id",
        "etm_g_z",
        "etm_g_hat",
        "accepted_three_domain_g",
        "n_tasks_observed",
        "task_pattern",
        "has_dd",
        "has_gradcpt",
        "has_flanker",
    ]
    if "etm_g_z_forced" in scored.columns:
        base_cols.append("etm_g_z_forced")
    if four is not None:
        base_cols.extend(
            [
                "etm_g4_z",
                "etm_g4_hat",
                "accepted_four_domain_g",
                "n_tasks_observed_four_domain",
                "task_pattern_four_domain",
                "has_four_dd",
                "has_four_gradcpt",
                "has_four_flanker",
                "has_four_emorecog",
            ]
        )
    extra_cols = [
        task_columns["DD"],
        task_columns["GRADCPT"],
        task_columns["FLANKER"],
        task_columns.get("EMORECOG"),
        "flanker_input_source",
        "etm_attention_exec_z",
        "ses_ea_proxy_z",
        "teacher_z",
        "ea_years",
        "sex_c",
        "yob_c",
        "fractional_yob",
        "yob",
        "birth_year",
    ]
    ordered = base_cols + [col for col in extra_cols if col and col in scored.columns and col not in base_cols]
    wide = scored.loc[:, ordered].copy()
    if "n_tasks_observed_four_domain" in wide.columns:
        scored_only = wide.loc[(wide["n_tasks_observed"] >= 1) | (wide["n_tasks_observed_four_domain"] >= 1)].copy()
    else:
        scored_only = wide.loc[wide["n_tasks_observed"] >= 1].copy()
    return wide, scored_only


def reference_standardization_table(result: FactorScoreResult) -> pd.DataFrame:
    rows = []
    n_cc = int(result.complete_case_mask.sum())
    for label, col, mean, sd in zip(result.labels, result.columns, result.means, result.sds):
        rows.append(
            {
                "task": label,
                "column": col,
                "complete_case_n": n_cc,
                "mean_complete_case": float(mean),
                "sd_complete_case_ddof1": float(sd),
            }
        )
    return pd.DataFrame(rows)


def task_missingness_counts(df: pd.DataFrame, task_columns: dict[str, str], labels: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for label in labels:
        col = task_columns[label]
        nonmissing = int(pd.to_numeric(df[col], errors="coerce").notna().sum())
        rows.append(
            {
                "task": label,
                "column": col,
                "n_total": int(len(df)),
                "n_nonmissing": nonmissing,
                "n_missing": int(len(df) - nonmissing),
            }
        )
    return pd.DataFrame(rows)


def pattern_counts(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append({"group_type": "overall", "group_value": "all", "n": int(len(scored))})
    for col in ("role", "n_tasks_observed", "task_pattern", "n_tasks_observed_four_domain", "task_pattern_four_domain"):
        if col not in scored.columns:
            continue
        for value, sub in scored.groupby(col, dropna=False):
            rows.append({"group_type": col, "group_value": str(value), "n": int(len(sub))})
    return pd.DataFrame(rows)


def loadings_table(result: FactorScoreResult, accepted: bool) -> pd.DataFrame:
    rows = []
    for label, col, loading, pca_loading in zip(result.labels, result.columns, result.loadings, result.pca_loadings):
        rows.append(
            {
                "task": label,
                "column": col,
                "loading": float(loading),
                "abs_loading": float(abs(loading)),
                "loading_positive": int(loading > 0),
                "loading_ge_0.20": int(loading >= 0.20),
                "pca_loading": float(pca_loading),
                "accepted_factor": int(accepted),
            }
        )
    return pd.DataFrame(rows)


def uniqueness_table(result: FactorScoreResult) -> pd.DataFrame:
    rows = []
    for label, col, psi, psi_score, floored in zip(
        result.labels, result.columns, result.uniquenesses, result.uniquenesses_for_scoring, result.floored_uniqueness
    ):
        rows.append(
            {
                "task": label,
                "column": col,
                "uniqueness": float(psi),
                "uniqueness_for_scoring": float(psi_score),
                "uniqueness_floored": int(floored),
            }
        )
    return pd.DataFrame(rows)


def factor_summary_table(result: FactorScoreResult, accepted: bool, reason: str, random_state: int) -> pd.DataFrame:
    rows = [
        {"parameter": "complete_case_n", "value": int(result.complete_case_mask.sum())},
        {"parameter": "fa_random_state", "value": random_state},
        {"parameter": "fa_n_iter", "value": result.fa_n_iter},
        {"parameter": "fa_loglike_last", "value": result.fa_loglike_last},
        {"parameter": "orientation_sign", "value": result.orientation_sign},
        {"parameter": "accepted_factor", "value": int(accepted)},
        {"parameter": "acceptance_reason", "value": reason},
        {"parameter": "pca_pc1_explained_variance_ratio", "value": result.pca_explained_variance_ratio},
        {"parameter": "complete_case_g_hat_mean", "value": result.g_hat_mean_complete_case},
        {"parameter": "complete_case_g_hat_sd_ddof1", "value": result.g_hat_sd_complete_case},
    ]
    for idx, message in enumerate(result.warning_messages, start=1):
        rows.append({"parameter": f"warning_{idx}", "value": message})
    return pd.DataFrame(rows)


def attention_summary(attention: FactorScoreResult, scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, col, loading, psi in zip(attention.labels, attention.columns, attention.loadings, attention.uniquenesses):
        rows.append(
            {
                "section": "loading",
                "metric": label,
                "column": col,
                "value": float(loading),
                "extra": f"uniqueness={float(psi):.8g}",
            }
        )
    for target in ("ses_ea_proxy_z", "teacher_z", "ea_years"):
        p, s, n = pearson_spearman(scored["etm_attention_exec_z"], scored[target])
        rows.append({"section": "external_validation", "metric": target, "column": "etm_attention_exec_z", "value": p, "extra": f"n={n};spearman={s}"})
    rows.append({"section": "reference", "metric": "complete_case_n", "column": "", "value": int(attention.complete_case_mask.sum()), "extra": ""})
    return pd.DataFrame(rows)


def compare_three_four_models(
    scored: pd.DataFrame,
    three: FactorScoreResult,
    four: FactorScoreResult | None,
    accepted_three: bool,
    accepted_four: bool | None,
) -> pd.DataFrame:
    rows = []
    model_defs = [
        ("three_domain", "etm_g_z", three, accepted_three),
    ]
    if four is not None:
        model_defs.append(("four_domain", "etm_g4_z", four, bool(accepted_four)))

    for model_name, score_col, result, accepted in model_defs:
        rows.append(
            {
                "section": "reference",
                "model": model_name,
                "metric": "complete_case_n",
                "target": "",
                "task": "",
                "value": int(result.complete_case_mask.sum()),
                "n": int(result.complete_case_mask.sum()),
                "accepted": int(accepted),
                "spearman": math.nan,
            }
        )
        rows.append(
            {
                "section": "reference",
                "model": model_name,
                "metric": "pc1_explained_variance_ratio",
                "target": "",
                "task": "",
                "value": result.pca_explained_variance_ratio,
                "n": int(result.complete_case_mask.sum()),
                "accepted": int(accepted),
                "spearman": math.nan,
            }
        )
        for label, col, loading, psi in zip(result.labels, result.columns, result.loadings, result.uniquenesses):
            rows.append(
                {
                    "section": "loading",
                    "model": model_name,
                    "metric": "fa_loading",
                    "target": "",
                    "task": label,
                    "column": col,
                    "value": float(loading),
                    "n": int(result.complete_case_mask.sum()),
                    "accepted": int(accepted),
                    "spearman": math.nan,
                    "extra": f"uniqueness={float(psi):.8g}",
                }
            )
        for target in ("ses_ea_proxy_z", "teacher_z", "ea_years"):
            p, s, n = pearson_spearman(scored[score_col], scored[target])
            rows.append(
                {
                    "section": "external_validation",
                    "model": model_name,
                    "metric": "pearson",
                    "target": target,
                    "task": "",
                    "column": score_col,
                    "value": p,
                    "n": n,
                    "accepted": int(accepted),
                    "spearman": s,
                }
            )

    if four is not None and "etm_g4_z" in scored.columns:
        overlap = np.isfinite(pd.to_numeric(scored["etm_g_z"], errors="coerce")) & np.isfinite(pd.to_numeric(scored["etm_g4_z"], errors="coerce"))
        p, s, n = pearson_spearman(scored.loc[overlap, "etm_g_z"], scored.loc[overlap, "etm_g4_z"])
        rows.append(
            {
                "section": "model_overlap",
                "model": "three_vs_four",
                "metric": "score_correlation",
                "target": "etm_g_z_vs_etm_g4_z",
                "task": "",
                "column": "etm_g_z,etm_g4_z",
                "value": p,
                "n": n,
                "accepted": int(bool(accepted_three) and bool(accepted_four)),
                "spearman": s,
            }
        )
        for target in ("ses_ea_proxy_z", "teacher_z", "ea_years"):
            for score_col, model_name in (("etm_g_z", "three_domain_overlap"), ("etm_g4_z", "four_domain_overlap")):
                p, s, n = pearson_spearman(scored.loc[overlap, score_col], scored.loc[overlap, target])
                rows.append(
                    {
                        "section": "external_validation_overlap",
                        "model": model_name,
                        "metric": "pearson",
                        "target": target,
                        "task": "",
                        "column": score_col,
                        "value": p,
                        "n": n,
                        "accepted": int(score_col == "etm_g_z" and accepted_three or score_col == "etm_g4_z" and bool(accepted_four)),
                        "spearman": s,
                    }
                )
    return pd.DataFrame(rows)


def write_outputs(
    args: argparse.Namespace,
    wide: pd.DataFrame,
    scored_only: pd.DataFrame,
    diagnostics: dict[str, pd.DataFrame],
) -> None:
    diag_dir = args.output_dir / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    wide_path = args.output_dir / "etm_general_factor_scores_wide.tsv"
    scored_path = args.output_dir / "etm_general_factor_scores_scored_only.tsv"
    wide.to_csv(wide_path, sep="\t", index=False)
    scored_only.to_csv(scored_path, sep="\t", index=False)
    print(f"Wrote {wide_path}", flush=True)
    print(f"Wrote {scored_path}", flush=True)

    for name, table in diagnostics.items():
        path = diag_dir / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        print(f"Wrote {path}", flush=True)

    if args.stage_aggregate:
        require(args.workspace_scrap_dir is not None, "--stage-aggregate requires --workspace-scrap-dir")
        args.workspace_scrap_dir.mkdir(parents=True, exist_ok=True)
        for path in diag_dir.glob("*.tsv"):
            dest = args.workspace_scrap_dir / path.name
            shutil.copyfile(path, dest)
            print(f"Staged aggregate diagnostic {dest}", flush=True)


def complete_case_corr_tables(main: FactorScoreResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    complete = corr_long(main.labels, main.standardized, main.complete_case_mask, "complete_case")
    pair = corr_long(main.labels, main.standardized, np.ones(main.standardized.shape[0], dtype=bool), "pair_available")
    return complete, pair


def main() -> None:
    args = parse_args()
    require(args.uniqueness_floor > 0, "--uniqueness-floor must be positive")
    require(args.min_complete_case_n >= 3, "--min-complete-case-n must be at least 3")
    if args.base_covar_file is None:
        args.base_covar_file = args.all_scores_file.parent / "base_covar.txt"

    task_header = read_header(args.task_score_file)
    flanker_col, flanker_source = resolve_flanker_column(task_header, args.flanker_input)
    task_columns = {label: DEFAULT_TASK_COLUMNS[label] for label in THREE_TASK_LABELS}
    task_columns["FLANKER"] = flanker_col
    four_domain_enabled = DEFAULT_TASK_COLUMNS["EMORECOG"] in task_header
    if four_domain_enabled:
        task_columns["EMORECOG"] = DEFAULT_TASK_COLUMNS["EMORECOG"]
    missing = [col for col in task_columns.values() if col not in task_header]
    require(not missing, f"task-score file missing selected task columns: {missing}")

    signature, signature_rows = build_run_signature(args, task_columns, flanker_source)
    check_idempotency(args, signature)

    print("=== ETM general factor scoring ===", flush=True)
    print(f"Task score file: {args.task_score_file}", flush=True)
    print(f"All scores file: {args.all_scores_file}", flush=True)
    print(f"Base covar file: {args.base_covar_file}", flush=True)
    print(f"Output dir: {args.output_dir}", flush=True)
    print(f"Selected task columns: {task_columns}", flush=True)
    print(f"Flanker input source: {flanker_source}", flush=True)
    print(f"Four-domain ETM-g enabled: {int(four_domain_enabled)}", flush=True)

    df = load_inputs(args, task_columns)
    print(f"SES-EA proxy cohort rows: {len(df)}", flush=True)
    print(df["role"].value_counts(dropna=False).to_string(), flush=True)

    main_result = fit_factor_and_score(
        df,
        THREE_TASK_LABELS,
        tuple(task_columns[label] for label in THREE_TASK_LABELS),
        min_n=args.min_complete_case_n,
        random_state=args.random_state,
        uniqueness_floor=args.uniqueness_floor,
    )
    complete_corr, pair_corr = complete_case_corr_tables(main_result)
    accepted, acceptance_reason = accepted_three_domain(main_result, complete_corr)
    print(f"Complete-case reference N: {int(main_result.complete_case_mask.sum())}", flush=True)
    print(f"FA loadings: {dict(zip(TASK_LABELS, main_result.loadings))}", flush=True)
    print(f"Accepted three-domain ETM-g: {int(accepted)} ({acceptance_reason})", flush=True)

    four_result = None
    complete_corr_four = None
    pair_corr_four = None
    accepted4 = None
    acceptance_reason4 = None
    if four_domain_enabled:
        four_result = fit_factor_and_score(
            df,
            FOUR_TASK_LABELS,
            tuple(task_columns[label] for label in FOUR_TASK_LABELS),
            min_n=args.min_complete_case_n,
            random_state=args.random_state,
            uniqueness_floor=args.uniqueness_floor,
        )
        complete_corr_four, pair_corr_four = complete_case_corr_tables(four_result)
        accepted4, acceptance_reason4 = accepted_four_domain(four_result, complete_corr_four)
        print(f"Four-domain complete-case reference N: {int(four_result.complete_case_mask.sum())}", flush=True)
        print(f"Four-domain FA loadings: {dict(zip(FOUR_TASK_LABELS, four_result.loadings))}", flush=True)
        print(f"Accepted four-domain ETM-g: {int(accepted4)} ({acceptance_reason4})", flush=True)

    attention_result = fit_factor_and_score(
        df,
        ("GRADCPT", "FLANKER"),
        (task_columns["GRADCPT"], task_columns["FLANKER"]),
        min_n=args.min_complete_case_n,
        random_state=args.random_state,
        uniqueness_floor=args.uniqueness_floor,
    )

    wide, scored_only = main_score_tables(
        df,
        task_columns,
        flanker_source,
        main_result,
        attention_result,
        accepted,
        acceptance_reason,
        args.force_three_domain_g,
        four_result,
        accepted4,
        acceptance_reason4,
    )
    diagnostic_df = wide.copy()
    if accepted:
        diagnostic_score_col = "etm_g_z"
    elif "etm_g_z_forced" in diagnostic_df.columns:
        diagnostic_score_col = "etm_g_z_forced"
    else:
        diagnostic_score_col = "etm_g_z_candidate_unaccepted"
        diagnostic_df[diagnostic_score_col] = main_result.g_z_candidate

    yob_col = choose_yob_column(wide)
    require(yob_col is not None, "no YOB-like column available for validation")

    params = pd.concat(
        [
            signature_rows,
            pd.DataFrame(
                [
                    {"parameter": "python_version", "value": platform.python_version()},
                    {"parameter": "pandas_version", "value": pd.__version__},
                    {"parameter": "numpy_version", "value": np.__version__},
                    {"parameter": "scipy_version", "value": scipy.__version__},
                    {"parameter": "sklearn_version", "value": sklearn_version},
                    {"parameter": "three_domain_task_column_order", "value": ",".join(task_columns[label] for label in THREE_TASK_LABELS)},
                    {
                        "parameter": "four_domain_task_column_order",
                        "value": ",".join(task_columns[label] for label in FOUR_TASK_LABELS if label in task_columns),
                    },
                    {"parameter": "four_domain_enabled", "value": int(four_domain_enabled)},
                    {"parameter": "flanker_input_source", "value": flanker_source},
                    {"parameter": "complete_case_reference_n", "value": int(main_result.complete_case_mask.sum())},
                    {"parameter": "complete_case_g_hat_mean", "value": main_result.g_hat_mean_complete_case},
                    {"parameter": "complete_case_g_hat_sd_ddof1", "value": main_result.g_hat_sd_complete_case},
                    {"parameter": "accepted_three_domain_g", "value": int(accepted)},
                    {"parameter": "acceptance_reason", "value": acceptance_reason},
                    {"parameter": "four_domain_complete_case_reference_n", "value": int(four_result.complete_case_mask.sum()) if four_result is not None else ""},
                    {"parameter": "accepted_four_domain_g", "value": int(accepted4) if accepted4 is not None else ""},
                    {"parameter": "four_domain_acceptance_reason", "value": acceptance_reason4 or ""},
                    {"parameter": "uniqueness_floor", "value": args.uniqueness_floor},
                    {"parameter": "orientation_rule", "value": "if sum(lambda_vec) < 0, multiply loadings by -1"},
                    {"parameter": "orientation_sign", "value": main_result.orientation_sign},
                    {"parameter": "yob_validation_column", "value": yob_col},
                ]
            ),
        ],
        ignore_index=True,
    )

    diagnostics = {
        "etm_g_reference_standardization": reference_standardization_table(main_result),
        "etm_g_task_missingness_counts": task_missingness_counts(df, task_columns, tuple(label for label in FOUR_TASK_LABELS if label in task_columns)),
        "etm_g_task_pattern_counts": pattern_counts(diagnostic_df),
        "etm_g_complete_case_task_correlations": complete_corr,
        "etm_g_pair_available_task_correlations": pair_corr,
        "etm_g_fa_loadings": loadings_table(main_result, accepted),
        "etm_g_uniquenesses": uniqueness_table(main_result),
        "etm_g_model_implied_correlations": matrix_long(TASK_LABELS, main_result.model_corr, "model_implied_correlation"),
        "etm_g_residual_correlations": matrix_long(TASK_LABELS, main_result.residual_corr, "residual_correlation"),
        "etm_g_scoring_weights_by_pattern": main_result.weights_by_pattern,
        "etm_g_score_distribution_by_pattern": score_distributions(diagnostic_df, diagnostic_score_col),
        "etm_g_age_sex_validation": age_sex_validation(diagnostic_df, diagnostic_score_col, yob_col),
        "etm_g_external_validation_correlations": external_validation(diagnostic_df, diagnostic_score_col),
        "etm_g_comparison_scores_summary": comparison_scores(diagnostic_df, main_result, attention_result.g_z_candidate, diagnostic_score_col),
        "etm_attention_exec_diagnostic_summary": attention_summary(attention_result, diagnostic_df),
        "etm_g_factor_summary": factor_summary_table(main_result, accepted, acceptance_reason, args.random_state),
        "etm_g_three_vs_four_comparison": compare_three_four_models(diagnostic_df, main_result, four_result, accepted, accepted4),
        "etm_g_reproducibility_params": params,
    }
    if four_result is not None and complete_corr_four is not None and pair_corr_four is not None:
        diagnostics.update(
            {
                "etm_g4_reference_standardization": reference_standardization_table(four_result),
                "etm_g4_complete_case_task_correlations": complete_corr_four,
                "etm_g4_pair_available_task_correlations": pair_corr_four,
                "etm_g4_fa_loadings": loadings_table(four_result, bool(accepted4)),
                "etm_g4_uniquenesses": uniqueness_table(four_result),
                "etm_g4_model_implied_correlations": matrix_long(FOUR_TASK_LABELS, four_result.model_corr, "model_implied_correlation"),
                "etm_g4_residual_correlations": matrix_long(FOUR_TASK_LABELS, four_result.residual_corr, "residual_correlation"),
                "etm_g4_scoring_weights_by_pattern": four_result.weights_by_pattern,
                "etm_g4_score_distribution_by_pattern": score_distributions(diagnostic_df, "etm_g4_z"),
                "etm_g4_age_sex_validation": age_sex_validation(diagnostic_df, "etm_g4_z", yob_col),
                "etm_g4_external_validation_correlations": external_validation(diagnostic_df, "etm_g4_z"),
                "etm_g4_factor_summary": factor_summary_table(four_result, bool(accepted4), acceptance_reason4 or "", args.random_state),
            }
        )
    write_outputs(args, wide, scored_only, diagnostics)

    print("\nETM-g external validation correlations:", flush=True)
    ext = diagnostics["etm_g_external_validation_correlations"]
    print(ext.loc[ext["group_type"] == "overall"].to_string(index=False), flush=True)
    if "etm_g4_external_validation_correlations" in diagnostics:
        print("\nFour-domain ETM-g external validation correlations:", flush=True)
        ext4 = diagnostics["etm_g4_external_validation_correlations"]
        print(ext4.loc[ext4["group_type"] == "overall"].to_string(index=False), flush=True)
        print("\nThree-vs-four ETM-g comparison:", flush=True)
        comp = diagnostics["etm_g_three_vs_four_comparison"]
        keep = comp["section"].isin(["loading", "external_validation", "external_validation_overlap", "model_overlap"])
        print(comp.loc[keep].to_string(index=False), flush=True)
    print("\n=== ETM general factor scoring complete ===", flush=True)


if __name__ == "__main__":
    main()
