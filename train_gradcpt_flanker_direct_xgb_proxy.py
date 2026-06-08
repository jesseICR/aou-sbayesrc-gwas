#!/usr/bin/env python3
"""Train a direct survey-to-GradCPT/Flanker XGBoost proxy.

This downstream helper keeps the SES-EA proxy fold structure, adds the
education-response item as a feature, trains scratch XGBoost models on a
missing-pattern-aware GradCPT/Flanker target, and writes the direct survey
proxy needed by the final no-teacher calibration wrapper.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import shlex
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import setup_ses_ea_proxy_gwas as ses


EDU_QID = 1585940
EDU_ANSWER_IDS = [1585941, 1585942, 1585943, 1585944, 1585945, 1585946, 1585947, 1585948]
MODEL_NAMES = [f"fold_{i}" for i in range(5)] + ["final_model"]
MIN_GROUP_N = 20


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def env_int(name: str, default: int) -> int:
    val = os.environ.get(name)
    return int(val) if val not in (None, "") else default


def env_float(name: str, default: float) -> float:
    val = os.environ.get(name)
    return float(val) if val not in (None, "") else default


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ses-ea-dir", type=Path, required=True)
    parser.add_argument("--all-scores-file", type=Path, required=True)
    parser.add_argument("--task-score-file", type=Path, required=True)
    parser.add_argument("--fine-tuned-score-file", type=Path, required=True)
    parser.add_argument("--ea-query", type=Path, required=True)
    parser.add_argument("--main-survey", type=Path, required=True)
    parser.add_argument("--bhp-survey", type=Path, required=True)
    parser.add_argument("--area-ses", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workspace-output-dir", type=Path)
    parser.add_argument("--workspace-scrap-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--stage-aggregate", action="store_true")
    parser.add_argument("--no-write-regenie-inputs", action="store_true")
    parser.add_argument("--seed", type=int, default=env_int("DIRECT_XGB_SEED", 2026))
    parser.add_argument("--threads", type=int, default=env_int("DIRECT_XGB_THREADS", max(1, os.cpu_count() or 1)))
    parser.add_argument("--eta", type=float, default=env_float("DIRECT_XGB_ETA", 0.05))
    parser.add_argument("--max-depth", type=int, default=env_int("DIRECT_XGB_MAX_DEPTH", 6))
    parser.add_argument("--min-child-weight", type=float, default=env_float("DIRECT_XGB_MIN_CHILD_WEIGHT", 20.0))
    parser.add_argument("--reg-lambda", type=float, default=env_float("DIRECT_XGB_LAMBDA", 1.0))
    parser.add_argument("--alpha", type=float, default=env_float("DIRECT_XGB_ALPHA", 0.0))
    parser.add_argument("--subsample", type=float, default=env_float("DIRECT_XGB_SUBSAMPLE", 0.8))
    parser.add_argument("--colsample-bytree", type=float, default=env_float("DIRECT_XGB_COLSAMPLE_BYTREE", 0.8))
    parser.add_argument("--num-boost-round", type=int, default=env_int("DIRECT_XGB_NUM_BOOST_ROUND", 2000))
    parser.add_argument("--early-stopping-rounds", type=int, default=env_int("DIRECT_XGB_EARLY_STOPPING_ROUNDS", 50))
    parser.add_argument("--cv-folds", type=int, default=env_int("DIRECT_XGB_CV_FOLDS", 4))
    parser.add_argument("--min-train-samples", type=int, default=env_int("DIRECT_XGB_MIN_TRAIN_SAMPLES", 1000))
    return parser.parse_args()


def sha256_file(path: Path, max_bytes: int | None = None) -> str | None:
    if not path.exists():
        return None
    if max_bytes is not None and path.stat().st_size > max_bytes:
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def path_record(path: Path, *, hash_limit: int | None = 64 * 1024 * 1024) -> dict[str, object]:
    digest = None if hash_limit is None else sha256_file(path, max_bytes=hash_limit)
    return {
        "path": str(path),
        "size_bytes": int(path.stat().st_size) if path.exists() else None,
        "mtime_ns": int(path.stat().st_mtime_ns) if path.exists() else None,
        "sha256": digest,
    }


def feature_columns_sha256(feature_columns: list[str]) -> str:
    return hashlib.sha256("\n".join(feature_columns).encode("utf-8")).hexdigest()


def write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: list[str]) -> None:
    pd.DataFrame(list(rows), columns=columns).to_csv(path, sep="\t", index=False)


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")


def zscore_from_mask(values: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, float, float]:
    vals = values[mask & np.isfinite(values)]
    require(len(vals) >= MIN_GROUP_N, "too few values for z-scoring")
    mean = float(np.mean(vals))
    sd = float(np.std(vals, ddof=1))
    require(math.isfinite(sd) and sd > 0, "invalid z-score reference SD")
    return (values - mean) / sd, mean, sd


def pearson(x: np.ndarray | pd.Series, y: np.ndarray | pd.Series) -> float:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if int(mask.sum()) < 3:
        return math.nan
    return float(np.corrcoef(x_arr[mask], y_arr[mask])[0, 1])


def spearman(x: np.ndarray | pd.Series, y: np.ndarray | pd.Series) -> float:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if int(mask.sum()) < 3:
        return math.nan
    xr = pd.Series(x_arr[mask]).rank(method="average").to_numpy(float)
    yr = pd.Series(y_arr[mask]).rank(method="average").to_numpy(float)
    return float(np.corrcoef(xr, yr)[0, 1])


def describe(values: np.ndarray) -> dict[str, object]:
    vals = values[np.isfinite(values)]
    if len(vals) == 0:
        return {"n": 0, "mean": math.nan, "sd": math.nan, "min": math.nan, "p01": math.nan, "p05": math.nan, "p50": math.nan, "p95": math.nan, "p99": math.nan, "max": math.nan}
    return {
        "n": int(len(vals)),
        "mean": float(np.mean(vals)),
        "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else math.nan,
        "min": float(np.min(vals)),
        "p01": float(np.quantile(vals, 0.01)),
        "p05": float(np.quantile(vals, 0.05)),
        "p50": float(np.quantile(vals, 0.50)),
        "p95": float(np.quantile(vals, 0.95)),
        "p99": float(np.quantile(vals, 0.99)),
        "max": float(np.max(vals)),
    }


def load_base_scores(args: argparse.Namespace) -> pd.DataFrame:
    all_scores = pd.read_csv(args.all_scores_file, sep="\t", dtype={"FID": str, "IID": str, "role": str, "fold_id": str})
    required = {"FID", "IID", "role", "fold_id", "teacher_z", "ses_ea_proxy_z", "final_model_train_allowed"}
    missing = required - set(all_scores.columns)
    require(not missing, f"all_scores missing columns: {sorted(missing)}")
    require(all_scores["IID"].is_unique, "all_scores has duplicate IID rows")
    for col in ["teacher_z", "ses_ea_proxy_z", "score_raw", "ea_years", "final_model_train_allowed"]:
        if col in all_scores.columns:
            all_scores[col] = numeric(all_scores, col)
    all_scores["final_model_train_allowed"] = all_scores["final_model_train_allowed"].fillna(0).astype(int)

    task_cols = ["IID", "gradcpt_perf_z_age_sex", "flanker_efficiency_z_age_sex"]
    task = pd.read_csv(args.task_score_file, sep="\t", dtype={"IID": str}, usecols=task_cols)
    require(task["IID"].is_unique, "task-score file has duplicate IID rows")
    for col in task_cols[1:]:
        task[col] = numeric(task, col)
    out = all_scores.merge(task, on="IID", how="left", validate="one_to_one")

    with args.fine_tuned_score_file.open() as handle:
        ft_header = handle.readline().rstrip("\n").split("\t")
    ft_cols = ["IID", "gradcpt_flanker_finetuned_ea_proxy_z"]
    if "gradcpt_flanker_teacher_ea_calibrated_proxy_z" in ft_header:
        ft_cols.append("gradcpt_flanker_teacher_ea_calibrated_proxy_z")
    ft = pd.read_csv(args.fine_tuned_score_file, sep="\t", dtype={"IID": str}, usecols=ft_cols)
    require(ft["IID"].is_unique, "fine-tuned score file has duplicate IID rows")
    for col in ft_cols[1:]:
        ft[col] = numeric(ft, col)
    out = out.merge(ft, on="IID", how="left", validate="one_to_one")
    return out


def build_gradcpt_flanker_targets(df: pd.DataFrame, diag_dir: Path) -> pd.DataFrame:
    g = df["gradcpt_perf_z_age_sex"].to_numpy(float)
    f = df["flanker_efficiency_z_age_sex"].to_numpy(float)
    has_g = np.isfinite(g)
    has_f = np.isfinite(f)
    both = has_g & has_f
    either = has_g | has_f
    require(int(both.sum()) >= MIN_GROUP_N, "too few GradCPT+Flanker complete cases")
    r = pearson(g[both], f[both])
    require(math.isfinite(r) and r > 0, f"GradCPT/Flanker correlation must be positive; got {r}")
    loading = math.sqrt(r)
    psi = 1.0 - r
    require(psi > 0, "invalid equal-loading uniqueness")

    g_z, g_mu, g_sd = zscore_from_mask(g, both)
    f_z, f_mu, f_sd = zscore_from_mask(f, both)
    x_all = np.column_stack([g_z, f_z])
    lambda_vec = np.array([loading, loading], dtype=float)
    psi_vec = np.array([psi, psi], dtype=float)

    g_hat = np.full(len(df), np.nan, dtype=float)
    weights_by_pattern: dict[str, np.ndarray] = {}
    for idx, row in enumerate(x_all):
        obs = np.isfinite(row)
        if int(obs.sum()) == 0:
            continue
        lam = lambda_vec[obs]
        ps = psi_vec[obs]
        sigma = np.outer(lam, lam) + np.diag(ps)
        weights = np.linalg.solve(sigma, lam)
        g_hat[idx] = float(row[obs].dot(weights))
        pattern = pattern_name(bool(obs[0]), bool(obs[1]))
        weights_full = np.zeros(2, dtype=float)
        weights_full[obs] = weights
        weights_by_pattern[pattern] = weights_full

    factor_z, factor_mu, factor_sd = zscore_from_mask(g_hat, both)
    mean_raw = np.full(len(df), np.nan, dtype=float)
    mean_raw[both] = (g[both] + f[both]) / 2.0
    mean_z, mean_mu, mean_sd = zscore_from_mask(mean_raw, both)

    df["has_gradcpt"] = has_g.astype(int)
    df["has_flanker"] = has_f.astype(int)
    df["has_gradcpt_or_flanker"] = either
    df["has_gradcpt_and_flanker"] = both
    df["task_pattern_two_domain"] = [pattern_name(bool(a), bool(b)) for a, b in zip(has_g, has_f)]
    df["gradcpt_flanker_factor_hat"] = g_hat
    df["gradcpt_flanker_factor_z"] = factor_z
    df["gradcpt_flanker_mean_raw"] = mean_raw
    df["gradcpt_flanker_mean_z"] = mean_z
    df["used_as_direct_xgb_label"] = np.isfinite(factor_z)

    summary_rows = [
        {"metric": "n_gradcpt", "value": int(has_g.sum())},
        {"metric": "n_flanker", "value": int(has_f.sum())},
        {"metric": "n_both", "value": int(both.sum())},
        {"metric": "n_either", "value": int(either.sum())},
        {"metric": "gradcpt_flanker_complete_case_r", "value": r},
        {"metric": "equal_loading_sqrt_r", "value": loading},
        {"metric": "equal_loading_uniqueness_1_minus_r", "value": psi},
        {"metric": "gradcpt_complete_case_mean", "value": g_mu},
        {"metric": "gradcpt_complete_case_sd", "value": g_sd},
        {"metric": "flanker_complete_case_mean", "value": f_mu},
        {"metric": "flanker_complete_case_sd", "value": f_sd},
        {"metric": "factor_hat_complete_case_mean", "value": factor_mu},
        {"metric": "factor_hat_complete_case_sd", "value": factor_sd},
        {"metric": "mean_raw_complete_case_mean", "value": mean_mu},
        {"metric": "mean_raw_complete_case_sd", "value": mean_sd},
    ]
    pd.DataFrame(summary_rows).to_csv(diag_dir / "gradcpt_flanker_direct_target_summary.tsv", sep="\t", index=False)

    weight_rows = []
    for pattern in ["GRADCPT_FLANKER", "GRADCPT_ONLY", "FLANKER_ONLY"]:
        weights = weights_by_pattern.get(pattern, np.array([math.nan, math.nan]))
        weight_rows.append({"task_pattern": pattern, "weight_gradcpt": weights[0], "weight_flanker": weights[1]})
    pd.DataFrame(weight_rows).to_csv(diag_dir / "gradcpt_flanker_factor_weights_by_pattern.tsv", sep="\t", index=False)

    pattern_rows = []
    for pattern, sub in df.groupby("task_pattern_two_domain", dropna=False):
        row = {"task_pattern": pattern, "n": int(len(sub))}
        row.update({f"factor_z_{k}": v for k, v in describe(sub["gradcpt_flanker_factor_z"].to_numpy(float)).items() if k != "n"})
        pattern_rows.append(row)
    pd.DataFrame(pattern_rows).to_csv(diag_dir / "gradcpt_flanker_direct_target_pattern_counts.tsv", sep="\t", index=False)
    return df


def pattern_name(has_gradcpt: bool, has_flanker: bool) -> str:
    if has_gradcpt and has_flanker:
        return "GRADCPT_FLANKER"
    if has_gradcpt:
        return "GRADCPT_ONLY"
    if has_flanker:
        return "FLANKER_ONLY"
    return "NONE"


def load_saved_feature_columns(path: Path) -> list[str]:
    feature_columns = json.loads(path.read_text())
    require(isinstance(feature_columns, list) and all(isinstance(x, str) for x in feature_columns), "invalid saved feature columns JSON")
    return feature_columns


def append_education_features(x_base: np.ndarray, base_columns: list[str], ea_query: Path, iids: list[str]) -> tuple[np.ndarray, list[str], pd.DataFrame]:
    ea = pd.read_csv(ea_query, dtype={"IID": str})
    require({"IID", "answer_concept_id", "ea_years"} <= set(ea.columns), "ea_query missing required columns")
    ea["answer_concept_id"] = pd.to_numeric(ea["answer_concept_id"], errors="coerce").astype("Int64")
    ea["ea_years"] = pd.to_numeric(ea["ea_years"], errors="coerce")
    ea = ea.drop_duplicates("IID", keep="last").set_index("IID")
    answer = ea.reindex(iids)["answer_concept_id"]
    missing = answer.isna()
    require(int(missing.sum()) == 0, f"{int(missing.sum())} cohort rows missing education answer feature")
    answer_ids = answer.astype(int).to_numpy()
    years = np.array([ses.EA_MAPPING[int(a)] for a in answer_ids], dtype=np.float32)

    extra_cols = ["q1585940_highest_grade_revised_ea_years_num"]
    extra_arrays = [years]
    for aid in EDU_ANSWER_IDS:
        extra_cols.append(f"q1585940_a{aid}_highest_grade")
        extra_arrays.append((answer_ids == aid).astype(np.float32))
    x_extra = np.vstack(extra_arrays).T.astype(np.float32, copy=False)
    x = np.column_stack([x_base, x_extra]).astype(np.float32, copy=False)
    columns = base_columns + extra_cols

    counts = pd.DataFrame({"answer_concept_id": answer_ids}).value_counts().reset_index(name="n")
    counts["feature_column"] = counts["answer_concept_id"].map({aid: f"q1585940_a{aid}_highest_grade" for aid in EDU_ANSWER_IDS})
    counts = counts.sort_values("answer_concept_id")
    counts["numeric_feature_column"] = "q1585940_highest_grade_revised_ea_years_num"
    counts["ea_years_feature_value"] = counts["answer_concept_id"].map(ses.EA_MAPPING)
    return x, columns, counts


def rebuild_feature_matrix_plus_education(args: argparse.Namespace, iids: list[str], diag_dir: Path) -> tuple[np.ndarray, list[str], str]:
    saved_columns = load_saved_feature_columns(args.ses_ea_dir / "xgboost_feature_columns.json")
    metadata = ses.load_local_metadata(str(args.metadata))
    feature_state = {
        "survey_seen": defaultdict(set),
        "survey_age": defaultdict(dict),
        "question_meta": {},
        "answered_qids": defaultdict(set),
        "selected_answers": defaultdict(lambda: defaultdict(set)),
        "answer_names": {},
        "numeric_values": defaultdict(dict),
        "ordinal_answer_values": {},
        "branch_ordinal_values": defaultdict(dict),
        "branch_numeric_values": defaultdict(dict),
        "pmi_missing_iids": defaultdict(set),
        "nonresponse": defaultdict(int),
        "question_encoding": {},
    }
    main_rows = ses.add_survey_feature_rows(str(args.main_survey), feature_state, metadata, bhp=False)
    bhp_rows = ses.add_survey_feature_rows(str(args.bhp_survey), feature_state, metadata, bhp=True)
    branch_recode_rows = ses.apply_lifestyle_branch_recodes(feature_state)
    sex_map = sex_map_from_all_scores(args.all_scores_file)
    area_ses = ses.load_area_ses(str(args.area_ses))
    x_base, base_columns = ses.build_feature_matrix(iids, sex_map, area_ses, feature_state)
    require(base_columns == saved_columns, "rebuilt base feature columns do not match SES-EA proxy feature contract")
    x, columns, edu_counts = append_education_features(x_base, base_columns, args.ea_query, iids)
    feature_hash = feature_columns_sha256(columns)

    pd.DataFrame({"feature_index": range(len(columns)), "feature": columns}).to_csv(
        args.output_dir / "xgboost_feature_columns.tsv", sep="\t", index=False
    )
    (args.output_dir / "xgboost_feature_columns.json").write_text(json.dumps(columns, indent=2) + "\n")
    edu_counts.to_csv(diag_dir / "direct_xgb_education_feature_counts.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"metric": "base_feature_columns", "value": len(base_columns)},
            {"metric": "added_education_feature_columns", "value": len(columns) - len(base_columns)},
            {"metric": "total_feature_columns", "value": len(columns)},
            {"metric": "base_feature_columns_sha256", "value": feature_columns_sha256(base_columns)},
            {"metric": "expanded_feature_columns_sha256", "value": feature_hash},
            {"metric": "main_survey_rows_read", "value": main_rows},
            {"metric": "bhp_survey_rows_read", "value": bhp_rows},
            {"metric": "lifestyle_branch_recode_rules_with_samples", "value": len(branch_recode_rows)},
            {"metric": "lifestyle_branch_recoded_sample_feature_cells", "value": sum(int(row["recoded_samples"]) for row in branch_recode_rows)},
        ]
    ).to_csv(diag_dir / "direct_xgb_feature_counts.tsv", sep="\t", index=False)
    return x, columns, feature_hash


def sex_map_from_all_scores(path: Path) -> dict[str, int]:
    df = pd.read_csv(path, sep="\t", dtype={"IID": str}, usecols=["IID"])
    covar_path = path.parent / "base_covar.txt"
    require(covar_path.exists(), f"missing base covariate file for sex feature: {covar_path}")
    covar = pd.read_csv(covar_path, sep="\t", dtype={"IID": str}, usecols=["IID", "sex_c"])
    covar["sex_c"] = pd.to_numeric(covar["sex_c"], errors="coerce")
    covar["genetic_sex_01"] = covar["sex_c"] + 0.5
    vals = covar.set_index("IID")["genetic_sex_01"].to_dict()
    missing = set(df["IID"]) - set(vals)
    require(not missing, f"base_covar missing sex rows for {len(missing)} IIDs")
    return {iid: int(round(float(vals[iid]))) for iid in df["IID"].astype(str)}


def xgb_params(args: argparse.Namespace) -> dict[str, object]:
    return {
        "objective": "reg:squarederror",
        "eta": args.eta,
        "max_depth": args.max_depth,
        "min_child_weight": args.min_child_weight,
        "subsample": args.subsample,
        "colsample_bytree": args.colsample_bytree,
        "lambda": args.reg_lambda,
        "alpha": args.alpha,
        "tree_method": "hist",
        "seed": args.seed,
        "nthread": args.threads,
        "eval_metric": "rmse",
    }


def dmatrix(x: np.ndarray, feature_columns: list[str], label: np.ndarray | None = None):
    import xgboost as xgb
    return xgb.DMatrix(x, label=label, feature_names=feature_columns, missing=np.nan)


def aligned_label(target: np.ndarray, proxy: np.ndarray, train_idx: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    y = target[train_idx]
    p = proxy[train_idx]
    mask = np.isfinite(y) & np.isfinite(p)
    require(int(mask.sum()) == len(train_idx), "nonfinite values in direct XGB training label/proxy scale")
    y_mean = float(np.mean(y))
    y_sd = float(np.std(y, ddof=1))
    p_mean = float(np.mean(p))
    p_sd = float(np.std(p, ddof=1))
    require(math.isfinite(y_sd) and y_sd > 0, "invalid target training SD")
    require(math.isfinite(p_sd) and p_sd > 0, "invalid proxy training SD")
    out = p_mean + p_sd * ((y - y_mean) / y_sd)
    return out, {
        "target_train_mean": y_mean,
        "target_train_sd": y_sd,
        "ses_ea_proxy_train_mean": p_mean,
        "ses_ea_proxy_train_sd": p_sd,
    }


def train_one_scratch_model(
    *,
    args: argparse.Namespace,
    model_name: str,
    x: np.ndarray,
    feature_columns: list[str],
    feature_hash: str,
    train_idx: np.ndarray,
    pred_idx: np.ndarray,
    target: np.ndarray,
    proxy: np.ndarray,
    selected_rounds_override: int | None = None,
):
    import xgboost as xgb

    require(len(train_idx) >= args.min_train_samples, f"{model_name} has too few direct XGB labels: {len(train_idx)}")
    y_aligned, align = aligned_label(target, proxy, train_idx)
    params = xgb_params(args)
    if selected_rounds_override is None:
        cv = xgb.cv(
            params,
            dmatrix(x[train_idx, :], feature_columns, label=y_aligned),
            num_boost_round=args.num_boost_round,
            nfold=min(args.cv_folds, max(2, len(train_idx) // 2)),
            early_stopping_rounds=args.early_stopping_rounds,
            seed=args.seed + stable_model_seed_offset(model_name),
            verbose_eval=False,
        )
        selected_rounds = int(len(cv))
        best_rmse = float(cv["test-rmse-mean"].iloc[-1]) if "test-rmse-mean" in cv.columns else math.nan
    else:
        selected_rounds = int(selected_rounds_override)
        best_rmse = math.nan
    model = xgb.train(
        params,
        dmatrix(x[train_idx, :], feature_columns, label=y_aligned),
        num_boost_round=selected_rounds,
        verbose_eval=False,
    )
    pred = model.predict(dmatrix(x[pred_idx, :], feature_columns))
    model.set_attr(
        role="direct_gradcpt_flanker_xgb",
        model_name=model_name,
        feature_columns_sha256=feature_hash,
        train_samples=str(len(train_idx)),
        prediction_samples=str(len(pred_idx)),
        selected_rounds=str(selected_rounds),
        target_column="gradcpt_flanker_factor_z",
        target_alignment="fold_specific_ses_ea_proxy_z_mean_sd",
        seed=str(args.seed),
        eta=str(args.eta),
        max_depth=str(args.max_depth),
        min_child_weight=str(args.min_child_weight),
        reg_lambda=str(args.reg_lambda),
    )
    row = {
        "model_name": model_name,
        "train_n": int(len(train_idx)),
        "predict_n": int(len(pred_idx)),
        "selected_rounds": int(selected_rounds),
        "best_cv_test_rmse": best_rmse,
        **align,
        "feature_columns_sha256": feature_hash,
    }
    return pred, row, model


def stable_model_seed_offset(name: str) -> int:
    if name.startswith("fold_"):
        return int(name.split("_")[1])
    return 99


def train_direct_xgb(args: argparse.Namespace, df: pd.DataFrame, x: np.ndarray, feature_columns: list[str], feature_hash: str, diag_dir: Path) -> pd.DataFrame:
    model_dir = args.output_dir / "xgboost_models"
    model_dir.mkdir(parents=True, exist_ok=True)
    target = df["gradcpt_flanker_factor_z"].to_numpy(float)
    proxy = df["ses_ea_proxy_z"].to_numpy(float)
    label_mask = np.isfinite(target) & np.isfinite(proxy)
    role = df["role"].astype(str)
    fold = df["fold_id"].astype(str)
    oof = role.eq("oof").to_numpy(bool)
    applied = role.eq("applied").to_numpy(bool)
    final_allowed = df["final_model_train_allowed"].astype(bool).to_numpy(bool)
    pred_raw = np.full(len(df), np.nan, dtype=float)
    pred_model = np.full(len(df), "", dtype=object)
    rows = []

    for k in range(5):
        train_mask = oof & ~fold.eq(str(k)).to_numpy(bool) & label_mask
        pred_mask = oof & fold.eq(str(k)).to_numpy(bool)
        pred, row, model = train_one_scratch_model(
            args=args,
            model_name=f"fold_{k}",
            x=x,
            feature_columns=feature_columns,
            feature_hash=feature_hash,
            train_idx=np.flatnonzero(train_mask),
            pred_idx=np.flatnonzero(pred_mask),
            target=target,
            proxy=proxy,
        )
        pred_raw[pred_mask] = pred
        pred_model[pred_mask] = f"fold_{k}_direct_xgb"
        model.save_model(str(model_dir / f"fold_{k}_direct_xgb.json"))
        rows.append(row)
        print(f"fold_{k}: train_n={row['train_n']} predict_n={row['predict_n']} selected_rounds={row['selected_rounds']}", flush=True)

    selected_final_rounds = int(np.median([int(row["selected_rounds"]) for row in rows]))
    final_train_mask = oof & final_allowed & label_mask
    final_pred_mask = applied
    pred, row, model = train_one_scratch_model(
        args=args,
        model_name="final_model",
        x=x,
        feature_columns=feature_columns,
        feature_hash=feature_hash,
        train_idx=np.flatnonzero(final_train_mask),
        pred_idx=np.flatnonzero(final_pred_mask),
        target=target,
        proxy=proxy,
        selected_rounds_override=selected_final_rounds,
    )
    pred_raw[final_pred_mask] = pred
    pred_model[final_pred_mask] = "final_model_direct_xgb"
    model.save_model(str(model_dir / "final_model_direct_xgb.json"))
    row["selected_rounds"] = selected_final_rounds
    rows.append(row)
    print(f"final_model: train_n={row['train_n']} predict_n={row['predict_n']} selected_rounds={row['selected_rounds']}", flush=True)

    require(np.all(np.isfinite(pred_raw)), "not every cohort row received a direct XGB prediction")
    mean = float(np.mean(pred_raw))
    sd = float(np.std(pred_raw, ddof=1))
    require(math.isfinite(sd) and sd > 0, "invalid direct XGB prediction SD")
    df["gradcpt_flanker_direct_xgb_proxy_raw"] = pred_raw
    df["gradcpt_flanker_direct_xgb_proxy_z"] = (pred_raw - mean) / sd
    df["direct_xgb_prediction_model_name"] = pred_model
    pd.DataFrame(rows).to_csv(diag_dir / "direct_xgb_round_selection.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(args.output_dir / "direct_xgb_model_manifest.tsv", sep="\t", index=False)
    return df


def add_four_variable_calibration(df: pd.DataFrame, args: argparse.Namespace, diag_dir: Path) -> pd.DataFrame:
    feature_cols = [
        "teacher_z",
        "ses_ea_proxy_z",
        "gradcpt_flanker_finetuned_ea_proxy_z",
        "gradcpt_flanker_direct_xgb_proxy_z",
    ]
    target_col = "gradcpt_flanker_mean_z"
    feature_mat = df[feature_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    target = pd.to_numeric(df[target_col], errors="coerce").to_numpy(float)
    features_finite = np.all(np.isfinite(feature_mat), axis=1)
    labels_finite = features_finite & np.isfinite(target)
    role = df["role"].astype(str)
    fold = df["fold_id"].astype(str)
    oof = role.eq("oof").to_numpy(bool)
    applied = role.eq("applied").to_numpy(bool)
    final_allowed = df["final_model_train_allowed"].astype(bool).to_numpy(bool)
    pred = np.full(len(df), np.nan, dtype=float)
    coef_rows = []

    scale_rows = []
    groups = {
        "calibration_labeled_all": labels_finite,
        "calibration_labeled_oof": labels_finite & oof,
        "calibration_labeled_final_allowed": labels_finite & oof & final_allowed,
        "calibration_labeled_applied": labels_finite & applied,
        "finite_features_full_cohort": features_finite,
    }
    for group_name, mask in groups.items():
        for col_idx, col in enumerate(feature_cols):
            row = {"group": group_name, "variable": col}
            row.update(describe(feature_mat[mask, col_idx]))
            scale_rows.append(row)
        row = {"group": group_name, "variable": target_col}
        row.update(describe(target[mask]))
        scale_rows.append(row)
    pd.DataFrame(scale_rows).to_csv(diag_dir / "four_variable_calibration_feature_scales.tsv", sep="\t", index=False)

    def fit_predict(train_mask: np.ndarray, pred_mask: np.ndarray, fit_name: str, predict_group: str) -> None:
        require(int(train_mask.sum()) >= args.min_train_samples, f"too few labels for {fit_name} calibration")
        design_train = np.column_stack([np.ones(int(train_mask.sum())), feature_mat[train_mask]])
        coef, *_ = np.linalg.lstsq(design_train, target[train_mask], rcond=None)
        if int(pred_mask.sum()) > 0:
            design_pred = np.column_stack([np.ones(int(pred_mask.sum())), feature_mat[pred_mask]])
            pred[pred_mask] = design_pred.dot(coef)
        row = {
            "fit": fit_name,
            "predict_group": predict_group,
            "train_n": int(train_mask.sum()),
            "predict_n": int(pred_mask.sum()),
            "intercept": float(coef[0]),
        }
        for col, value in zip(feature_cols, coef[1:]):
            row[f"coef_{col}"] = float(value)
        coef_rows.append(row)

    for k in range(5):
        train_mask = oof & ~fold.eq(str(k)).to_numpy(bool) & labels_finite
        pred_mask = oof & fold.eq(str(k)).to_numpy(bool) & features_finite
        fit_predict(train_mask, pred_mask, f"oof_fold_{k}_train_other_folds", f"oof_fold_{k}")

    train_mask = oof & final_allowed & labels_finite
    pred_mask = applied & features_finite
    fit_predict(train_mask, pred_mask, "applied_model_train_kinholdout_oof", "applied")

    require(np.all(np.isfinite(pred[features_finite])), "some finite-feature rows did not receive four-variable calibration")
    mean = float(np.mean(pred[features_finite]))
    sd = float(np.std(pred[features_finite], ddof=1))
    require(math.isfinite(sd) and sd > 0, "invalid four-variable calibrated SD")
    df["gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_raw"] = pred
    df["gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_z"] = (pred - mean) / sd
    coef_df = pd.DataFrame(coef_rows)
    coef_df["calibrated_raw_mean_full_cohort"] = mean
    coef_df["calibrated_raw_sd_full_cohort"] = sd
    coef_df.to_csv(diag_dir / "four_variable_calibration_coefficients.tsv", sep="\t", index=False)
    return df


def correlation_rows(df: pd.DataFrame, diag_dir: Path) -> None:
    score_cols = [
        "ses_ea_proxy_z",
        "gradcpt_flanker_finetuned_ea_proxy_z",
        "gradcpt_flanker_direct_xgb_proxy_z",
        "gradcpt_flanker_teacher_ea_calibrated_proxy_z",
        "gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_z",
    ]
    target_cols = [
        "gradcpt_flanker_mean_z",
        "gradcpt_flanker_factor_z",
        "gradcpt_perf_z_age_sex",
        "flanker_efficiency_z_age_sex",
        "teacher_z",
        "ses_ea_proxy_z",
    ]
    role = df["role"].astype(str)
    groups = {
        "combined_both_target": np.isfinite(df["gradcpt_flanker_mean_z"].to_numpy(float)),
        "oof_both_target": role.eq("oof").to_numpy(bool) & np.isfinite(df["gradcpt_flanker_mean_z"].to_numpy(float)),
        "applied_both_target": role.eq("applied").to_numpy(bool) & np.isfinite(df["gradcpt_flanker_mean_z"].to_numpy(float)),
        "combined_either_target": np.isfinite(df["gradcpt_flanker_factor_z"].to_numpy(float)),
        "oof_either_target": role.eq("oof").to_numpy(bool) & np.isfinite(df["gradcpt_flanker_factor_z"].to_numpy(float)),
        "applied_either_target": role.eq("applied").to_numpy(bool) & np.isfinite(df["gradcpt_flanker_factor_z"].to_numpy(float)),
        "full_cohort": np.ones(len(df), dtype=bool),
    }
    rows = []
    for group_name, mask in groups.items():
        for score_col in score_cols:
            if score_col not in df.columns:
                continue
            x = df[score_col].to_numpy(float)
            for target_col in target_cols:
                if target_col not in df.columns:
                    continue
                y = df[target_col].to_numpy(float)
                finite = mask & np.isfinite(x) & np.isfinite(y)
                rows.append(
                    {
                        "group": group_name,
                        "score": score_col,
                        "target": target_col,
                        "n": int(finite.sum()),
                        "pearson_r": pearson(x[finite], y[finite]),
                        "spearman_r": spearman(x[finite], y[finite]),
                    }
                )
    pd.DataFrame(rows).to_csv(diag_dir / "direct_xgb_and_calibration_correlations.tsv", sep="\t", index=False)


def score_distributions(df: pd.DataFrame, diag_dir: Path) -> None:
    cols = [
        "gradcpt_flanker_direct_xgb_proxy_z",
        "gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_z",
        "gradcpt_flanker_factor_z",
        "gradcpt_flanker_mean_z",
    ]
    groups = {
        "full_cohort": np.ones(len(df), dtype=bool),
        "role_oof": df["role"].eq("oof").to_numpy(bool),
        "role_applied": df["role"].eq("applied").to_numpy(bool),
        "both_gradcpt_flanker": df["has_gradcpt_and_flanker"].to_numpy(bool),
        "either_gradcpt_flanker": df["has_gradcpt_or_flanker"].to_numpy(bool),
        "gradcpt_only": df["task_pattern_two_domain"].eq("GRADCPT_ONLY").to_numpy(bool),
        "flanker_only": df["task_pattern_two_domain"].eq("FLANKER_ONLY").to_numpy(bool),
    }
    rows = []
    for group_name, mask in groups.items():
        for col in cols:
            row = {"group": group_name, "variable": col}
            row.update(describe(df.loc[mask, col].to_numpy(float)))
            rows.append(row)
    pd.DataFrame(rows).to_csv(diag_dir / "direct_xgb_score_distributions.tsv", sep="\t", index=False)


def write_outputs(args: argparse.Namespace, df: pd.DataFrame, diag_dir: Path) -> None:
    out_cols = [
        "FID",
        "IID",
        "role",
        "fold_id",
        "final_model_train_allowed",
        "teacher_z",
        "ses_ea_proxy_z",
        "gradcpt_perf_z_age_sex",
        "flanker_efficiency_z_age_sex",
        "has_gradcpt",
        "has_flanker",
        "has_gradcpt_or_flanker",
        "has_gradcpt_and_flanker",
        "task_pattern_two_domain",
        "gradcpt_flanker_factor_hat",
        "gradcpt_flanker_factor_z",
        "gradcpt_flanker_mean_raw",
        "gradcpt_flanker_mean_z",
        "used_as_direct_xgb_label",
        "gradcpt_flanker_finetuned_ea_proxy_z",
        "gradcpt_flanker_direct_xgb_proxy_raw",
        "gradcpt_flanker_direct_xgb_proxy_z",
        "direct_xgb_prediction_model_name",
        "gradcpt_flanker_teacher_ea_calibrated_proxy_z",
        "gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_raw",
        "gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_z",
    ]
    existing = [col for col in out_cols if col in df.columns]
    df[existing].to_csv(args.output_dir / "gradcpt_flanker_direct_xgb_proxy_scores_wide.tsv", sep="\t", index=False)
    df.loc[df["role"].eq("oof"), existing].to_csv(args.output_dir / "gradcpt_flanker_direct_xgb_proxy_oof_scores.tsv", sep="\t", index=False)
    df.loc[df["role"].eq("applied"), existing].to_csv(args.output_dir / "gradcpt_flanker_direct_xgb_proxy_applied_scores.tsv", sep="\t", index=False)

    if not args.no_write_regenie_inputs:
        phen_col = "gradcpt_flanker_direct_xgb_proxy_z"
        df[["FID", "IID", phen_col]].to_csv(args.output_dir / "phen.txt", sep="\t", index=False)
        for name in ["base_covar.txt", "covar.txt"]:
            src = args.ses_ea_dir / name
            if src.exists():
                shutil.copy2(src, args.output_dir / name)
        df[["FID", "IID"]].to_csv(args.output_dir / "training_iids.txt", sep="\t", index=False, header=False)

    aggregate_files = [
        "direct_xgb_round_selection.tsv",
        "direct_xgb_and_calibration_correlations.tsv",
        "direct_xgb_score_distributions.tsv",
        "four_variable_calibration_coefficients.tsv",
        "four_variable_calibration_feature_scales.tsv",
        "gradcpt_flanker_direct_target_summary.tsv",
        "gradcpt_flanker_direct_target_pattern_counts.tsv",
        "gradcpt_flanker_factor_weights_by_pattern.tsv",
        "direct_xgb_feature_counts.tsv",
        "direct_xgb_education_feature_counts.tsv",
    ]
    pd.DataFrame({"aggregate_diagnostic": aggregate_files}).to_csv(diag_dir / "direct_xgb_aggregate_diagnostic_manifest.tsv", sep="\t", index=False)


def copy_outputs(args: argparse.Namespace, diag_dir: Path) -> None:
    if args.workspace_output_dir is not None:
        output_uri = workspace_path_to_gs_uri(args.workspace_output_dir)
        if output_uri:
            gcloud_cp_dir_contents(args.output_dir, output_uri)
        else:
            args.workspace_output_dir.mkdir(parents=True, exist_ok=True)
            for path in args.output_dir.iterdir():
                dest = args.workspace_output_dir / path.name
                copy_path_retry(path, dest)
        print(f"Copied direct XGB outputs to {args.workspace_output_dir}", flush=True)

    if args.stage_aggregate:
        require(args.workspace_scratch_dir is not None, "--stage-aggregate requires --workspace-scrap-dir")
        scratch_uri = workspace_path_to_gs_uri(args.workspace_scratch_dir)
        if scratch_uri:
            gcloud_cp_files(
                list(diag_dir.glob("*.tsv"))
                + [p for p in [args.output_dir / "direct_xgb_model_manifest.tsv", args.output_dir / "xgboost_feature_columns.tsv"] if p.exists()],
                scratch_uri,
            )
        else:
            args.workspace_scratch_dir.mkdir(parents=True, exist_ok=True)
            for path in diag_dir.glob("*.tsv"):
                shutil.copy2(path, args.workspace_scratch_dir / path.name)
            for name in ["direct_xgb_model_manifest.tsv", "xgboost_feature_columns.tsv"]:
                src = args.output_dir / name
                if src.exists():
                    shutil.copy2(src, args.workspace_scratch_dir / name)
        print(f"Staged aggregate diagnostics to {args.workspace_scratch_dir}", flush=True)


def workspace_bucket_uri_from_mount() -> str | None:
    try:
        output = subprocess.check_output(["mount"], text=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    needle = " on /home/jupyter/workspace/workspace-bucket "
    for line in output.splitlines():
        if needle in line:
            source = line.split()[0]
            return source if source.startswith("gs://") else f"gs://{source}"
    return None


def workspace_path_to_gs_uri(path: Path) -> str | None:
    mount_root = Path("/home/jupyter/workspace/workspace-bucket")
    try:
        rel = path.resolve().relative_to(mount_root)
    except ValueError:
        return None
    bucket = workspace_bucket_uri_from_mount()
    if not bucket:
        return None
    return f"{bucket.rstrip('/')}/{str(rel).strip('/')}"


def gcloud_cp_dir_contents(src_dir: Path, dest_uri: str) -> None:
    cmd = (
        "gcloud storage cp --recursive "
        f"{shlex.quote(str(src_dir))}/* "
        f"{shlex.quote(dest_uri.rstrip('/') + '/')}"
    )
    subprocess.run(cmd, shell=True, check=True)


def gcloud_cp_files(paths: list[Path], dest_uri: str) -> None:
    if not paths:
        return
    quoted = " ".join(shlex.quote(str(path)) for path in paths)
    cmd = f"gcloud storage cp {quoted} {shlex.quote(dest_uri.rstrip('/') + '/')}"
    subprocess.run(cmd, shell=True, check=True)


def copy_path_retry(src: Path, dest: Path, attempts: int = 4) -> None:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            if src.is_dir():
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(src, dest)
            else:
                shutil.copy2(src, dest)
            return
        except OSError as exc:
            last_error = exc
            print(f"Copy attempt {attempt}/{attempts} failed for {dest}: {exc}", flush=True)
            time.sleep(min(30, 2 * attempt))
    raise RuntimeError(f"Failed to copy {src} to {dest}") from last_error


def build_run_signature(args: argparse.Namespace) -> str:
    payload = {
        "paths": {
            "all_scores_file": path_record(args.all_scores_file),
            "task_score_file": path_record(args.task_score_file),
            "fine_tuned_score_file": path_record(args.fine_tuned_score_file),
            "ea_query": path_record(args.ea_query),
            "main_survey": path_record(args.main_survey, hash_limit=None),
            "bhp_survey": path_record(args.bhp_survey, hash_limit=None),
            "area_ses": path_record(args.area_ses),
            "metadata": path_record(args.metadata),
            "setup_script": path_record(Path(__file__).with_name("setup_ses_ea_proxy_gwas.py")),
            "direct_xgb_script": path_record(Path(__file__)),
        },
        "params": {
            "seed": args.seed,
            "threads": args.threads,
            "eta": args.eta,
            "max_depth": args.max_depth,
            "min_child_weight": args.min_child_weight,
            "lambda": args.reg_lambda,
            "alpha": args.alpha,
            "subsample": args.subsample,
            "colsample_bytree": args.colsample_bytree,
            "num_boost_round": args.num_boost_round,
            "early_stopping_rounds": args.early_stopping_rounds,
            "cv_folds": args.cv_folds,
            "min_train_samples": args.min_train_samples,
            "education_feature": "revised_numeric_plus_one_hot",
            "target": "gradcpt_flanker_factor_z_missing_pattern_shrunk",
            "label_alignment": "fold_specific_ses_ea_proxy_z_mean_sd",
        },
    }
    blob = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def maybe_skip(args: argparse.Namespace, run_signature: str) -> bool:
    manifest = args.output_dir / "direct_xgb_runtime_manifest.json"
    expected = [
        args.output_dir / "gradcpt_flanker_direct_xgb_proxy_scores_wide.tsv",
        args.output_dir / "direct_xgb_model_manifest.tsv",
        args.output_dir / "xgboost_models" / "final_model_direct_xgb.json",
        manifest,
    ]
    if args.force:
        if args.output_dir.exists():
            shutil.rmtree(args.output_dir)
        return False
    if all(path.exists() and path.stat().st_size > 0 for path in expected):
        try:
            old = json.loads(manifest.read_text())
        except json.JSONDecodeError:
            return False
        if old.get("run_signature") == run_signature:
            print(f"Existing direct XGB outputs match inputs/parameters; skipping: {args.output_dir}", flush=True)
            return True
    return False


def write_runtime_manifest(args: argparse.Namespace, run_signature: str, df: pd.DataFrame, feature_columns: list[str], feature_hash: str, elapsed: float) -> None:
    manifest = {
        "created_utc": pd.Timestamp.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "elapsed_seconds": elapsed,
        "run_signature": run_signature,
        "cohort_rows": int(len(df)),
        "feature_columns": int(len(feature_columns)),
        "feature_columns_sha256": feature_hash,
        "used_as_direct_xgb_label": int(df["used_as_direct_xgb_label"].sum()),
        "both_gradcpt_flanker": int(df["has_gradcpt_and_flanker"].sum()),
        "final_model_train_allowed": int(df["final_model_train_allowed"].sum()),
        "final_model_used_as_direct_label": int((df["used_as_direct_xgb_label"] & df["final_model_train_allowed"].astype(bool)).sum()),
        "primary_phenotype_z_col": "gradcpt_flanker_direct_xgb_proxy_z",
        "diagnostic_teacher_including_calibration_z_col": "gradcpt_flanker_teacher_ea_finetuned_direct_calibrated_proxy_z",
        "params": xgb_params(args),
        "num_boost_round": args.num_boost_round,
        "early_stopping_rounds": args.early_stopping_rounds,
        "cv_folds": args.cv_folds,
        "education_feature": "revised_numeric_plus_one_hot",
        "label_alignment": "fold_specific_ses_ea_proxy_z_mean_sd",
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
    }
    try:
        import xgboost as xgb
        manifest["xgboost_version"] = xgb.__version__
    except Exception:
        manifest["xgboost_version"] = ""
    (args.output_dir / "direct_xgb_runtime_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = parse_args()
    start = time.time()
    run_signature = build_run_signature(args)
    if maybe_skip(args, run_signature):
        copy_outputs(args, args.output_dir / "diagnostics")
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    diag_dir = args.output_dir / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)

    print("=== Direct GradCPT/Flanker scratch XGBoost proxy ===", flush=True)
    print(f"output_dir={args.output_dir}", flush=True)
    print(f"run_signature={run_signature}", flush=True)
    print(f"xgb_params={xgb_params(args)}", flush=True)

    df = load_base_scores(args)
    df = build_gradcpt_flanker_targets(df, diag_dir)
    x, feature_columns, feature_hash = rebuild_feature_matrix_plus_education(args, df["IID"].astype(str).tolist(), diag_dir)
    require(x.shape[0] == len(df), "feature matrix row count mismatch")
    df = train_direct_xgb(args, df, x, feature_columns, feature_hash, diag_dir)
    df = add_four_variable_calibration(df, args, diag_dir)
    correlation_rows(df, diag_dir)
    score_distributions(df, diag_dir)
    write_outputs(args, df, diag_dir)
    elapsed = time.time() - start
    write_runtime_manifest(args, run_signature, df, feature_columns, feature_hash, elapsed)
    copy_outputs(args, diag_dir)
    print(
        "Direct XGB complete: "
        f"cohort_rows={len(df)} "
        f"direct_labels={int(df['used_as_direct_xgb_label'].sum())} "
        f"both_target={int(df['has_gradcpt_and_flanker'].sum())} "
        f"feature_columns={len(feature_columns)} "
        f"elapsed_seconds={elapsed:.1f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
