#!/usr/bin/env python3
"""Fine-tune saved SES-EA proxy boosters toward ETM-derived cognitive targets.

This command is downstream of the SES-EA proxy and ETM scoring commands.  It
does not query BigQuery, run GWAS, or modify the original SES-EA proxy outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import setup_ses_ea_proxy_gwas as ses


HASH_MAX_BYTES = 100 * 1024 * 1024
MIN_GROUP_N = 20
MODEL_NAMES = [f"fold_{i}" for i in range(5)] + ["final_model"]
TASK_COLS = [
    "dd_patience_z_age_sex",
    "gradcpt_perf_z_age_sex",
    "flanker_efficiency_z_age_sex",
    "emorecog_perf_z_age_sex",
]
SCORE_PREFIX_BY_TARGET = {
    "strong-task-g4": "g4_finetuned_ea_proxy",
    "all-etm-g4": "g4_finetuned_ea_proxy_all_etm_target",
    "gradcpt-flanker-mean": "gradcpt_flanker_finetuned_ea_proxy",
}
CALIBRATED_PREFIX_BY_TARGET = {
    "gradcpt-flanker-mean": "gradcpt_flanker_teacher_ea_calibrated_proxy",
}


def env_float(name: str, default: float) -> float:
    return float(os.environ.get(name, default))


def env_int(name: str, default: int) -> int:
    return int(os.environ.get(name, default))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ses-ea-dir", type=Path, required=True)
    parser.add_argument("--all-scores-file", type=Path, required=True)
    parser.add_argument("--task-score-file", type=Path, required=True)
    parser.add_argument("--etm-g-file", type=Path, default=None)
    parser.add_argument("--ea-query", type=Path, required=True)
    parser.add_argument("--main-survey", type=Path, required=True)
    parser.add_argument("--bhp-survey", type=Path, required=True)
    parser.add_argument("--area-ses", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workspace-output-dir", type=Path, default=None)
    parser.add_argument("--workspace-scrap-dir", type=Path, default=None)
    parser.add_argument(
        "--target",
        choices=("strong-task-g4", "all-etm-g4", "gradcpt-flanker-mean"),
        default="strong-task-g4",
    )
    parser.add_argument("--stage-aggregate", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-write-regenie-inputs", action="store_true")
    parser.add_argument("--eta", type=float, default=env_float("G4_FINETUNE_ETA", 0.01))
    parser.add_argument("--max-depth", type=int, default=env_int("G4_FINETUNE_MAX_DEPTH", 3))
    parser.add_argument("--min-child-weight", type=float, default=env_float("G4_FINETUNE_MIN_CHILD_WEIGHT", 10.0))
    parser.add_argument("--reg-lambda", type=float, default=env_float("G4_FINETUNE_LAMBDA", 2.0))
    parser.add_argument("--alpha", type=float, default=env_float("G4_FINETUNE_ALPHA", 0.0))
    parser.add_argument("--max-rounds", type=int, default=env_int("G4_FINETUNE_MAX_ROUNDS", 500))
    parser.add_argument("--early-stopping-rounds", type=int, default=env_int("G4_FINETUNE_EARLY_STOPPING_ROUNDS", 25))
    parser.add_argument("--valid-fraction", type=float, default=env_float("G4_FINETUNE_VALID_FRACTION", 0.20))
    parser.add_argument("--min-train-samples", type=int, default=env_int("G4_FINETUNE_MIN_TRAIN_SAMPLES", 1000))
    parser.add_argument("--seed", type=int, default=env_int("G4_FINETUNE_SEED", 2026))
    parser.add_argument("--threads", type=int, default=env_int("G4_FINETUNE_THREADS", max(1, os.cpu_count() or 1)))
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"ERROR: {message}")


def target_uses_etm_g(target_mode: str) -> bool:
    return target_mode in ("strong-task-g4", "all-etm-g4")


def read_header(path: Path, sep: str = "\t") -> list[str]:
    require(path.exists() and path.stat().st_size > 0, f"missing or empty file: {path}")
    with path.open() as handle:
        return handle.readline().rstrip("\n").split(sep)


def sha256_file(path: Path) -> str | None:
    if path.stat().st_size > HASH_MAX_BYTES:
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_fingerprint(path: Path) -> dict[str, object]:
    require(path.exists() and path.stat().st_size > 0, f"missing or empty file: {path}")
    stat = path.stat()
    fingerprint = {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "sha256": sha256_file(path),
    }
    if fingerprint["sha256"] is None:
        fingerprint["mtime_ns"] = int(stat.st_mtime_ns)
    return fingerprint


def feature_columns_sha256(feature_columns: list[str]) -> str:
    return hashlib.sha256("\n".join(feature_columns).encode("utf-8")).hexdigest()


def numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")


def pearson(x: Iterable[float], y: Iterable[float]) -> tuple[float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < MIN_GROUP_N:
        return math.nan, n
    return float(np.corrcoef(x[mask], y[mask])[0, 1]), n


def spearman(x: Iterable[float], y: Iterable[float]) -> tuple[float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < MIN_GROUP_N:
        return math.nan, n
    xr = pd.Series(x[mask]).rank(method="average").to_numpy()
    yr = pd.Series(y[mask]).rank(method="average").to_numpy()
    return float(np.corrcoef(xr, yr)[0, 1]), n


def summarize(values: Iterable[float]) -> dict[str, float | int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    out: dict[str, float | int] = {"n": int(x.size)}
    if x.size < MIN_GROUP_N:
        for key in ("mean", "sd", "min", "p01", "p05", "p50", "p95", "p99", "max"):
            out[key] = math.nan
        return out
    qs = np.quantile(x, [0.01, 0.05, 0.50, 0.95, 0.99])
    out.update(
        {
            "mean": float(np.mean(x)),
            "sd": float(np.std(x, ddof=1)) if x.size > 1 else math.nan,
            "min": float(np.min(x)),
            "p01": float(qs[0]),
            "p05": float(qs[1]),
            "p50": float(qs[2]),
            "p95": float(qs[3]),
            "p99": float(qs[4]),
            "max": float(np.max(x)),
        }
    )
    return out


def load_base_tables(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_scores = pd.read_csv(
        args.all_scores_file,
        sep="\t",
        dtype={"FID": str, "IID": str, "role": str, "fold_id": str},
    )
    required = {"FID", "IID", "role", "fold_id", "ea_years", "teacher_z", "score_raw", "ses_ea_proxy_z"}
    missing = required - set(all_scores.columns)
    require(not missing, f"{args.all_scores_file} missing columns: {sorted(missing)}")
    require(all_scores["IID"].is_unique, "all_scores.tsv has duplicate IID rows")
    for col in ("ea_years", "teacher_z", "score_raw", "ses_ea_proxy_z"):
        all_scores[col] = numeric_series(all_scores, col)
    if "final_model_train_allowed" not in all_scores.columns:
        all_scores["final_model_train_allowed"] = all_scores["role"].eq("oof").astype(int)
    all_scores["final_model_train_allowed"] = pd.to_numeric(
        all_scores["final_model_train_allowed"], errors="coerce"
    ).fillna(0).astype(int)

    covar_file = args.ses_ea_dir / "covar.txt"
    require(covar_file.exists(), f"missing covariate file: {covar_file}")
    covars = pd.read_csv(covar_file, sep="\t", dtype={"FID": str, "IID": str})
    require(covars["IID"].is_unique, "covar.txt has duplicate IID rows")
    return all_scores, covars


def load_etm_tables(args: argparse.Namespace, all_scores: pd.DataFrame) -> pd.DataFrame:
    task_header = set(read_header(args.task_score_file))
    task_cols = ["IID"] + [c for c in TASK_COLS if c in task_header]
    missing_task = {"gradcpt_perf_z_age_sex", "flanker_efficiency_z_age_sex"} - set(task_cols)
    require(not missing_task, f"{args.task_score_file} missing columns: {sorted(missing_task)}")
    task = pd.read_csv(args.task_score_file, sep="\t", dtype={"IID": str}, usecols=task_cols)
    require(task["IID"].is_unique, "task-score file has duplicate IID rows")

    out = all_scores.merge(task, on="IID", how="left", validate="one_to_one")
    if args.etm_g_file is not None and args.etm_g_file.exists() and args.etm_g_file.stat().st_size > 0:
        etm_header = set(read_header(args.etm_g_file))
        etm_cols = ["IID"]
        for optional in ("etm_g4_z", "n_tasks_observed_four_domain", "task_pattern_four_domain", "etm_g_z"):
            if optional in etm_header:
                etm_cols.append(optional)
        if len(etm_cols) > 1:
            etm = pd.read_csv(args.etm_g_file, sep="\t", dtype={"IID": str}, usecols=etm_cols)
            require(etm["IID"].is_unique, "ETM-g file has duplicate IID rows")
            out = out.merge(etm, on="IID", how="left", validate="one_to_one")
    if target_uses_etm_g(args.target):
        require(args.etm_g_file is not None, f"{args.target} requires --etm-g-file")
        require("etm_g4_z" in out.columns, f"{args.etm_g_file} missing etm_g4_z")
    for col in ["etm_g4_z", "etm_g_z"] + TASK_COLS:
        if col in out.columns:
            out[col] = numeric_series(out, col)
    out["has_etm_g4"] = np.isfinite(out["etm_g4_z"].to_numpy(float)) if "etm_g4_z" in out.columns else False
    out["has_gradcpt_or_flanker"] = (
        np.isfinite(out["gradcpt_perf_z_age_sex"].to_numpy(float))
        | np.isfinite(out["flanker_efficiency_z_age_sex"].to_numpy(float))
    )
    out["has_gradcpt_and_flanker"] = (
        np.isfinite(out["gradcpt_perf_z_age_sex"].to_numpy(float))
        & np.isfinite(out["flanker_efficiency_z_age_sex"].to_numpy(float))
    )
    both = out["has_gradcpt_and_flanker"].to_numpy(bool)
    mean_raw = np.full(len(out), np.nan, dtype=float)
    mean_raw[both] = (
        out.loc[both, "gradcpt_perf_z_age_sex"].to_numpy(float)
        + out.loc[both, "flanker_efficiency_z_age_sex"].to_numpy(float)
    ) / 2.0
    mean_z = np.full(len(out), np.nan, dtype=float)
    if int(both.sum()) >= MIN_GROUP_N:
        m = float(np.mean(mean_raw[both]))
        sd = float(np.std(mean_raw[both], ddof=1))
        require(math.isfinite(sd) and sd > 0, "invalid GradCPT/Flanker mean SD")
        mean_z[both] = (mean_raw[both] - m) / sd
    out["gradcpt_flanker_mean_raw"] = mean_raw
    out["gradcpt_flanker_mean_z"] = mean_z
    if args.target == "gradcpt-flanker-mean":
        out["used_as_g4_finetune_label"] = out["has_gradcpt_and_flanker"]
    elif args.target == "all-etm-g4":
        out["used_as_g4_finetune_label"] = out["has_etm_g4"]
    else:
        out["used_as_g4_finetune_label"] = out["has_etm_g4"] & out["has_gradcpt_or_flanker"]
    out["used_as_finetune_label"] = out["used_as_g4_finetune_label"]
    if "final_model_train_allowed" not in out.columns:
        out["final_model_train_allowed"] = out["role"].eq("oof").astype(int)
    out["final_model_train_allowed"] = pd.to_numeric(
        out["final_model_train_allowed"], errors="coerce"
    ).fillna(0).astype(int)
    return out


def sex_map_from_base_covar(path: Path, iids: Iterable[str]) -> dict[str, int]:
    covar = pd.read_csv(path, sep="\t", dtype={"IID": str}, usecols=["IID", "sex_c"])
    covar["sex_c"] = pd.to_numeric(covar["sex_c"], errors="coerce")
    sex_map = {}
    for row in covar.itertuples(index=False):
        if not math.isfinite(float(row.sex_c)):
            continue
        sex = int(round(float(row.sex_c) + 0.5))
        if sex not in (0, 1):
            raise ValueError(f"Cannot derive sex_01 from sex_c={row.sex_c} for IID={row.IID}")
        sex_map[str(row.IID)] = sex
    missing = set(iids) - set(sex_map)
    require(not missing, f"base_covar.txt missing sex_c for {len(missing)} cohort IIDs")
    return sex_map


def rebuild_feature_matrix(args: argparse.Namespace, iids: list[str], saved_feature_columns: list[str]) -> np.ndarray:
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
        "nonresponse": Counter(),
        "question_encoding": {},
    }
    main_rows = ses.add_survey_feature_rows(str(args.main_survey), feature_state, metadata, bhp=False)
    bhp_rows = ses.add_survey_feature_rows(str(args.bhp_survey), feature_state, metadata, bhp=True)
    branch_rows = ses.apply_lifestyle_branch_recodes(feature_state)
    area_ses = ses.load_area_ses(str(args.area_ses))
    sex_map = sex_map_from_base_covar(args.ses_ea_dir / "base_covar.txt", iids)
    x, rebuilt_columns = ses.build_feature_matrix(iids, sex_map, area_ses, feature_state)
    if rebuilt_columns != saved_feature_columns:
        mismatches = []
        for idx, (got, expected) in enumerate(zip(rebuilt_columns, saved_feature_columns)):
            if got != expected:
                mismatches.append(f"{idx}: rebuilt={got!r} expected={expected!r}")
            if len(mismatches) >= 5:
                break
        raise SystemExit(
            "ERROR: rebuilt feature columns do not match xgboost_feature_columns.json. "
            f"rebuilt={len(rebuilt_columns)} expected={len(saved_feature_columns)} "
            f"first_mismatches={mismatches}"
        )
    print(
        f"Rebuilt feature matrix: rows={x.shape[0]} columns={x.shape[1]} "
        f"main_rows={main_rows} bhp_rows={bhp_rows} branch_rules={len(branch_rows)}",
        flush=True,
    )
    return x


def load_feature_contract(args: argparse.Namespace) -> tuple[list[str], str, pd.DataFrame]:
    feature_json = args.ses_ea_dir / "xgboost_feature_columns.json"
    manifest_path = args.ses_ea_dir / "xgboost_model_manifest.tsv"
    require(feature_json.exists(), f"missing {feature_json}")
    require(manifest_path.exists(), f"missing {manifest_path}")
    feature_columns = json.loads(feature_json.read_text())
    require(isinstance(feature_columns, list) and all(isinstance(x, str) for x in feature_columns), "invalid feature columns JSON")
    feature_hash = feature_columns_sha256(feature_columns)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    require("feature_columns_sha256" in manifest.columns, "model manifest missing feature_columns_sha256")
    manifest_hashes = set(manifest["feature_columns_sha256"].dropna().astype(str))
    require(manifest_hashes == {feature_hash}, f"feature hash mismatch: JSON={feature_hash}, manifest={sorted(manifest_hashes)}")
    return feature_columns, feature_hash, manifest


def load_booster(path: Path, feature_hash: str):
    import xgboost as xgb

    require(path.exists() and path.stat().st_size > 0, f"missing base booster: {path}")
    booster = xgb.Booster()
    booster.load_model(str(path))
    attrs = booster.attributes()
    require(attrs.get("feature_columns_sha256") == feature_hash, f"{path} has wrong feature_columns_sha256 attribute")
    return booster


def dmatrix(x: np.ndarray, feature_columns: list[str], label: np.ndarray | None = None):
    import xgboost as xgb

    return xgb.DMatrix(x, label=label, feature_names=feature_columns, missing=np.nan)


def stable_model_seed(base_seed: int, model_name: str) -> int:
    digest = hashlib.sha256(model_name.encode("utf-8")).hexdigest()
    return base_seed + int(digest[:8], 16) % 100000


def train_valid_split(indices: np.ndarray, valid_fraction: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    require(0.0 < valid_fraction < 0.5, "--valid-fraction must be >0 and <0.5")
    rng = np.random.default_rng(seed)
    shuffled = np.array(indices, dtype=np.int64)
    rng.shuffle(shuffled)
    valid_n = max(1, int(round(len(shuffled) * valid_fraction)))
    valid = shuffled[:valid_n]
    train = shuffled[valid_n:]
    require(len(train) > 0 and len(valid) > 0, "internal train/validation split is empty")
    return train, valid


def align_target_to_base_scale(
    target: np.ndarray,
    base_pred: np.ndarray,
    components: np.ndarray | None = None,
    component_names: list[str] | None = None,
) -> tuple[np.ndarray, dict[str, float]]:
    base_mean = float(np.mean(base_pred))
    base_sd = float(np.std(base_pred, ddof=1))
    require(math.isfinite(base_sd) and base_sd > 0, "invalid base raw-prediction SD for label alignment")
    if components is not None:
        require(components.ndim == 2 and components.shape[0] == len(target), "invalid component target matrix")
        names = component_names or [f"component_{i}" for i in range(components.shape[1])]
        require(len(names) == components.shape[1], "component_names length does not match components")
        aligned_components = []
        stats: dict[str, float] = {
            "base_train_pred_mean": base_mean,
            "base_train_pred_sd": base_sd,
        }
        for j, name in enumerate(names):
            comp = components[:, j].astype(float)
            comp_mean = float(np.mean(comp))
            comp_sd = float(np.std(comp, ddof=1))
            require(math.isfinite(comp_sd) and comp_sd > 0, f"invalid {name} SD for label alignment")
            aligned_component = base_mean + base_sd * ((comp - comp_mean) / comp_sd)
            aligned_components.append(aligned_component)
            stats[f"{name}_train_mean"] = comp_mean
            stats[f"{name}_train_sd"] = comp_sd
        component_mean = np.mean(np.column_stack(aligned_components), axis=1)
        pre_mean = float(np.mean(component_mean))
        pre_sd = float(np.std(component_mean, ddof=1))
        require(math.isfinite(pre_sd) and pre_sd > 0, "invalid combined component-label SD")
        aligned = base_mean + base_sd * ((component_mean - pre_mean) / pre_sd)
        stats.update(
            {
                "target_train_mean": float(np.mean(target)),
                "target_train_sd": float(np.std(target, ddof=1)),
                "component_mean_aligned_pre_zmatch_mean": pre_mean,
                "component_mean_aligned_pre_zmatch_sd": pre_sd,
                "aligned_label_mean": float(np.mean(aligned)),
                "aligned_label_sd": float(np.std(aligned, ddof=1)),
            }
        )
        return aligned.astype(np.float32), stats

    target_mean = float(np.mean(target))
    target_sd = float(np.std(target, ddof=1))
    require(math.isfinite(target_sd) and target_sd > 0, "invalid target SD for label alignment")
    aligned = base_mean + base_sd * ((target - target_mean) / target_sd)
    stats = {
        "base_train_pred_mean": base_mean,
        "base_train_pred_sd": base_sd,
        "target_train_mean": target_mean,
        "target_train_sd": target_sd,
        "etm_g4_train_mean": target_mean,
        "etm_g4_train_sd": target_sd,
        "aligned_label_mean": float(np.mean(aligned)),
        "aligned_label_sd": float(np.std(aligned, ddof=1)),
    }
    return aligned.astype(np.float32), stats


def xgb_params(args: argparse.Namespace) -> dict[str, object]:
    return {
        "objective": "reg:squarederror",
        "eta": args.eta,
        "max_depth": args.max_depth,
        "min_child_weight": args.min_child_weight,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "lambda": args.reg_lambda,
        "alpha": args.alpha,
        "tree_method": "hist",
        "eval_metric": "rmse",
        "seed": args.seed,
        "nthread": args.threads,
    }


def selected_rounds_from_eval(evals_result: dict[str, dict[str, list[float]]]) -> tuple[int, float]:
    valid = evals_result.get("valid", {}).get("rmse", [])
    require(len(valid) > 0, "XGBoost did not record validation RMSE")
    arr = np.asarray(valid, dtype=float)
    idx = int(np.nanargmin(arr))
    return idx + 1, float(arr[idx])


def fine_tune_one_model(
    *,
    args: argparse.Namespace,
    model_name: str,
    base_model_path: Path,
    x: np.ndarray,
    feature_columns: list[str],
    feature_hash: str,
    train_idx: np.ndarray,
    pred_idx: np.ndarray,
    target: np.ndarray,
    target_col: str,
    label_components: np.ndarray | None = None,
    label_component_names: list[str] | None = None,
    selected_rounds_override: int | None = None,
) -> tuple[np.ndarray, dict[str, object], object]:
    import xgboost as xgb

    require(len(train_idx) >= args.min_train_samples, f"{model_name} has too few fine-tuning labels: {len(train_idx)}")
    base = load_booster(base_model_path, feature_hash)
    base_pred_train = base.predict(dmatrix(x[train_idx, :], feature_columns))
    component_train = label_components[train_idx, :] if label_components is not None else None
    y_aligned, align_stats = align_target_to_base_scale(
        target[train_idx],
        base_pred_train,
        components=component_train,
        component_names=label_component_names,
    )

    if selected_rounds_override is None:
        internal_train, internal_valid = train_valid_split(
            np.arange(len(train_idx)), args.valid_fraction, stable_model_seed(args.seed, model_name)
        )
        evals_result: dict[str, dict[str, list[float]]] = {}
        tune = xgb.train(
            xgb_params(args),
            dmatrix(x[train_idx[internal_train], :], feature_columns, label=y_aligned[internal_train]),
            num_boost_round=args.max_rounds,
            evals=[
                (dmatrix(x[train_idx[internal_train], :], feature_columns, label=y_aligned[internal_train]), "train"),
                (dmatrix(x[train_idx[internal_valid], :], feature_columns, label=y_aligned[internal_valid]), "valid"),
            ],
            early_stopping_rounds=args.early_stopping_rounds,
            evals_result=evals_result,
            xgb_model=base,
            verbose_eval=False,
        )
        selected_rounds, best_valid_rmse = selected_rounds_from_eval(evals_result)
        del tune
    else:
        selected_rounds = int(selected_rounds_override)
        best_valid_rmse = math.nan

    base_full = load_booster(base_model_path, feature_hash)
    finetuned = xgb.train(
        xgb_params(args),
        dmatrix(x[train_idx, :], feature_columns, label=y_aligned),
        num_boost_round=selected_rounds,
        xgb_model=base_full,
        verbose_eval=False,
    )
    finetuned.set_attr(
        base_model_file=str(base_model_path),
        fine_tune_target=target_col,
        target_mode=args.target,
        aligned_label="base_raw_scale_zmatched_target",
        feature_columns_sha256=feature_hash,
        fine_tune_samples=str(len(train_idx)),
        prediction_samples=str(len(pred_idx)),
        selected_rounds=str(selected_rounds),
        seed=str(args.seed),
        eta=str(args.eta),
        max_depth=str(args.max_depth),
        min_child_weight=str(args.min_child_weight),
        reg_lambda=str(args.reg_lambda),
    )
    pred = finetuned.predict(dmatrix(x[pred_idx, :], feature_columns))
    row: dict[str, object] = {
        "model_name": model_name,
        "base_model_file": str(base_model_path),
        "fine_tuned_model_file": f"xgboost_models/{model_name}_g4_finetuned.json",
        "target_mode": args.target,
        "fine_tune_samples": len(train_idx),
        "prediction_samples": len(pred_idx),
        "selected_rounds": selected_rounds,
        "best_valid_rmse_aligned_scale": best_valid_rmse,
        "feature_columns_sha256": feature_hash,
    }
    row.update(align_stats)
    return pred, row, finetuned


def build_run_signature(args: argparse.Namespace, feature_hash: str) -> tuple[str, pd.DataFrame]:
    paths = {
        "all_scores_file": args.all_scores_file,
        "task_score_file": args.task_score_file,
        "feature_columns_json": args.ses_ea_dir / "xgboost_feature_columns.json",
        "model_manifest": args.ses_ea_dir / "xgboost_model_manifest.tsv",
        "base_covar": args.ses_ea_dir / "base_covar.txt",
        "covar": args.ses_ea_dir / "covar.txt",
        "ea_query": args.ea_query,
        "main_survey": args.main_survey,
        "bhp_survey": args.bhp_survey,
        "area_ses": args.area_ses,
        "metadata": args.metadata,
        "fine_tune_script": Path(__file__).resolve(),
        "setup_script": Path(ses.__file__).resolve(),
    }
    if target_uses_etm_g(args.target):
        require(args.etm_g_file is not None, f"{args.target} requires --etm-g-file")
        paths["etm_g_file"] = args.etm_g_file
    for name in MODEL_NAMES:
        paths[f"base_model_{name}"] = args.ses_ea_dir / "xgboost_models" / f"{name}.json"
    payload = {
        "paths": {name: file_fingerprint(path) for name, path in paths.items()},
        "target_mode": args.target,
        "target_column": target_column_for_mode(args.target),
        "score_prefix": SCORE_PREFIX_BY_TARGET[args.target],
        "feature_columns_sha256": feature_hash,
        "params": {
            "eta": args.eta,
            "max_depth": args.max_depth,
            "min_child_weight": args.min_child_weight,
            "lambda": args.reg_lambda,
            "alpha": args.alpha,
            "max_rounds": args.max_rounds,
            "early_stopping_rounds": args.early_stopping_rounds,
            "valid_fraction": args.valid_fraction,
            "min_train_samples": args.min_train_samples,
            "seed": args.seed,
            "threads": args.threads,
            "write_regenie_inputs": not args.no_write_regenie_inputs,
        },
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    signature = hashlib.sha256(blob.encode("utf-8")).hexdigest()
    rows = [{"parameter": "run_signature", "value": signature}]
    rows += [{"parameter": key, "value": json.dumps(value, sort_keys=True)} for key, value in payload.items()]
    return signature, pd.DataFrame(rows)


def existing_signature(output_dir: Path) -> str | None:
    path = output_dir / "g4_finetuned_runtime_manifest.json"
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return None
    value = payload.get("run_signature")
    return str(value) if value else None


def target_column_for_mode(target_mode: str) -> str:
    if target_mode in ("strong-task-g4", "all-etm-g4"):
        return "etm_g4_z"
    if target_mode == "gradcpt-flanker-mean":
        return "gradcpt_flanker_mean_z"
    raise ValueError(f"Unsupported target mode: {target_mode}")


def label_component_columns_for_mode(target_mode: str) -> list[str] | None:
    if target_mode == "gradcpt-flanker-mean":
        return ["gradcpt_perf_z_age_sex", "flanker_efficiency_z_age_sex"]
    return None


def expected_outputs(output_dir: Path, score_prefix: str, write_regenie_inputs: bool) -> list[Path]:
    paths = [
        output_dir / "g4_finetuned_ea_proxy_scores_wide.tsv",
        output_dir / "g4_finetuned_ea_proxy_oof_scores.tsv",
        output_dir / "g4_finetuned_ea_proxy_applied_scores.tsv",
        output_dir / "g4_finetuned_model_manifest.tsv",
        output_dir / "g4_finetuned_runtime_manifest.json",
    ] + [output_dir / "xgboost_models" / f"{name}_g4_finetuned.json" for name in MODEL_NAMES]
    if write_regenie_inputs:
        paths.extend(
            [
                output_dir / "phen.txt",
                output_dir / f"{score_prefix}.phen.txt",
                output_dir / "base_covar.txt",
                output_dir / "covar.txt",
                output_dir / "training_iids.txt",
            ]
        )
    return paths


def should_skip(args: argparse.Namespace, signature: str, score_prefix: str) -> bool:
    if args.force:
        return False
    paths = expected_outputs(args.output_dir, score_prefix, not args.no_write_regenie_inputs)
    if not any(path.exists() for path in paths):
        return False
    saved = existing_signature(args.output_dir)
    if saved == signature and all(path.exists() and path.stat().st_size > 0 for path in paths):
        print(f"Existing G4 fine-tuned outputs match inputs/parameters; skipping: {args.output_dir}", flush=True)
        return True
    raise SystemExit(
        "ERROR: G4 fine-tuned outputs already exist but do not match current inputs/parameters. "
        "Rerun with --force to overwrite."
    )


def correlation_rows(df: pd.DataFrame, groups: dict[str, pd.Series], score_cols: list[str], target_cols: list[str]) -> pd.DataFrame:
    rows = []
    for group_name, mask in groups.items():
        mask_arr = mask.to_numpy(bool) if isinstance(mask, pd.Series) else np.asarray(mask, dtype=bool)
        for score_col in score_cols:
            for target_col in target_cols:
                if score_col not in df.columns or target_col not in df.columns:
                    continue
                p, n = pearson(df.loc[mask_arr, score_col], df.loc[mask_arr, target_col])
                s, _ = spearman(df.loc[mask_arr, score_col], df.loc[mask_arr, target_col])
                rows.append(
                    {
                        "group": group_name,
                        "score": score_col,
                        "target": target_col,
                        "n": n,
                        "pearson": p,
                        "spearman": s,
                    }
                )
    return pd.DataFrame(rows)


def distribution_rows(df: pd.DataFrame, groups: dict[str, pd.Series], score_cols: list[str]) -> pd.DataFrame:
    rows = []
    for group_name, mask in groups.items():
        mask_arr = mask.to_numpy(bool) if isinstance(mask, pd.Series) else np.asarray(mask, dtype=bool)
        for score_col in score_cols:
            if score_col not in df.columns:
                continue
            row = {"group": group_name, "score": score_col}
            row.update(summarize(df.loc[mask_arr, score_col]))
            rows.append(row)
    return pd.DataFrame(rows)


def add_fold_safe_linear_calibration(
    df: pd.DataFrame,
    *,
    args: argparse.Namespace,
    fine_tuned_z_col: str,
    target_col: str,
    diag_dir: Path,
) -> tuple[str, str, pd.DataFrame]:
    require(args.target in CALIBRATED_PREFIX_BY_TARGET, "linear calibration is only defined for supported targets")
    prefix = CALIBRATED_PREFIX_BY_TARGET[args.target]
    raw_col = f"{prefix}_raw"
    z_col = f"{prefix}_z"
    feature_cols = ["teacher_z", "ses_ea_proxy_z", fine_tuned_z_col]
    missing = [
        col for col in feature_cols + [
            target_col, "role", "fold_id", "used_as_finetune_label", "final_model_train_allowed",
        ]
        if col not in df.columns
    ]
    require(not missing, f"missing columns for linear calibration: {missing}")

    feature_mat = df[feature_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    target = pd.to_numeric(df[target_col], errors="coerce").to_numpy(float)
    features_finite = np.all(np.isfinite(feature_mat), axis=1)
    labels_finite = features_finite & np.isfinite(target) & df["used_as_finetune_label"].to_numpy(bool)
    role = df["role"].astype(str)
    fold_str = df["fold_id"].astype(str)
    oof = role.eq("oof").to_numpy(bool)
    applied = role.eq("applied").to_numpy(bool)
    final_model_train_allowed = df["final_model_train_allowed"].astype(bool).to_numpy(bool)
    pred = np.full(len(df), np.nan, dtype=float)
    coef_rows: list[dict[str, object]] = []

    scale_rows: list[dict[str, object]] = []
    scale_groups = {
        "linear_calibration_target_labeled_all": labels_finite,
        "linear_calibration_target_labeled_oof": labels_finite & oof,
        "linear_calibration_target_labeled_final_model_allowed": labels_finite & oof & final_model_train_allowed,
        "linear_calibration_target_labeled_applied": labels_finite & applied,
        "finite_feature_full_cohort": features_finite,
    }
    scale_values = {col: feature_mat[:, idx] for idx, col in enumerate(feature_cols)}
    scale_values[target_col] = target
    for group_name, group_mask in scale_groups.items():
        for col, values in scale_values.items():
            vals = values[group_mask & np.isfinite(values)]
            row = {
                "group": group_name,
                "variable": col,
                "n": int(len(vals)),
                "mean": float(np.mean(vals)) if len(vals) else math.nan,
                "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else math.nan,
                "min": float(np.min(vals)) if len(vals) else math.nan,
                "p50": float(np.quantile(vals, 0.50)) if len(vals) else math.nan,
                "max": float(np.max(vals)) if len(vals) else math.nan,
            }
            scale_rows.append(row)
    pd.DataFrame(scale_rows).to_csv(diag_dir / f"{prefix}_linear_calibration_feature_scales.tsv", sep="\t", index=False)

    def fit_predict(train_mask: np.ndarray, pred_mask: np.ndarray, fit_name: str, predict_group: str) -> None:
        require(int(train_mask.sum()) >= args.min_train_samples, f"too few labels for {fit_name} linear calibration")
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

    for fold in range(5):
        train_mask = oof & ~fold_str.eq(str(fold)).to_numpy(bool) & labels_finite
        pred_mask = oof & fold_str.eq(str(fold)).to_numpy(bool) & features_finite
        fit_predict(train_mask, pred_mask, f"oof_fold_{fold}_train_other_folds", f"oof_fold_{fold}")

    train_mask = oof & labels_finite & final_model_train_allowed
    pred_mask = applied & features_finite
    fit_predict(train_mask, pred_mask, "applied_model_train_all_oof", "applied")

    require(np.all(np.isfinite(pred[features_finite])), "some finite-feature rows did not receive a calibrated prediction")
    mean = float(np.nanmean(pred))
    sd = float(np.nanstd(pred, ddof=1))
    require(math.isfinite(sd) and sd > 0, "invalid calibrated prediction SD")
    df[raw_col] = pred
    df[z_col] = (pred - mean) / sd
    coef_df = pd.DataFrame(coef_rows)
    coef_df["calibrated_raw_mean_full_cohort"] = mean
    coef_df["calibrated_raw_sd_full_cohort"] = sd
    coef_df.to_csv(diag_dir / f"{prefix}_linear_calibration_coefficients.tsv", sep="\t", index=False)
    return raw_col, z_col, coef_df


def target_overlap_counts(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "role",
        "fold_id",
        "has_etm_g4",
        "has_gradcpt_or_flanker",
        "has_gradcpt_and_flanker",
        "used_as_finetune_label",
        "final_model_train_allowed",
        "used_as_g4_finetune_label",
    ]
    if "task_pattern_four_domain" in df.columns:
        group_cols.append("task_pattern_four_domain")
    out = df.groupby(group_cols, dropna=False).size().reset_index(name="n")
    return out


def feature_importance(base_final, finetuned_final) -> pd.DataFrame:
    base_gain = base_final.get_score(importance_type="gain")
    base_cover = base_final.get_score(importance_type="cover")
    fine_gain = finetuned_final.get_score(importance_type="gain")
    fine_cover = finetuned_final.get_score(importance_type="cover")
    features = sorted(set(base_gain) | set(base_cover) | set(fine_gain) | set(fine_cover))
    rows = []
    for feature in features:
        bg = float(base_gain.get(feature, 0.0))
        fg = float(fine_gain.get(feature, 0.0))
        rows.append(
            {
                "feature": feature,
                "base_gain": bg,
                "finetuned_gain": fg,
                "base_cover": float(base_cover.get(feature, 0.0)),
                "finetuned_cover": float(fine_cover.get(feature, 0.0)),
                "gain_delta": fg - bg,
                "gain_ratio": fg / bg if bg > 0 else math.nan,
            }
        )
    rows.sort(key=lambda r: (-float(r["finetuned_gain"]), r["feature"]))
    return pd.DataFrame(rows)


def write_regenie_inputs(
    args: argparse.Namespace,
    output_dir: Path,
    score_prefix: str,
    wide: pd.DataFrame,
    *,
    primary_z_col: str | None = None,
    extra_z_cols: list[str] | None = None,
) -> None:
    z_col = primary_z_col or f"{score_prefix}_z"
    phen_path = output_dir / "phen.txt"
    with phen_path.open("w") as handle:
        handle.write(f"FID\tIID\t{z_col}\n")
        for row in wide[["FID", "IID", z_col]].itertuples(index=False):
            handle.write(f"{row.FID}\t{row.IID}\t{float(getattr(row, z_col)):.10g}\n")
    primary_prefix = z_col[:-2] if z_col.endswith("_z") else z_col
    shutil.copyfile(phen_path, output_dir / f"{primary_prefix}.phen.txt")
    for extra_col in extra_z_cols or []:
        if extra_col == z_col:
            continue
        extra_prefix = extra_col[:-2] if extra_col.endswith("_z") else extra_col
        with (output_dir / f"{extra_prefix}.phen.txt").open("w") as handle:
            handle.write(f"FID\tIID\t{extra_col}\n")
            for row in wide[["FID", "IID", extra_col]].itertuples(index=False):
                handle.write(f"{row.FID}\t{row.IID}\t{float(getattr(row, extra_col)):.10g}\n")
    shutil.copyfile(args.ses_ea_dir / "base_covar.txt", output_dir / "base_covar.txt")
    shutil.copyfile(args.ses_ea_dir / "covar.txt", output_dir / "covar.txt")
    with (output_dir / "training_iids.txt").open("w") as handle:
        for row in wide[["FID", "IID"]].itertuples(index=False):
            handle.write(f"{row.FID} {row.IID}\n")


def copy_outputs(output_dir: Path, workspace_output_dir: Path | None, workspace_scratch_dir: Path | None, stage_aggregate: bool) -> None:
    workspace_mount = Path(os.environ.get("WORKSPACE_BUCKET_MOUNT", "/home/jupyter/workspace/workspace-bucket"))
    workspace_uri = os.environ.get("WORKSPACE_BUCKET_URI", "")
    google_project = os.environ.get("GOOGLE_PROJECT", "")

    def workspace_path_to_uri(path: Path) -> str | None:
        if not workspace_uri:
            return None
        try:
            rel = path.resolve().relative_to(workspace_mount.resolve())
        except ValueError:
            return None
        return workspace_uri.rstrip("/") + "/" + str(rel)

    def copy_file_retry(src: Path, dest: Path, attempts: int = 4) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        last_error: Exception | None = None
        for attempt in range(1, attempts + 1):
            try:
                dest_uri = workspace_path_to_uri(dest)
                if dest_uri is not None:
                    cmd = ["gcloud", "storage", "cp", str(src), dest_uri]
                    if google_project:
                        cmd.append(f"--billing-project={google_project}")
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                else:
                    if dest.exists():
                        dest.unlink()
                    shutil.copyfile(src, dest)
                return
            except (OSError, subprocess.CalledProcessError) as exc:
                last_error = exc
                detail = exc.stderr.strip() if isinstance(exc, subprocess.CalledProcessError) and exc.stderr else str(exc)
                print(f"Copy attempt {attempt}/{attempts} failed for {dest}: {detail}", flush=True)
                time.sleep(2 * attempt)
        raise last_error if last_error is not None else RuntimeError(f"failed to copy {src} to {dest}")

    if workspace_output_dir is not None:
        workspace_output_dir.mkdir(parents=True, exist_ok=True)
        for path in output_dir.iterdir():
            if path.name == "diagnostics":
                continue
            dest = workspace_output_dir / path.name
            if path.is_dir():
                if dest.exists():
                    shutil.rmtree(dest)
                dest.mkdir(parents=True, exist_ok=True)
                for child in sorted(path.iterdir()):
                    if child.is_file():
                        copy_file_retry(child, dest / child.name)
            else:
                copy_file_retry(path, dest)
        print(f"Copied G4 fine-tuned outputs to {workspace_output_dir}", flush=True)
    if stage_aggregate:
        require(workspace_scratch_dir is not None, "--stage-aggregate requires --workspace-scrap-dir")
        workspace_scratch_dir.mkdir(parents=True, exist_ok=True)
        diag_dir = output_dir / "diagnostics"
        for path in diag_dir.glob("*.tsv"):
            copy_file_retry(path, workspace_scratch_dir / path.name)
        print(f"Staged aggregate diagnostics to {workspace_scratch_dir}", flush=True)


def main() -> None:
    start = time.time()
    args = parse_args()
    require(args.max_rounds >= 1, "--max-rounds must be positive")
    require(args.early_stopping_rounds >= 1, "--early-stopping-rounds must be positive")
    score_prefix = SCORE_PREFIX_BY_TARGET[args.target]
    target_col = target_column_for_mode(args.target)
    component_cols = label_component_columns_for_mode(args.target)

    feature_columns, feature_hash, base_manifest = load_feature_contract(args)
    signature, signature_rows = build_run_signature(args, feature_hash)
    if should_skip(args, signature, score_prefix):
        return
    if args.force and args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    model_dir = args.output_dir / "xgboost_models"
    diag_dir = args.output_dir / "diagnostics"
    model_dir.mkdir(parents=True, exist_ok=True)
    diag_dir.mkdir(parents=True, exist_ok=True)

    print("=== G4-finetuned SES-EA proxy ===", flush=True)
    print(f"Target mode: {args.target}", flush=True)
    print(f"Score prefix: {score_prefix}", flush=True)
    print(f"Feature columns: {len(feature_columns)} hash={feature_hash}", flush=True)

    all_scores, covars = load_base_tables(args)
    df = load_etm_tables(args, all_scores)
    df = df.merge(covars, on=["FID", "IID"], how="left", validate="one_to_one", suffixes=("", "_covar"))
    x = rebuild_feature_matrix(args, df["IID"].astype(str).tolist(), feature_columns)
    require(x.shape[0] == len(df), "feature matrix row count does not match score table")
    require(feature_columns_sha256(feature_columns) == feature_hash, "feature hash changed after rebuild")

    import xgboost as xgb

    require(target_col in df.columns, f"target column {target_col} was not built")
    target = df[target_col].to_numpy(float)
    target_available = pd.Series(np.isfinite(target), index=df.index)
    label_components = None
    if component_cols is not None:
        missing_components = [col for col in component_cols if col not in df.columns]
        require(not missing_components, f"missing component target columns: {missing_components}")
        label_components = df[component_cols].to_numpy(float)
    pred_raw = np.full(len(df), np.nan, dtype=float)
    pred_model_name = np.full(len(df), "", dtype=object)
    model_rows = []
    round_rows = []
    saved_boosters = {}
    label_mask = df["used_as_finetune_label"].to_numpy(bool)

    for fold in range(5):
        train_mask = (df["role"].eq("oof") & ~df["fold_id"].astype(str).eq(str(fold))).to_numpy(bool) & label_mask
        pred_mask = (df["role"].eq("oof") & df["fold_id"].astype(str).eq(str(fold))).to_numpy(bool)
        train_idx = np.flatnonzero(train_mask)
        pred_idx = np.flatnonzero(pred_mask)
        base_path = args.ses_ea_dir / "xgboost_models" / f"fold_{fold}.json"
        pred, row, booster = fine_tune_one_model(
            args=args,
            model_name=f"fold_{fold}",
            base_model_path=base_path,
            x=x,
            feature_columns=feature_columns,
            feature_hash=feature_hash,
            train_idx=train_idx,
            pred_idx=pred_idx,
            target=target,
            target_col=target_col,
            label_components=label_components,
            label_component_names=component_cols,
        )
        pred_raw[pred_idx] = pred
        pred_model_name[pred_idx] = f"fold_{fold}_g4_finetuned"
        model_path = model_dir / f"fold_{fold}_g4_finetuned.json"
        booster.save_model(str(model_path))
        model_rows.append(row)
        round_rows.append(row.copy())
        saved_boosters[f"fold_{fold}"] = booster
        print(
            f"fold_{fold}: fine_tune_n={len(train_idx)} predict_n={len(pred_idx)} "
            f"selected_rounds={row['selected_rounds']}",
            flush=True,
        )

    selected_final_rounds = int(np.median([int(r["selected_rounds"]) for r in round_rows]))
    final_model_train_allowed = df["final_model_train_allowed"].astype(bool).to_numpy(bool)
    final_train_mask = df["role"].eq("oof").to_numpy(bool) & label_mask & final_model_train_allowed
    final_pred_mask = df["role"].eq("applied").to_numpy(bool)
    final_pred, final_row, final_booster = fine_tune_one_model(
        args=args,
        model_name="final_model",
        base_model_path=args.ses_ea_dir / "xgboost_models" / "final_model.json",
        x=x,
        feature_columns=feature_columns,
        feature_hash=feature_hash,
        train_idx=np.flatnonzero(final_train_mask),
        pred_idx=np.flatnonzero(final_pred_mask),
        target=target,
        target_col=target_col,
        label_components=label_components,
        label_component_names=component_cols,
        selected_rounds_override=selected_final_rounds,
    )
    pred_raw[final_pred_mask] = final_pred
    pred_model_name[final_pred_mask] = "final_model_g4_finetuned"
    final_booster.save_model(str(model_dir / "final_model_g4_finetuned.json"))
    model_rows.append(final_row)
    round_rows.append(final_row.copy())
    saved_boosters["final_model"] = final_booster
    print(
        f"final_model: fine_tune_n={int(final_train_mask.sum())} predict_n={int(final_pred_mask.sum())} "
        f"selected_rounds={selected_final_rounds}",
        flush=True,
    )

    require(np.all(np.isfinite(pred_raw)), "not every cohort row received a fine-tuned prediction")
    pred_mean = float(np.mean(pred_raw))
    pred_sd = float(np.std(pred_raw, ddof=1))
    require(math.isfinite(pred_sd) and pred_sd > 0, "invalid fine-tuned prediction SD")
    z_col = f"{score_prefix}_z"
    raw_col = f"{score_prefix}_raw"
    df[raw_col] = pred_raw
    df[z_col] = (pred_raw - pred_mean) / pred_sd
    df["prediction_model_name"] = pred_model_name
    primary_z_col = z_col
    extra_phenotype_z_cols: list[str] = []
    calibrated_raw_col = None
    calibrated_z_col = None
    calibration_rows = None
    if args.target in CALIBRATED_PREFIX_BY_TARGET:
        calibrated_raw_col, calibrated_z_col, calibration_rows = add_fold_safe_linear_calibration(
            df,
            args=args,
            fine_tuned_z_col=z_col,
            target_col=target_col,
            diag_dir=diag_dir,
        )
        primary_z_col = calibrated_z_col
        extra_phenotype_z_cols.append(z_col)

    base_score_col = "ses_ea_proxy_z"
    output_cols = [
        "FID",
        "IID",
        "role",
        "fold_id",
        "ea_years",
        "teacher_z",
        "ses_ea_proxy_z",
        "score_raw",
        "gradcpt_flanker_mean_raw",
        "gradcpt_flanker_mean_z",
        "has_etm_g4",
        "has_gradcpt_or_flanker",
        "has_gradcpt_and_flanker",
        "used_as_g4_finetune_label",
        "used_as_finetune_label",
        "final_model_train_allowed",
        raw_col,
        z_col,
        "prediction_model_name",
    ]
    for col in ["etm_g4_z", "etm_g_z", "n_tasks_observed_four_domain", "task_pattern_four_domain"]:
        if col in df.columns and col not in output_cols:
            output_cols.append(col)
    for col in [calibrated_raw_col, calibrated_z_col]:
        if col and col not in output_cols:
            output_cols.append(col)
    for col in TASK_COLS:
        if col in df.columns and col not in output_cols:
            output_cols.append(col)
    wide = df.loc[:, output_cols].copy()
    wide.to_csv(args.output_dir / "g4_finetuned_ea_proxy_scores_wide.tsv", sep="\t", index=False)
    wide.loc[wide["role"].eq("oof")].to_csv(args.output_dir / "g4_finetuned_ea_proxy_oof_scores.tsv", sep="\t", index=False)
    wide.loc[wide["role"].eq("applied")].to_csv(args.output_dir / "g4_finetuned_ea_proxy_applied_scores.tsv", sep="\t", index=False)
    if not args.no_write_regenie_inputs:
        write_regenie_inputs(
            args,
            args.output_dir,
            score_prefix,
            wide,
            primary_z_col=primary_z_col,
            extra_z_cols=extra_phenotype_z_cols,
        )

    model_manifest = pd.DataFrame(model_rows)
    model_manifest.to_csv(args.output_dir / "g4_finetuned_model_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(round_rows).to_csv(diag_dir / "g4_finetuned_round_selection.tsv", sep="\t", index=False)

    score_cols = [base_score_col, z_col]
    if calibrated_z_col:
        score_cols.append(calibrated_z_col)
    target_cols = [
        target_col,
        "etm_g4_z",
        "gradcpt_flanker_mean_z",
        "teacher_z",
        "ea_years",
        "dd_patience_z_age_sex",
        "gradcpt_perf_z_age_sex",
        "flanker_efficiency_z_age_sex",
        "emorecog_perf_z_age_sex",
        "yob_c",
        "sex_c",
    ] + [f"PC{i}_AVG" for i in range(1, 11)]
    target_cols = list(dict.fromkeys(target_cols))
    groups = {
        "oof_target_overall": df["role"].eq("oof") & target_available,
        "applied_target_holdout": df["role"].eq("applied") & target_available,
        "combined_target_holdout": target_available,
        "combined_all_scored": pd.Series(True, index=df.index),
        "primary_finetune_label_set_diagnostic_only": df["used_as_finetune_label"],
        "final_model_finetune_label_set_diagnostic_only": df["used_as_finetune_label"] & df["final_model_train_allowed"].astype(bool),
    }
    for fold in range(5):
        groups[f"fold_{fold}_target_heldout"] = df["role"].eq("oof") & df["fold_id"].astype(str).eq(str(fold)) & target_available
    corr = correlation_rows(df, groups, score_cols, [c for c in target_cols if c in df.columns])
    corr.to_csv(diag_dir / "g4_finetuned_before_after_correlations.tsv", sep="\t", index=False)
    covar_targets = [c for c in ["yob_c", "sex_c"] + [f"PC{i}_AVG" for i in range(1, 11)] if c in df.columns]
    correlation_rows(df, groups, score_cols, covar_targets).to_csv(
        diag_dir / "g4_finetuned_covariate_correlations.tsv", sep="\t", index=False
    )
    correlation_rows(
        df,
        {"role_oof": df["role"].eq("oof"), "role_applied": df["role"].eq("applied"), "combined": pd.Series(True, index=df.index)},
        score_cols,
        [c for c in target_cols if c in df.columns],
    ).to_csv(diag_dir / "g4_finetuned_group_correlations.tsv", sep="\t", index=False)

    dist_groups = {
        "all_scored": pd.Series(True, index=df.index),
        "role_oof": df["role"].eq("oof"),
        "role_applied": df["role"].eq("applied"),
        "has_target": target_available,
        "no_target": ~target_available,
        "has_etm_g4": df["has_etm_g4"],
        "no_etm_g4": ~df["has_etm_g4"],
        "primary_finetune_target": df["used_as_finetune_label"],
        "etm_but_no_gradcpt_or_flanker": df["has_etm_g4"] & ~df["has_gradcpt_or_flanker"],
        "non_etm_samples": ~df["has_etm_g4"],
    }
    distribution_raw_cols = [raw_col] + ([calibrated_raw_col] if calibrated_raw_col else [])
    distribution_rows(df, dist_groups, score_cols + distribution_raw_cols).to_csv(
        diag_dir / "g4_finetuned_score_distributions.tsv", sep="\t", index=False
    )
    distribution_rows(
        df,
        {"target_labeled": target_available, "no_target": ~target_available},
        score_cols + distribution_raw_cols,
    ).to_csv(diag_dir / "g4_finetuned_no_target_distribution_shift.tsv", sep="\t", index=False)

    target_overlap_counts(df).to_csv(diag_dir / "g4_finetuned_target_overlap_counts.tsv", sep="\t", index=False)
    if "task_pattern_four_domain" in df.columns:
        pattern_groups = {
            f"pattern_{pattern}": df["task_pattern_four_domain"].fillna("NONE").eq(pattern)
            for pattern in sorted(df["task_pattern_four_domain"].fillna("NONE").unique())
        }
        correlation_rows(df, pattern_groups, score_cols, ["etm_g4_z", "teacher_z", "ea_years"]).to_csv(
            diag_dir / "g4_finetuned_validation_by_task_pattern.tsv", sep="\t", index=False
        )
    else:
        pd.DataFrame().to_csv(diag_dir / "g4_finetuned_validation_by_task_pattern.tsv", sep="\t", index=False)

    base_final = load_booster(args.ses_ea_dir / "xgboost_models" / "final_model.json", feature_hash)
    feature_importance(base_final, final_booster).to_csv(
        diag_dir / "g4_finetuned_feature_importance.tsv", sep="\t", index=False
    )
    model_manifest.to_csv(diag_dir / "g4_finetuned_model_manifest.tsv", sep="\t", index=False)

    runtime_manifest = {
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "run_signature": signature,
        "target_mode": args.target,
        "score_prefix": score_prefix,
        "score_raw_mean_full_cohort": pred_mean,
        "score_raw_sd_full_cohort": pred_sd,
        "feature_columns_sha256": feature_hash,
        "feature_columns": len(feature_columns),
        "xgboost_version": xgb.__version__,
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "numpy_version": np.__version__,
        "params": xgb_params(args),
        "primary_phenotype_z_col": primary_z_col,
        "calibrated_raw_col": calibrated_raw_col,
        "calibrated_z_col": calibrated_z_col,
        "linear_calibration_rows": 0 if calibration_rows is None else int(len(calibration_rows)),
        "base_manifest_rows": int(len(base_manifest)),
        "fine_tuned_models": len(model_rows),
        "cohort_rows": int(len(df)),
        "has_etm_g4": int(df["has_etm_g4"].sum()),
        "target_column": target_col,
        "used_as_finetune_label": int(df["used_as_finetune_label"].sum()),
        "used_as_g4_finetune_label": int(df["used_as_g4_finetune_label"].sum()),
        "final_model_train_allowed": int(df["final_model_train_allowed"].sum()),
        "final_model_used_as_finetune_label": int((df["used_as_finetune_label"] & df["final_model_train_allowed"].astype(bool)).sum()),
    }
    (args.output_dir / "g4_finetuned_runtime_manifest.json").write_text(json.dumps(runtime_manifest, indent=2, sort_keys=True) + "\n")
    pd.concat(
        [
            signature_rows,
            pd.DataFrame(
                [
                    {"parameter": "target_mode", "value": args.target},
                    {"parameter": "score_prefix", "value": score_prefix},
                    {"parameter": "cohort_rows", "value": len(df)},
                    {"parameter": "has_etm_g4", "value": int(df["has_etm_g4"].sum())},
                    {"parameter": "target_column", "value": target_col},
                    {"parameter": "used_as_finetune_label", "value": int(df["used_as_finetune_label"].sum())},
                    {"parameter": "used_as_g4_finetune_label", "value": int(df["used_as_g4_finetune_label"].sum())},
                    {"parameter": "final_model_train_allowed", "value": int(df["final_model_train_allowed"].sum())},
                    {
                        "parameter": "final_model_used_as_finetune_label",
                        "value": int((df["used_as_finetune_label"] & df["final_model_train_allowed"].astype(bool)).sum()),
                    },
                    {"parameter": "score_raw_mean_full_cohort", "value": pred_mean},
                    {"parameter": "score_raw_sd_full_cohort", "value": pred_sd},
                    {"parameter": "primary_phenotype_z_col", "value": primary_z_col},
                    {"parameter": "calibrated_raw_col", "value": calibrated_raw_col or ""},
                    {"parameter": "calibrated_z_col", "value": calibrated_z_col or ""},
                ]
            ),
        ],
        ignore_index=True,
    ).to_csv(args.output_dir / "g4_finetuned_params.tsv", sep="\t", index=False)

    log_lines = [
        "=== G4-finetuned SES-EA proxy ===",
        f"target_mode={args.target}",
        f"cohort_rows={len(df)}",
        f"has_etm_g4={int(df['has_etm_g4'].sum())}",
        f"target_column={target_col}",
        f"used_as_finetune_label={int(df['used_as_finetune_label'].sum())}",
        f"final_model_train_allowed={int(df['final_model_train_allowed'].sum())}",
        f"final_model_used_as_finetune_label={int((df['used_as_finetune_label'] & df['final_model_train_allowed'].astype(bool)).sum())}",
        f"feature_columns_sha256={feature_hash}",
        f"score_raw_mean_full_cohort={pred_mean:.10g}",
        f"score_raw_sd_full_cohort={pred_sd:.10g}",
        f"primary_phenotype_z_col={primary_z_col}",
        f"calibrated_z_col={calibrated_z_col or ''}",
    ]
    (args.output_dir / "g4_finetuned_ea_proxy_log.txt").write_text("\n".join(log_lines) + "\n")
    copy_outputs(args.output_dir, args.workspace_output_dir, args.workspace_scrap_dir, args.stage_aggregate)
    print("\nBefore/after overall correlations:", flush=True)
    keep = corr["group"].isin(["combined_target_holdout", "applied_target_holdout"]) & corr["target"].eq(target_col)
    print(corr.loc[keep].to_string(index=False), flush=True)
    print("\n=== G4 fine-tuning complete ===", flush=True)


if __name__ == "__main__":
    main()
