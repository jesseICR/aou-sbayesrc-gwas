"""Classify European samples from this pipeline's ADMIXTURE K=6 output.

This is the Step 6 fallback used when AoU hard ancestry predictions are
available but AoU RYE admixture fractions are not. It writes the same downstream
European keep-list consumed by PCA/GWAS steps, plus aggregate hard-call overlap
summaries when AoU ancestry predictions are supplied.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


OURS_COMPONENTS = [
    ("European", "European"),
    ("East_Asian", "East_Asian"),
    ("American", "American"),
    ("African", "African"),
    ("South_Asian", "South_Asian"),
    ("Oceanian", "Oceanian"),
]
PRED_ORDER = ["eur", "afr", "amr", "eas", "sas", "mid", "missing"]
OVERLAP_ORDER = ["both", "aou_eur_only", "ours_eur_only", "neither"]


def require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")


def read_ours(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"FID": str, "IID": str})
    require_columns(
        frame,
        {"FID", "IID", "European", "East_Asian", "American", "African", "South_Asian", "Oceanian"},
        path,
    )
    for column, _label in OURS_COMPONENTS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["FID"] = frame["FID"].astype(str)
    frame["IID"] = frame["IID"].astype(str)
    return frame


def read_aou_pred(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype={"research_id": str, "ancestry_pred": str, "ancestry_pred_other": str},
        usecols=["research_id", "ancestry_pred", "ancestry_pred_other"],
    )
    require_columns(frame, {"research_id", "ancestry_pred", "ancestry_pred_other"}, path)
    frame["research_id"] = frame["research_id"].astype(str)
    frame["ancestry_pred"] = frame["ancestry_pred"].str.lower()
    frame["ancestry_pred_other"] = frame["ancestry_pred_other"].str.lower()
    return frame


def classify_european(
    frame: pd.DataFrame,
    eur_min: float,
    afr_max: float,
    amr_max: float,
    eas_max: float,
    oce_max: float,
) -> pd.Series:
    return (
        (frame["European"] >= eur_min)
        & (frame["African"] <= afr_max)
        & (frame["American"] <= amr_max)
        & (frame["East_Asian"] <= eas_max)
        & (frame["Oceanian"] <= oce_max)
    )


def write_keep(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame[["FID", "IID"]].to_csv(path, sep="\t", header=False, index=False)


def safe_fraction(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return math.nan
    return float(numerator / denominator)


def ordered_counts(series: pd.Series, order: list[str]) -> pd.Series:
    counts = series.value_counts(dropna=False)
    ordered = counts.reindex(order, fill_value=0)
    extras = counts.drop(labels=[key for key in order if key in counts.index], errors="ignore")
    return pd.concat([ordered, extras])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ours-admixture", type=Path, required=True)
    parser.add_argument("--aou-ancestry-pred", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--europeans-dir", type=Path, required=True)
    parser.add_argument("--eur-min", type=float, default=0.8)
    parser.add_argument("--afr-max", type=float, default=0.1)
    parser.add_argument("--amr-max", type=float, default=0.1)
    parser.add_argument("--eas-max", type=float, default=0.1)
    parser.add_argument("--oce-max", type=float, default=0.1)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.europeans_dir.mkdir(parents=True, exist_ok=True)

    ours = read_ours(args.ours_admixture)
    aou_pred = read_aou_pred(args.aou_ancestry_pred)

    analysis = ours.copy()
    if aou_pred is not None:
        analysis = analysis.merge(aou_pred, left_on="IID", right_on="research_id", how="left", validate="many_to_one")
        analysis["aou_pred_missing"] = analysis["ancestry_pred"].isna()
        analysis["ancestry_pred"] = analysis["ancestry_pred"].fillna("missing")
        analysis["ancestry_pred_other"] = analysis["ancestry_pred_other"].fillna("missing")
    else:
        analysis["aou_pred_missing"] = True
        analysis["ancestry_pred"] = "missing"
        analysis["ancestry_pred_other"] = "missing"

    analysis["aou_is_european"] = analysis["ancestry_pred"] == "eur"
    analysis["ours_is_european"] = classify_european(
        analysis, args.eur_min, args.afr_max, args.amr_max, args.eas_max, args.oce_max
    )
    analysis["european_set_group"] = "neither"
    analysis.loc[analysis["aou_is_european"] & analysis["ours_is_european"], "european_set_group"] = "both"
    analysis.loc[analysis["aou_is_european"] & ~analysis["ours_is_european"], "european_set_group"] = "aou_eur_only"
    analysis.loc[~analysis["aou_is_european"] & analysis["ours_is_european"], "european_set_group"] = "ours_eur_only"
    analysis["european_set_group"] = pd.Categorical(analysis["european_set_group"], OVERLAP_ORDER, ordered=True)

    ours_eur = analysis[analysis["ours_is_european"]].copy()
    aou_eur_only = analysis[analysis["aou_is_european"] & ~analysis["ours_is_european"]].copy()
    ours_eur_only = analysis[~analysis["aou_is_european"] & analysis["ours_is_european"]].copy()

    write_keep(args.europeans_dir / "classified_european_iids.txt", ours_eur)
    write_keep(args.output_dir / "ours_classified_european_iids.txt", ours_eur)
    write_keep(args.output_dir / "aou_eur_not_ours_iids.txt", aou_eur_only)
    write_keep(args.output_dir / "ours_eur_not_aou_iids.txt", ours_eur_only)

    n = len(analysis)
    aou_pred_samples = int((~analysis["aou_pred_missing"]).sum())
    aou_eur_count = int(analysis["aou_is_european"].sum())
    ours_eur_count = int(analysis["ours_is_european"].sum())
    both_count = int((analysis["aou_is_european"] & analysis["ours_is_european"]).sum())
    aou_only_count = len(aou_eur_only)
    ours_only_count = len(ours_eur_only)
    neither_count = int((~analysis["aou_is_european"] & ~analysis["ours_is_european"]).sum())
    union = both_count + aou_only_count + ours_only_count

    summary_rows = [
        ("comparison_mode", "classification_only_no_aou_rye"),
        ("ours_admixture_samples", n),
        ("aou_hard_pred_samples", aou_pred_samples),
        ("samples_in_classifier", n),
        ("samples_missing_aou_hard_pred", int(analysis["aou_pred_missing"].sum())),
        ("aou_eur_pred_count", aou_eur_count),
        ("ours_european_count", ours_eur_count),
        ("european_both", both_count),
        ("aou_eur_only", aou_only_count),
        ("ours_eur_only", ours_only_count),
        ("neither_european", neither_count),
        ("european_union", union),
        ("european_jaccard_both_over_union", safe_fraction(both_count, union)),
        ("ours_precision_if_aou_eur_is_reference", safe_fraction(both_count, ours_eur_count)),
        ("ours_recall_if_aou_eur_is_reference", safe_fraction(both_count, aou_eur_count)),
        ("ours_eur_min", args.eur_min),
        ("ours_afr_max", args.afr_max),
        ("ours_amr_max", args.amr_max),
        ("ours_eas_max", args.eas_max),
        ("ours_oce_max", args.oce_max),
    ]
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(
        args.output_dir / "aou_vs_ours.summary.tsv", sep="\t", index=False
    )

    overlap = ordered_counts(analysis["european_set_group"].astype(str), OVERLAP_ORDER).rename_axis(
        "european_set_group"
    )
    overlap = overlap.reset_index(name="count")
    overlap["fraction_of_samples"] = overlap["count"] / n if n else np.nan
    overlap.to_csv(args.output_dir / "european_set_overlap_summary.tsv", sep="\t", index=False)

    pred_counts = ordered_counts(analysis["ancestry_pred"], PRED_ORDER).rename_axis("ancestry_pred").reset_index(
        name="count"
    )
    pred_counts.to_csv(args.output_dir / "aou_pred_counts.tsv", sep="\t", index=False)

    pred_group_rows: list[dict[str, object]] = []
    for pred, group in analysis.groupby("ancestry_pred", dropna=False):
        pred_group_rows.append(
            {
                "ancestry_pred": pred,
                "count": len(group),
                "ours_european_count": int(group["ours_is_european"].sum()),
                "ours_european_rate": float(group["ours_is_european"].mean()),
                "ours_European_mean": float(group["European"].mean()),
                "ours_East_Asian_mean": float(group["East_Asian"].mean()),
                "ours_American_mean": float(group["American"].mean()),
                "ours_African_mean": float(group["African"].mean()),
                "ours_South_Asian_mean": float(group["South_Asian"].mean()),
                "ours_Oceanian_mean": float(group["Oceanian"].mean()),
            }
        )
    pd.DataFrame(pred_group_rows).to_csv(
        args.output_dir / "aou_pred_vs_ours_european_summary.tsv", sep="\t", index=False
    )

    print(f"wrote {args.europeans_dir / 'classified_european_iids.txt'}")
    print(f"ours European samples: {ours_eur_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
