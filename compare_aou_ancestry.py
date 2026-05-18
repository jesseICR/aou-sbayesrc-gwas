"""Compare AoU-provided ancestry estimates to this pipeline's ADMIXTURE K=6.

The outputs are aggregate QC/comparison tables plus plots. Individual-level
derived files stay in the AoU workspace bucket.
"""

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


AOU_COMPONENTS = [
    ("aou_eur", "European"),
    ("aou_eas", "East_Asian"),
    ("aou_amr", "American"),
    ("aou_afr", "African"),
    ("aou_sas", "South_Asian"),
    ("aou_mid", "Middle_Eastern"),
]
OURS_COMPONENTS = [
    ("ours_European", "European"),
    ("ours_East_Asian", "East_Asian"),
    ("ours_American", "American"),
    ("ours_African", "African"),
    ("ours_South_Asian", "South_Asian"),
    ("ours_Oceanian", "Oceanian"),
]
COMMON_COMPONENTS = [
    ("aou_eur", "ours_European", "European"),
    ("aou_eas", "ours_East_Asian", "East_Asian"),
    ("aou_amr", "ours_American", "American"),
    ("aou_afr", "ours_African", "African"),
    ("aou_sas", "ours_South_Asian", "South_Asian"),
]
AOU_PRED_ORDER = ["eur", "afr", "amr", "eas", "sas", "mid", "missing"]
OVERLAP_ORDER = ["both", "aou_eur_only", "ours_eur_only", "neither"]


def require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")


def parse_thresholds(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise ValueError("--mid-thresholds cannot be empty")
    bad = [value for value in values if value < 0 or value > 1]
    if bad:
        raise ValueError(f"--mid-thresholds must be within [0, 1], got {bad}")
    return sorted(set(values), reverse=True)


def read_ours(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"FID": str, "IID": str})
    require_columns(
        frame,
        {"FID", "IID", "European", "East_Asian", "American", "African", "South_Asian", "Oceanian"},
        path,
    )
    frame = frame.rename(
        columns={
            "European": "ours_European",
            "East_Asian": "ours_East_Asian",
            "American": "ours_American",
            "African": "ours_African",
            "South_Asian": "ours_South_Asian",
            "Oceanian": "ours_Oceanian",
        }
    )
    for column, _label in OURS_COMPONENTS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["IID"] = frame["IID"].astype(str)
    frame["ours_order"] = np.arange(len(frame))
    return frame


def read_aou_q(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"research_id": str})
    require_columns(frame, {"eur", "eas", "amr", "afr", "sas", "mid", "research_id"}, path)
    frame = frame.rename(
        columns={
            "eur": "aou_eur",
            "eas": "aou_eas",
            "amr": "aou_amr",
            "afr": "aou_afr",
            "sas": "aou_sas",
            "mid": "aou_mid",
        }
    )
    for column, _label in AOU_COMPONENTS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["research_id"] = frame["research_id"].astype(str)
    return frame


def read_aou_pred(path: Path) -> pd.DataFrame:
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


def quantile(series: pd.Series, q: float) -> float:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    if clean.empty:
        return math.nan
    return float(clean.quantile(q))


def distribution_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for source, components in [("AoU_RYE", AOU_COMPONENTS), ("ours", OURS_COMPONENTS)]:
        for column, label in components:
            values = pd.to_numeric(frame[column], errors="coerce")
            rows.append(
                {
                    "source": source,
                    "component": label,
                    "n_nonmissing": int(values.notna().sum()),
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=1)),
                    "min": float(values.min()),
                    "p001": quantile(values, 0.001),
                    "p01": quantile(values, 0.01),
                    "p05": quantile(values, 0.05),
                    "p10": quantile(values, 0.10),
                    "p25": quantile(values, 0.25),
                    "p50": quantile(values, 0.50),
                    "p75": quantile(values, 0.75),
                    "p90": quantile(values, 0.90),
                    "p95": quantile(values, 0.95),
                    "p99": quantile(values, 0.99),
                    "p999": quantile(values, 0.999),
                    "max": float(values.max()),
                }
            )
    return pd.DataFrame(rows)


def component_pair_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for aou_col, ours_col, label in COMMON_COMPONENTS:
        data = frame[[aou_col, ours_col]].dropna()
        diff = data[ours_col] - data[aou_col]
        abs_diff = diff.abs()
        rows.append(
            {
                "component": label,
                "n": len(data),
                "pearson_r": float(data[aou_col].corr(data[ours_col], method="pearson")),
                "spearman_r": float(data[aou_col].corr(data[ours_col], method="spearman")),
                "mean_aou": float(data[aou_col].mean()),
                "mean_ours": float(data[ours_col].mean()),
                "mean_diff_ours_minus_aou": float(diff.mean()),
                "median_diff_ours_minus_aou": float(diff.median()),
                "mae": float(abs_diff.mean()),
                "rmse": float(np.sqrt(np.mean(np.square(diff)))),
                "median_abs_diff": float(abs_diff.median()),
                "count_abs_diff_gt_0.02": int((abs_diff > 0.02).sum()),
                "count_abs_diff_gt_0.05": int((abs_diff > 0.05).sum()),
                "count_abs_diff_gt_0.10": int((abs_diff > 0.10).sum()),
                "pct_abs_diff_gt_0.02": float((abs_diff > 0.02).mean()),
                "pct_abs_diff_gt_0.05": float((abs_diff > 0.05).mean()),
                "pct_abs_diff_gt_0.10": float((abs_diff > 0.10).mean()),
            }
        )
    return pd.DataFrame(rows)


def dominant_component(frame: pd.DataFrame, components: list[tuple[str, str]]) -> pd.Series:
    columns = [column for column, _label in components]
    labels = {column: label for column, label in components}
    return frame[columns].idxmax(axis=1).map(labels)


def classify_european(
    frame: pd.DataFrame,
    eur_min: float,
    afr_max: float,
    amr_max: float,
    eas_max: float,
    oce_max: float,
) -> pd.Series:
    return (
        (frame["ours_European"] >= eur_min)
        & (frame["ours_African"] <= afr_max)
        & (frame["ours_American"] <= amr_max)
        & (frame["ours_East_Asian"] <= eas_max)
        & (frame["ours_Oceanian"] <= oce_max)
    )


def component_long_summary(frame: pd.DataFrame, group_col: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group_name, group in frame.groupby(group_col, dropna=False, observed=False):
        for source, components in [("AoU_RYE", AOU_COMPONENTS), ("ours", OURS_COMPONENTS)]:
            for column, label in components:
                values = pd.to_numeric(group[column], errors="coerce")
                rows.append(
                    {
                        "group": group_name,
                        "source": source,
                        "component": label,
                        "n": int(values.notna().sum()),
                        "mean": float(values.mean()) if values.notna().any() else math.nan,
                        "median": quantile(values, 0.50),
                        "p05": quantile(values, 0.05),
                        "p95": quantile(values, 0.95),
                    }
                )
    return pd.DataFrame(rows)


def mid_threshold_summary(frame: pd.DataFrame, thresholds: list[float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []
    total = len(frame)
    for threshold in thresholds:
        subset = frame[frame["aou_mid"] >= threshold]
        wide: dict[str, object] = {
            "aou_mid_threshold": threshold,
            "count": len(subset),
            "fraction_of_samples": len(subset) / total if total else math.nan,
        }
        for source, components in [("aou", AOU_COMPONENTS), ("ours", OURS_COMPONENTS)]:
            for column, label in components:
                values = pd.to_numeric(subset[column], errors="coerce")
                wide[f"{source}_{label}_mean"] = float(values.mean()) if values.notna().any() else math.nan
                wide[f"{source}_{label}_median"] = quantile(values, 0.50)
                long_rows.append(
                    {
                        "aou_mid_threshold": threshold,
                        "source": "AoU_RYE" if source == "aou" else "ours",
                        "component": label,
                        "n": int(values.notna().sum()),
                        "mean": float(values.mean()) if values.notna().any() else math.nan,
                        "median": quantile(values, 0.50),
                        "p05": quantile(values, 0.05),
                        "p95": quantile(values, 0.95),
                    }
                )
        wide_rows.append(wide)
    return pd.DataFrame(wide_rows), pd.DataFrame(long_rows)


def mid_bin_summary(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    bins = [-np.inf, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, np.inf]
    labels = ["<0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50", "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90", ">=0.90"]
    out = frame.copy()
    out["aou_mid_bin"] = pd.cut(out["aou_mid"], bins=bins, labels=labels, right=False)
    counts = out["aou_mid_bin"].value_counts(sort=False).rename_axis("aou_mid_bin").reset_index(name="count")
    total = len(out)
    counts["fraction_of_samples"] = counts["count"] / total if total else math.nan
    long = component_long_summary(out, "aou_mid_bin").rename(columns={"group": "aou_mid_bin"})
    return counts, long


def write_keep(path: Path, frame: pd.DataFrame) -> None:
    keep = frame.sort_values("ours_order")[["FID", "IID"]].copy()
    keep.to_csv(path, sep="\t", header=False, index=False)


def plot_component_scatter(frame: pd.DataFrame, aou_col: str, ours_col: str, label: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(frame[aou_col], frame[ours_col], s=0.25, c="black", alpha=1.0, linewidths=0)
    ax.plot([-0.01, 1.01], [-0.01, 1.01], color="#d62728", linewidth=1)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlabel(f"AoU RYE {label}")
    ax.set_ylabel(f"Our ADMIXTURE {label}")
    ax.set_title(f"{label}: AoU vs ours")
    ax.grid(True, linewidth=0.3, alpha=0.35)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_mid_vs_ours(frame: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), sharex=True, sharey=True)
    for ax, (ours_col, label) in zip(axes.ravel(), OURS_COMPONENTS):
        ax.scatter(frame["aou_mid"], frame[ours_col], s=0.25, c="black", alpha=1.0, linewidths=0)
        ax.set_title(label)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-0.01, 1.01)
        ax.grid(True, linewidth=0.3, alpha=0.35)
    for ax in axes[-1, :]:
        ax.set_xlabel("AoU RYE Middle_Eastern")
    for ax in axes[:, 0]:
        ax.set_ylabel("Our ADMIXTURE fraction")
    fig.suptitle("AoU Middle Eastern fraction vs our K=6 components", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fraction_distributions(frame: pd.DataFrame, components: list[tuple[str, str]], title: str, path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), sharex=True)
    bins = np.linspace(0, 1, 51)
    for ax, (column, label) in zip(axes.ravel(), components):
        ax.hist(frame[column].dropna(), bins=bins, histtype="stepfilled", color="#4c78a8", alpha=0.8)
        ax.set_title(label)
        ax.set_xlim(0, 1)
        ax.grid(True, axis="y", linewidth=0.3, alpha=0.35)
    for ax in axes[-1, :]:
        ax.set_xlabel("Ancestry fraction")
    for ax in axes[:, 0]:
        ax.set_ylabel("Samples")
    fig.suptitle(title, y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_heatmap(matrix: pd.DataFrame, title: str, path: Path, fmt: str = ".2f", cmap: str = "viridis") -> None:
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(matrix.columns)), max(4, 0.45 * len(matrix.index))))
    image = ax.imshow(matrix.values.astype(float), cmap=cmap, aspect="auto")
    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
    ax.set_title(title)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iat[i, j]
            if pd.notna(value):
                ax.text(j, i, format(float(value), fmt), ha="center", va="center", fontsize=8, color="white")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_overlap_counts(overlap: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(overlap["european_set_group"], overlap["count"], color=["#4c78a8", "#f58518", "#54a24b", "#b279a2"])
    ax.set_ylabel("Samples")
    ax.set_title("AoU EUR hard call vs our European threshold")
    ax.tick_params(axis="x", rotation=20)
    for idx, row in overlap.iterrows():
        ax.text(idx, row["count"], f"{int(row['count']):,}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_mid_threshold_means(summary: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = summary["aou_mid_threshold"]
    for _column, label in OURS_COMPONENTS:
        mean_col = f"ours_{label}_mean"
        ax.plot(x, summary[mean_col], marker="o", linewidth=1.2, label=label)
    ax.invert_xaxis()
    ax.set_xlabel("AoU Middle_Eastern fraction threshold")
    ax.set_ylabel("Mean of our ADMIXTURE fraction")
    ax.set_title("Our ancestry composition among samples with high AoU MID")
    ax.set_ylim(-0.01, 1.01)
    ax.grid(True, linewidth=0.3, alpha=0.35)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_discordant_mean_composition(frame: pd.DataFrame, path: Path) -> None:
    groups = [
        ("aou_eur_only", "AoU EUR, ours not EUR"),
        ("ours_eur_only", "Ours EUR, AoU not EUR"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharey=True)
    for col, (group_key, group_label) in enumerate(groups):
        group = frame[frame["european_set_group"].astype(str) == group_key]
        for row, (source_label, components, color) in enumerate(
            [
                ("AoU RYE fractions", AOU_COMPONENTS, "#4c78a8"),
                ("Our ADMIXTURE fractions", OURS_COMPONENTS, "#54a24b"),
            ]
        ):
            ax = axes[row, col]
            labels = [label.replace("_", " ") for _column, label in components]
            means = [float(group[column].mean()) if len(group) else math.nan for column, _label in components]
            ax.bar(labels, means, color=color, alpha=0.85)
            ax.set_ylim(0, 1)
            ax.set_title(f"{group_label}\n{source_label}; n={len(group):,}")
            ax.tick_params(axis="x", rotation=45)
            ax.grid(True, axis="y", linewidth=0.3, alpha=0.35)
            for idx, value in enumerate(means):
                if pd.notna(value):
                    ax.text(idx, value + 0.015, f"{value:.2f}", ha="center", va="bottom", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("Mean ancestry fraction")
    fig.suptitle("Mean ancestry composition of European-call discordant samples", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_discordant_component_boxplots(frame: pd.DataFrame, path: Path) -> None:
    groups = [
        ("aou_eur_only", "AoU EUR, ours not EUR"),
        ("ours_eur_only", "Ours EUR, AoU not EUR"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=True)
    for col, (group_key, group_label) in enumerate(groups):
        group = frame[frame["european_set_group"].astype(str) == group_key]
        for row, (source_label, components) in enumerate(
            [("AoU RYE", AOU_COMPONENTS), ("ours", OURS_COMPONENTS)]
        ):
            ax = axes[row, col]
            labels = [label.replace("_", " ") for _column, label in components]
            data = [group[column].dropna().to_numpy() for column, _label in components]
            boxes = ax.boxplot(data, tick_labels=labels, showfliers=False, whis=(5, 95), patch_artist=True)
            for patch in boxes["boxes"]:
                patch.set_facecolor("#d9e8f5")
            ax.set_ylim(-0.01, 1.01)
            ax.set_title(f"{group_label}: {source_label}; n={len(group):,}")
            ax.tick_params(axis="x", rotation=45)
            ax.grid(True, axis="y", linewidth=0.3, alpha=0.35)
    for ax in axes[:, 0]:
        ax.set_ylabel("Ancestry fraction")
    fig.suptitle("Component distributions in European-call discordant samples", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_discordant_european_fraction_histograms(frame: pd.DataFrame, eur_min: float, path: Path) -> None:
    groups = [
        ("aou_eur_only", "AoU EUR, ours not EUR"),
        ("ours_eur_only", "Ours EUR, AoU not EUR"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
    bins = np.linspace(0, 1, 51)
    for ax, (group_key, group_label) in zip(axes, groups):
        group = frame[frame["european_set_group"].astype(str) == group_key]
        ax.hist(group["aou_eur"].dropna(), bins=bins, histtype="step", linewidth=1.5, label="AoU RYE European")
        ax.hist(group["ours_European"].dropna(), bins=bins, histtype="step", linewidth=1.5, label="Our European")
        ax.axvline(eur_min, color="#d62728", linestyle="--", linewidth=1, label=f"ours EUR >= {eur_min:g}")
        ax.set_title(f"{group_label}; n={len(group):,}")
        ax.set_xlabel("European ancestry fraction")
        ax.set_ylabel("Samples")
        ax.set_xlim(0, 1)
        ax.grid(True, axis="y", linewidth=0.3, alpha=0.35)
        ax.legend(fontsize=8)
    fig.suptitle("European-fraction distributions in discordant sets", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ours-admixture", type=Path, required=True)
    parser.add_argument("--aou-admixture-q", type=Path, required=True)
    parser.add_argument("--aou-ancestry-pred", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--europeans-dir", type=Path, required=True)
    parser.add_argument("--eur-min", type=float, default=0.8)
    parser.add_argument("--afr-max", type=float, default=0.1)
    parser.add_argument("--amr-max", type=float, default=0.1)
    parser.add_argument("--eas-max", type=float, default=0.1)
    parser.add_argument("--oce-max", type=float, default=0.1)
    parser.add_argument("--mid-thresholds", default="0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.europeans_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    mid_thresholds = parse_thresholds(args.mid_thresholds)
    ours = read_ours(args.ours_admixture)
    aou_q = read_aou_q(args.aou_admixture_q)
    aou_pred = read_aou_pred(args.aou_ancestry_pred)

    merged = ours.merge(aou_q, left_on="IID", right_on="research_id", how="outer", indicator="ours_aou_fraction_merge")
    merged["research_id"] = merged["research_id"].fillna(merged["IID"])
    merged = merged.merge(aou_pred, on="research_id", how="left", validate="many_to_one")
    merged["aou_pred_missing"] = merged["ancestry_pred"].isna()
    merged["ancestry_pred"] = merged["ancestry_pred"].fillna("missing")
    merged["ancestry_pred_other"] = merged["ancestry_pred_other"].fillna("missing")

    analysis = merged[merged["ours_aou_fraction_merge"] == "both"].copy()
    if analysis.empty:
        raise ValueError("No overlapping samples between our ADMIXTURE TSV and AoU RYE Q file")

    analysis["aou_is_european"] = analysis["ancestry_pred"] == "eur"
    analysis["ours_is_european"] = classify_european(
        analysis, args.eur_min, args.afr_max, args.amr_max, args.eas_max, args.oce_max
    )
    analysis["ours_dominant_component"] = dominant_component(analysis, OURS_COMPONENTS)
    analysis["aou_dominant_component"] = dominant_component(analysis, AOU_COMPONENTS)
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
    aou_eur_count = int(analysis["aou_is_european"].sum())
    ours_eur_count = int(analysis["ours_is_european"].sum())
    both_count = int(((analysis["aou_is_european"]) & (analysis["ours_is_european"])).sum())
    aou_only_count = len(aou_eur_only)
    ours_only_count = len(ours_eur_only)
    neither_count = int((~analysis["aou_is_european"] & ~analysis["ours_is_european"]).sum())
    union = both_count + aou_only_count + ours_only_count

    summary_rows = [
        ("ours_admixture_samples", len(ours)),
        ("aou_fraction_samples", len(aou_q)),
        ("aou_hard_pred_samples", len(aou_pred)),
        ("samples_in_all_sources", n),
        ("samples_only_in_ours_admixture", int((merged["ours_aou_fraction_merge"] == "left_only").sum())),
        ("samples_only_in_aou_fraction_file", int((merged["ours_aou_fraction_merge"] == "right_only").sum())),
        ("samples_missing_aou_hard_pred", int(analysis["aou_pred_missing"].sum())),
        ("aou_eur_pred_count", aou_eur_count),
        ("ours_european_count", ours_eur_count),
        ("european_both", both_count),
        ("aou_eur_only", aou_only_count),
        ("ours_eur_only", ours_only_count),
        ("neither_european", neither_count),
        ("european_union", union),
        ("european_jaccard_both_over_union", both_count / union if union else math.nan),
        ("ours_precision_if_aou_eur_is_reference", both_count / ours_eur_count if ours_eur_count else math.nan),
        ("ours_recall_if_aou_eur_is_reference", both_count / aou_eur_count if aou_eur_count else math.nan),
        ("ours_eur_min", args.eur_min),
        ("ours_afr_max", args.afr_max),
        ("ours_amr_max", args.amr_max),
        ("ours_eas_max", args.eas_max),
        ("ours_oce_max", args.oce_max),
        ("aou_mid_thresholds", ",".join(f"{value:g}" for value in mid_thresholds)),
    ]
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(args.output_dir / "aou_vs_ours.summary.tsv", sep="\t", index=False)

    pair_metrics = component_pair_metrics(analysis)
    pair_metrics.to_csv(args.output_dir / "component_pair_metrics.tsv", sep="\t", index=False)

    aou_cols = [column for column, _label in AOU_COMPONENTS]
    ours_cols = [column for column, _label in OURS_COMPONENTS]
    corr = analysis[aou_cols + ours_cols].corr().loc[aou_cols, ours_cols]
    corr.index = [label for _column, label in AOU_COMPONENTS]
    corr.columns = [label for _column, label in OURS_COMPONENTS]
    corr.to_csv(args.output_dir / "component_correlation_matrix.tsv", sep="\t")

    distribution_summary(analysis).to_csv(args.output_dir / "fraction_distribution_summary.tsv", sep="\t", index=False)

    pred_counts = analysis["ancestry_pred"].value_counts(dropna=False).rename_axis("ancestry_pred").reset_index(name="count")
    pred_counts["fraction_of_samples"] = pred_counts["count"] / n
    pred_counts.to_csv(args.output_dir / "aou_pred_counts.tsv", sep="\t", index=False)

    pred_vs_dominant = pd.crosstab(analysis["ancestry_pred"], analysis["ours_dominant_component"])
    pred_vs_dominant = pred_vs_dominant.reindex(index=[x for x in AOU_PRED_ORDER if x in pred_vs_dominant.index], fill_value=0)
    pred_vs_dominant.to_csv(args.output_dir / "aou_pred_vs_ours_dominant_counts.tsv", sep="\t")
    pred_vs_dominant_pct = pred_vs_dominant.div(pred_vs_dominant.sum(axis=1).replace(0, np.nan), axis=0)
    pred_vs_dominant_pct.to_csv(args.output_dir / "aou_pred_vs_ours_dominant_rowpct.tsv", sep="\t")

    pred_group_rows: list[dict[str, object]] = []
    for pred, group in analysis.groupby("ancestry_pred", dropna=False):
        row: dict[str, object] = {
            "ancestry_pred": pred,
            "count": len(group),
            "ours_european_count": int(group["ours_is_european"].sum()),
            "ours_european_rate": float(group["ours_is_european"].mean()),
        }
        for column, label in OURS_COMPONENTS:
            row[f"ours_{label}_mean"] = float(group[column].mean())
        for column, label in AOU_COMPONENTS:
            row[f"aou_{label}_mean"] = float(group[column].mean())
        pred_group_rows.append(row)
    pd.DataFrame(pred_group_rows).to_csv(args.output_dir / "aou_pred_vs_ours_european_summary.tsv", sep="\t", index=False)

    overlap = (
        analysis["european_set_group"]
        .value_counts(sort=False)
        .rename_axis("european_set_group")
        .reset_index(name="count")
    )
    overlap["fraction_of_samples"] = overlap["count"] / n
    overlap.to_csv(args.output_dir / "european_set_overlap_summary.tsv", sep="\t", index=False)

    set_component_summary = component_long_summary(analysis, "european_set_group")
    set_component_summary.to_csv(args.output_dir / "european_set_group_ancestry_summary.tsv", sep="\t", index=False)
    set_component_summary[
        set_component_summary["group"].isin(["aou_eur_only", "ours_eur_only"])
    ].to_csv(args.output_dir / "discordant_set_component_summary.tsv", sep="\t", index=False)

    set_by_pred = pd.crosstab(analysis["european_set_group"], analysis["ancestry_pred"]).reindex(OVERLAP_ORDER, fill_value=0)
    set_by_pred.to_csv(args.output_dir / "european_set_group_counts_by_aou_pred.tsv", sep="\t")
    set_by_dom = pd.crosstab(analysis["european_set_group"], analysis["ours_dominant_component"]).reindex(OVERLAP_ORDER, fill_value=0)
    set_by_dom.to_csv(args.output_dir / "european_set_group_counts_by_ours_dominant.tsv", sep="\t")

    mid_threshold_wide, mid_threshold_long = mid_threshold_summary(analysis, mid_thresholds)
    mid_threshold_wide.to_csv(args.output_dir / "aou_mid_threshold_summary.tsv", sep="\t", index=False)
    mid_threshold_long.to_csv(args.output_dir / "aou_mid_threshold_component_summary.tsv", sep="\t", index=False)
    mid_bin_counts, mid_bin_long = mid_bin_summary(analysis)
    mid_bin_counts.to_csv(args.output_dir / "aou_mid_bin_counts.tsv", sep="\t", index=False)
    mid_bin_long.to_csv(args.output_dir / "aou_mid_bin_component_summary.tsv", sep="\t", index=False)

    for aou_col, ours_col, label in COMMON_COMPONENTS:
        plot_component_scatter(analysis, aou_col, ours_col, label, plots_dir / f"scatter_aou_vs_ours_{label}.png")
    plot_mid_vs_ours(analysis, plots_dir / "scatter_aou_mid_vs_ours_components.png")
    plot_fraction_distributions(analysis, AOU_COMPONENTS, "AoU RYE ancestry fraction distributions", plots_dir / "distribution_aou_rye_fractions.png")
    plot_fraction_distributions(analysis, OURS_COMPONENTS, "Our ADMIXTURE K=6 ancestry fraction distributions", plots_dir / "distribution_ours_admixture_fractions.png")
    plot_heatmap(corr, "Pearson correlation: AoU RYE fractions vs ours", plots_dir / "component_correlation_heatmap.png", cmap="coolwarm")
    plot_heatmap(pred_vs_dominant_pct.fillna(0), "AoU hard ancestry vs our dominant component", plots_dir / "aou_pred_vs_ours_dominant_rowpct_heatmap.png", cmap="viridis")
    plot_overlap_counts(overlap, plots_dir / "european_set_overlap_counts.png")
    plot_mid_threshold_means(mid_threshold_wide, plots_dir / "aou_mid_threshold_ours_means.png")
    plot_discordant_mean_composition(analysis, plots_dir / "discordant_mean_composition.png")
    plot_discordant_component_boxplots(analysis, plots_dir / "discordant_component_boxplots.png")
    plot_discordant_european_fraction_histograms(analysis, args.eur_min, plots_dir / "discordant_european_fraction_histograms.png")

    print(f"samples in all sources: {n}")
    print(f"AoU EUR hard calls: {aou_eur_count}")
    print(f"ours European threshold calls: {ours_eur_count}")
    print(f"both: {both_count}")
    print(f"AoU EUR only: {aou_only_count}")
    print(f"ours European only: {ours_only_count}")
    print(f"neither: {neither_count}")
    print(f"wrote {args.output_dir}")
    print(f"wrote {args.europeans_dir / 'classified_european_iids.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
