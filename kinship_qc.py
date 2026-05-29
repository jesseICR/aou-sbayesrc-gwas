"""Compare this pipeline's KING kinship table to AoU-provided relatedness."""

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


KINSHIP_BINS = [
    ("duplicate_or_mz", 0.3535, math.inf),
    ("first_degree", 0.1767, 0.3535),
    ("second_degree", 0.0884, 0.1767),
    ("third_degree", 0.0442, 0.0884),
]


def pair_key(iid1: str, iid2: str) -> tuple[str, str]:
    return (iid1, iid2) if iid1 <= iid2 else (iid2, iid1)


def percentile(sorted_values: list[float], p: float) -> float:
    if not sorted_values:
        return math.nan
    idx = min(len(sorted_values) - 1, max(0, int(round((len(sorted_values) - 1) * p))))
    return sorted_values[idx]


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n == 0:
        return math.nan
    mx = sum(xs) / n
    my = sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n)
    if sx == 0 or sy == 0:
        return math.nan
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
    return cov / (sx * sy)


def read_kin0(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    records: dict[tuple[str, str], dict[str, float]] = {}
    with path.open() as handle:
        header = handle.readline().strip().split()
        cols = {name.lstrip("#"): idx for idx, name in enumerate(header)}
        required = {"IID1", "IID2", "HETHET", "IBS0", "KINSHIP"}
        missing = required.difference(cols)
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            iid1 = parts[cols["IID1"]]
            iid2 = parts[cols["IID2"]]
            records[pair_key(iid1, iid2)] = {
                "kinship": float(parts[cols["KINSHIP"]]),
                "ibs0": float(parts[cols["IBS0"]]),
                "hethet": float(parts[cols["HETHET"]]),
            }
    return records


def read_aou_relatedness(path: Path) -> dict[tuple[str, str], float]:
    records: dict[tuple[str, str], float] = {}
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cols = {name: idx for idx, name in enumerate(header)}
        required = {"i.s", "j.s", "kin"}
        missing = required.difference(cols)
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            records[pair_key(parts[cols["i.s"]], parts[cols["j.s"]])] = float(parts[cols["kin"]])
    return records


def bin_label(value: float, king_filter: float) -> str:
    for label, lo, hi in KINSHIP_BINS:
        if lo <= value < hi:
            return label
    if king_filter <= value < 0.0442:
        return "below_third_degree"
    if value < king_filter:
        return "below_king_filter"
    return "out_of_range"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ours-kin0", type=Path, required=True)
    parser.add_argument("--aou-relatedness", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--king-table-filter", type=float, default=0.035)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading ours KING table: {args.ours_kin0}")
    ours = read_kin0(args.ours_kin0)
    print(f"  ours pairs: {len(ours):,}")
    print(f"Loading AoU relatedness: {args.aou_relatedness}")
    aou = read_aou_relatedness(args.aou_relatedness)
    print(f"  AoU pairs: {len(aou):,}")

    overlap_keys = sorted(set(ours) & set(aou))
    if not overlap_keys:
        raise ValueError("No overlapping relatedness pairs between ours and AoU")

    rows = []
    ours_vals = []
    aou_vals = []
    diffs = []
    abs_diffs = []
    for key in overlap_keys:
        ours_kin = ours[key]["kinship"]
        aou_kin = aou[key]
        diff = ours_kin - aou_kin
        rows.append((key[0], key[1], ours_kin, aou_kin, diff, abs(diff), ours[key]["ibs0"], ours[key]["hethet"]))
        ours_vals.append(ours_kin)
        aou_vals.append(aou_kin)
        diffs.append(diff)
        abs_diffs.append(abs(diff))

    sorted_abs = sorted(abs_diffs)
    corr = pearson(ours_vals, aou_vals)
    mean_signed = sum(diffs) / len(diffs)
    mean_abs = sum(abs_diffs) / len(abs_diffs)

    ours_only = set(ours) - set(aou)
    aou_only = set(aou) - set(ours)
    summary_rows: list[tuple[str, object]] = [
        ("ours_total_pairs", len(ours)),
        ("aou_total_pairs", len(aou)),
        ("overlapping_pairs", len(overlap_keys)),
        ("ours_only_pairs", len(ours_only)),
        ("aou_only_pairs", len(aou_only)),
        ("king_table_filter", args.king_table_filter),
        ("pearson_r", corr),
        ("mean_signed_diff_ours_minus_aou", mean_signed),
        ("mean_abs_diff", mean_abs),
        ("median_abs_diff", percentile(sorted_abs, 0.50)),
        ("p90_abs_diff", percentile(sorted_abs, 0.90)),
        ("p95_abs_diff", percentile(sorted_abs, 0.95)),
        ("p99_abs_diff", percentile(sorted_abs, 0.99)),
        ("max_abs_diff", max(sorted_abs)),
    ]

    bin_names = ["duplicate_or_mz", "first_degree", "second_degree", "third_degree", "below_third_degree"]
    bin_rows = []
    for label in bin_names:
        both = sum(1 for _i, _j, ours_kin, aou_kin, *_rest in rows if bin_label(aou_kin, args.king_table_filter) == label)
        ours_only_count = sum(1 for key in ours_only if bin_label(ours[key]["kinship"], args.king_table_filter) == label)
        aou_only_count = sum(1 for key in aou_only if bin_label(aou[key], args.king_table_filter) == label)
        bin_rows.append((label, both, ours_only_count, aou_only_count))
        summary_rows.append((f"bin_{label}_overlap_by_aou_kinship", both))
        summary_rows.append((f"bin_{label}_ours_only", ours_only_count))
        summary_rows.append((f"bin_{label}_aou_only", aou_only_count))

    with (args.output_dir / "kinship_comparison_summary.tsv").open("w") as out:
        out.write("metric\tvalue\n")
        for key, value in summary_rows:
            out.write(f"{key}\t{value}\n")

    with (args.output_dir / "kinship_comparison_pairs.tsv").open("w") as out:
        out.write("iid1\tiid2\tours_kinship\taou_kinship\tkinship_diff\tkinship_abs_diff\tours_ibs0\tours_hethet\n")
        for row in rows:
            out.write("\t".join(str(value) for value in row) + "\n")

    with (args.output_dir / "kinship_pair_counts_by_bin.tsv").open("w") as out:
        out.write("bin\toverlap_by_aou_kinship\tours_only\taou_only\n")
        for row in bin_rows:
            out.write("\t".join(str(value) for value in row) + "\n")

    with (args.output_dir / "kinship_comparison_summary.txt").open("w") as out:
        out.write("=== Kinship Comparison: ours vs AoU ===\n\n")
        out.write(f"Ours total pairs:   {len(ours):,}\n")
        out.write(f"AoU total pairs:    {len(aou):,}\n")
        out.write(f"Overlapping pairs:  {len(overlap_keys):,}\n")
        out.write(f"Ours-only pairs:    {len(ours_only):,}\n")
        out.write(f"AoU-only pairs:     {len(aou_only):,}\n\n")
        out.write("--- Kinship Difference (ours - AoU) ---\n")
        out.write(f"Mean absolute difference:   {mean_abs:.6f}\n")
        out.write(f"Median absolute difference: {percentile(sorted_abs, 0.50):.6f}\n")
        out.write(f"90th percentile abs diff:   {percentile(sorted_abs, 0.90):.6f}\n")
        out.write(f"95th percentile abs diff:   {percentile(sorted_abs, 0.95):.6f}\n")
        out.write(f"99th percentile abs diff:   {percentile(sorted_abs, 0.99):.6f}\n")
        out.write(f"Max absolute difference:    {max(sorted_abs):.6f}\n")
        out.write(f"Mean signed difference:     {mean_signed:.6f}\n")
        out.write(f"Pearson correlation:        {corr:.6f}\n\n")
        out.write("--- Pair Counts by Bin ---\n")
        out.write("Bin\tOverlap_by_AoU_kinship\tOurs_only\tAoU_only\n")
        for row in bin_rows:
            out.write("\t".join(str(value) for value in row) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(aou_vals, ours_vals, s=2, c="black", alpha=0.5, linewidths=0)
    lo = min(min(aou_vals), min(ours_vals)) - 0.01
    hi = max(max(aou_vals), max(ours_vals)) + 0.01
    axes[0].plot([lo, hi], [lo, hi], color="red", linestyle="--", linewidth=1)
    axes[0].set_xlim(lo, hi)
    axes[0].set_ylim(lo, hi)
    axes[0].set_xlabel("AoU kinship")
    axes[0].set_ylabel("Our KING kinship")
    axes[0].set_title(f"Kinship: ours vs AoU\nr={corr:.4f}, n={len(overlap_keys):,}")
    axes[0].grid(True, linewidth=0.3, alpha=0.35)

    axes[1].hist(diffs, bins=100, color="#4c78a8", alpha=0.85)
    axes[1].axvline(0, color="red", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Kinship difference (ours - AoU)")
    axes[1].set_ylabel("Pairs")
    axes[1].set_title("Difference distribution")
    axes[1].grid(True, axis="y", linewidth=0.3, alpha=0.35)
    fig.tight_layout()
    fig.savefig(args.output_dir / "kinship_comparison_plots.png", dpi=180)
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(7, 5))
    means = [(x + y) / 2 for x, y in zip(ours_vals, aou_vals)]
    sd = math.sqrt(sum((d - mean_signed) ** 2 for d in diffs) / len(diffs))
    ax.scatter(means, diffs, s=2, c="black", alpha=0.5, linewidths=0)
    ax.axhline(mean_signed, color="red", linewidth=1, label=f"mean={mean_signed:.4f}")
    ax.axhline(mean_signed + 1.96 * sd, color="gray", linestyle="--", linewidth=1, label="+/-1.96 SD")
    ax.axhline(mean_signed - 1.96 * sd, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("Mean kinship")
    ax.set_ylabel("Kinship difference (ours - AoU)")
    ax.set_title("Bland-Altman: kinship")
    ax.grid(True, linewidth=0.3, alpha=0.35)
    ax.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(args.output_dir / "kinship_bland_altman.png", dpi=180)
    plt.close(fig2)

    print(f"overlapping pairs: {len(overlap_keys):,}")
    print(f"pearson r: {corr:.6f}")
    print(f"mean abs diff: {mean_abs:.6f}")
    print(f"wrote {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
