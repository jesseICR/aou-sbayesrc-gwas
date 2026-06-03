#!/usr/bin/env python3
"""Build final GWAS Step 1/Step 2 variant filters from PLINK QC metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path
from typing import IO


def require_columns(fieldnames: list[str] | None, columns: set[str], path: Path) -> None:
    missing = columns.difference(fieldnames or [])
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")


def read_liftover(path: Path) -> tuple[dict[str, tuple[str, str, float]], Counter]:
    mapping: dict[str, tuple[str, str, float]] = {}
    counts: Counter = Counter()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader.fieldnames, {"ID", "A1_hg38", "A2_hg38", "A1Freq"}, path)
        for row in reader:
            rsid = row["ID"]
            if not rsid or rsid in mapping:
                counts["duplicate_or_empty_id"] += 1
                continue
            try:
                freq = float(row["A1Freq"])
            except (TypeError, ValueError):
                counts["bad_a1freq"] += 1
                continue
            a1 = row["A1_hg38"].upper()
            a2 = row["A2_hg38"].upper()
            if not a1 or not a2 or math.isnan(freq):
                counts["bad_allele_or_freq"] += 1
                continue
            mapping[rsid] = (a1, a2, freq)
            counts["loaded"] += 1
    return mapping, counts


def read_vmiss(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, {"ID", "F_MISS"}, path)
        for row in reader:
            try:
                out[row["ID"]] = float(row["F_MISS"])
            except (TypeError, ValueError):
                out[row["ID"]] = math.nan
    return out


def liftover_alt_freq(
    liftover: dict[str, tuple[str, str, float]], rsid: str, ref: str, alt: str
) -> tuple[float | None, str]:
    hit = liftover.get(rsid)
    if hit is None:
        return None, "missing_from_liftover"
    a1, a2, freq = hit
    ref = ref.upper()
    alt = alt.upper()
    if a1 == alt and a2 == ref:
        return freq, "A1_is_alt"
    if a1 == ref and a2 == alt:
        return 1.0 - freq, "A1_is_ref"
    return None, "allele_mismatch"


def fmt_float(value: float | None) -> str:
    if value is None or math.isnan(value):
        return ""
    return f"{value:.10g}"


def open_text(path: Path, mode: str = "wt") -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def write_filter_steps(path: Path, label: str, counts: dict[str, int], thresholds: dict[str, float]) -> None:
    total = counts["total"]
    remaining = total
    rows: list[tuple[int, str, int, int, int, str]] = [
        (0, "source_variants", total, 0, total, f"{label} source variant set"),
    ]

    for step, metric, description in [
        (1, "liftover_alt_freq_available_and_alleles_match", "drop variants without unambiguous SBayesRC/snp.info ALT frequency"),
        (2, f"abs_fit_pca_alt_freq_diff_le_{thresholds['af_diff_max']:g}", "drop variants whose fit_pca ALT frequency differs from SBayesRC/snp.info"),
        (3, f"fit_pca_maf_ge_{thresholds['maf_min']:g}", "drop variants with low MAF in fit_pca_iids"),
        (4, f"classified_eur_geno_le_{thresholds['geno_max']:g}", "drop variants with high missingness in classified European samples"),
    ]:
        dropped = counts[f"dropped_step_{step}"]
        before = remaining
        remaining -= dropped
        rows.append((step, metric, before, dropped, remaining, description))

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["step", "filter", "input_variants", "dropped_this_step", "remaining", "description"])
        writer.writerows(rows)


def process_metric_set(
    *,
    label: str,
    acount_path: Path,
    vmiss_path: Path,
    liftover: dict[str, tuple[str, str, float]],
    af_diff_max: float,
    maf_min: float,
    geno_max: float,
    extract_path: Path,
    qc_path: Path,
    passing_af_writer: csv.writer | None,
    per_set_af_path: Path | None,
) -> dict[str, int | float | str]:
    vmiss = read_vmiss(vmiss_path)
    extract_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    if per_set_af_path is not None:
        per_set_af_path.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter = Counter()
    alignment_counts: Counter = Counter()
    min_mac = math.inf
    min_maf = math.inf
    max_missing = -math.inf
    max_abs_diff = -math.inf

    qc_fields = [
        "chrom",
        "rsid",
        "ref",
        "alt",
        "fit_pca_alt_freq",
        "fit_pca_maf",
        "fit_pca_mac",
        "fit_pca_obs_ct",
        "classified_eur_missing_rate",
        "liftover_alt_freq",
        "abs_alt_freq_diff",
        "liftover_alignment",
        "filter_step",
        "pass_gwas_filter",
    ]
    af_fields = [
        "chrom",
        "rsid",
        "ref",
        "alt",
        "fit_pca_alt_freq",
        "fit_pca_maf",
        "fit_pca_mac",
        "fit_pca_obs_ct",
        "liftover_alt_freq",
        "abs_alt_freq_diff",
    ]

    per_set_af_handle: IO[str] | None = None
    per_set_af_writer: csv.writer | None = None
    if per_set_af_path is not None:
        per_set_af_handle = open_text(per_set_af_path)
        per_set_af_writer = csv.writer(per_set_af_handle, delimiter="\t")
        per_set_af_writer.writerow(af_fields)

    with acount_path.open(newline="") as acount_handle, extract_path.open("w") as extract_handle, open_text(qc_path) as qc_handle:
        reader = csv.DictReader(acount_handle, delimiter="\t")
        require_columns(reader.fieldnames, {"#CHROM", "ID", "REF", "ALT", "ALT_CTS", "OBS_CT"}, acount_path)
        qc_writer = csv.DictWriter(qc_handle, delimiter="\t", fieldnames=qc_fields)
        qc_writer.writeheader()

        for row in reader:
            chrom = row["#CHROM"]
            rsid = row["ID"]
            ref = row["REF"].upper()
            alt = row["ALT"].upper()
            counts["total"] += 1

            try:
                alt_ct = float(row["ALT_CTS"])
                obs_ct = float(row["OBS_CT"])
            except (TypeError, ValueError):
                alt_ct = math.nan
                obs_ct = math.nan

            if obs_ct > 0 and not math.isnan(alt_ct):
                fit_af = alt_ct / obs_ct
                fit_maf = min(fit_af, 1.0 - fit_af)
                fit_mac = min(alt_ct, obs_ct - alt_ct)
            else:
                fit_af = math.nan
                fit_maf = math.nan
                fit_mac = math.nan

            missing = vmiss.get(rsid, math.nan)
            lift_af, alignment = liftover_alt_freq(liftover, rsid, ref, alt)
            alignment_counts[alignment] += 1
            if lift_af is None or math.isnan(fit_af):
                abs_diff = math.nan
            else:
                abs_diff = abs(fit_af - lift_af)

            if lift_af is None or math.isnan(abs_diff):
                filter_step = "liftover_missing_or_allele_mismatch"
                counts["dropped_step_1"] += 1
            elif abs_diff > af_diff_max:
                filter_step = "af_diff_gt_threshold"
                counts["dropped_step_2"] += 1
            elif math.isnan(fit_maf) or fit_maf < maf_min:
                filter_step = "fit_pca_maf_lt_threshold"
                counts["dropped_step_3"] += 1
            elif math.isnan(missing) or missing > geno_max:
                filter_step = "classified_eur_geno_gt_threshold"
                counts["dropped_step_4"] += 1
            else:
                filter_step = "pass"
                counts["pass"] += 1
                extract_handle.write(f"{rsid}\n")
                min_mac = min(min_mac, fit_mac)
                min_maf = min(min_maf, fit_maf)
                max_missing = max(max_missing, missing)
                max_abs_diff = max(max_abs_diff, abs_diff)
                af_row = [
                    chrom,
                    rsid,
                    ref,
                    alt,
                    fmt_float(fit_af),
                    fmt_float(fit_maf),
                    fmt_float(fit_mac),
                    fmt_float(obs_ct),
                    fmt_float(lift_af),
                    fmt_float(abs_diff),
                ]
                if passing_af_writer is not None:
                    passing_af_writer.writerow(af_row)
                if per_set_af_writer is not None:
                    per_set_af_writer.writerow(af_row)

            qc_writer.writerow(
                {
                    "chrom": chrom,
                    "rsid": rsid,
                    "ref": ref,
                    "alt": alt,
                    "fit_pca_alt_freq": fmt_float(fit_af),
                    "fit_pca_maf": fmt_float(fit_maf),
                    "fit_pca_mac": fmt_float(fit_mac),
                    "fit_pca_obs_ct": fmt_float(obs_ct),
                    "classified_eur_missing_rate": fmt_float(missing),
                    "liftover_alt_freq": fmt_float(lift_af),
                    "abs_alt_freq_diff": fmt_float(abs_diff),
                    "liftover_alignment": alignment,
                    "filter_step": filter_step,
                    "pass_gwas_filter": str(filter_step == "pass"),
                }
            )

    if per_set_af_handle is not None:
        per_set_af_handle.close()

    summary: dict[str, int | float | str] = {
        "label": label,
        "total_variants": counts["total"],
        "dropped_liftover_missing_or_allele_mismatch": counts["dropped_step_1"],
        "dropped_af_diff_gt_threshold": counts["dropped_step_2"],
        "dropped_fit_pca_maf_lt_threshold": counts["dropped_step_3"],
        "dropped_classified_eur_geno_gt_threshold": counts["dropped_step_4"],
        "passing_variants": counts["pass"],
        "min_fit_pca_mac_passing": "" if counts["pass"] == 0 else min_mac,
        "min_fit_pca_maf_passing": "" if counts["pass"] == 0 else min_maf,
        "max_classified_eur_missing_rate_passing": "" if counts["pass"] == 0 else max_missing,
        "max_abs_alt_freq_diff_passing": "" if counts["pass"] == 0 else max_abs_diff,
    }
    for status, count in sorted(alignment_counts.items()):
        summary[f"alignment_{status}"] = count
    return summary


def write_summary(path: Path, rows: list[dict[str, int | float | str]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_params(path: Path, args: argparse.Namespace, liftover_counts: Counter) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["parameter", "value"])
        writer.writerow(["step1_af_diff_max", args.step1_af_diff_max])
        writer.writerow(["step1_maf_min", args.step1_maf_min])
        writer.writerow(["step1_geno_max", args.step1_geno_max])
        writer.writerow(["step2_af_diff_max", args.step2_af_diff_max])
        writer.writerow(["step2_maf_min", args.step2_maf_min])
        writer.writerow(["step2_geno_max", args.step2_geno_max])
        writer.writerow(["liftover_file", args.liftover])
        writer.writerow(["liftover_loaded_rows", liftover_counts["loaded"]])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--liftover", type=Path, required=True)
    parser.add_argument("--metrics-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--step1-af-diff-max", type=float, required=True)
    parser.add_argument("--step1-maf-min", type=float, required=True)
    parser.add_argument("--step1-geno-max", type=float, required=True)
    parser.add_argument("--step2-af-diff-max", type=float, required=True)
    parser.add_argument("--step2-maf-min", type=float, required=True)
    parser.add_argument("--step2-geno-max", type=float, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    step1_dir = args.output_dir / "step1_direct"
    step2_dir = args.output_dir / "step2_wgs"
    step2_extract_dir = step2_dir / "extracts"
    step2_qc_dir = step2_dir / "qc"
    step2_af_dir = step2_dir / "fit_pca_af"
    for directory in [step1_dir, step2_extract_dir, step2_qc_dir, step2_af_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    print(f"Loading SBayesRC liftover frequencies: {args.liftover}")
    liftover, liftover_counts = read_liftover(args.liftover)
    print(f"Loaded {liftover_counts['loaded']} liftover rows")

    af_fields = [
        "chrom",
        "rsid",
        "ref",
        "alt",
        "fit_pca_alt_freq",
        "fit_pca_maf",
        "fit_pca_mac",
        "fit_pca_obs_ct",
        "liftover_alt_freq",
        "abs_alt_freq_diff",
    ]

    rows: list[dict[str, int | float | str]] = []
    direct_summary = process_metric_set(
        label="step1_direct",
        acount_path=args.metrics_dir / "direct" / "direct_hq.fit_pca.acount",
        vmiss_path=args.metrics_dir / "direct" / "direct_hq.our_eur.vmiss",
        liftover=liftover,
        af_diff_max=args.step1_af_diff_max,
        maf_min=args.step1_maf_min,
        geno_max=args.step1_geno_max,
        extract_path=step1_dir / "chr1_22_merged_gwas_step1.extract.txt",
        qc_path=step1_dir / "chr1_22_merged_gwas_step1.variant_qc.tsv",
        passing_af_writer=None,
        per_set_af_path=step1_dir / "chr1_22_merged_gwas_step1.fit_pca_alt_freqs.tsv",
    )
    rows.append(direct_summary)
    write_filter_steps(
        step1_dir / "gwas_step1_direct.filter_steps.tsv",
        "Step 1 direct bfile",
        {
            "total": int(direct_summary["total_variants"]),
            "dropped_step_1": int(direct_summary["dropped_liftover_missing_or_allele_mismatch"]),
            "dropped_step_2": int(direct_summary["dropped_af_diff_gt_threshold"]),
            "dropped_step_3": int(direct_summary["dropped_fit_pca_maf_lt_threshold"]),
            "dropped_step_4": int(direct_summary["dropped_classified_eur_geno_gt_threshold"]),
        },
        {"af_diff_max": args.step1_af_diff_max, "maf_min": args.step1_maf_min, "geno_max": args.step1_geno_max},
    )

    combined_af_path = step2_af_dir / "gwas_step2_fit_pca_alt_freqs_passing.tsv.gz"
    with gzip.open(combined_af_path, "wt", newline="") as combined_handle:
        combined_writer = csv.writer(combined_handle, delimiter="\t")
        combined_writer.writerow(af_fields)
        genome = Counter()
        genome_min_mac = math.inf
        genome_min_maf = math.inf
        genome_max_missing = -math.inf
        genome_max_diff = -math.inf

        for chrom in range(1, 23):
            label = f"chr{chrom}"
            print(f"Processing {label}")
            summary = process_metric_set(
                label=label,
                acount_path=args.metrics_dir / "wgs" / f"{label}.fit_pca.acount",
                vmiss_path=args.metrics_dir / "wgs" / f"{label}.our_eur.vmiss",
                liftover=liftover,
                af_diff_max=args.step2_af_diff_max,
                maf_min=args.step2_maf_min,
                geno_max=args.step2_geno_max,
                extract_path=step2_extract_dir / f"{label}.extract.txt",
                qc_path=step2_qc_dir / f"{label}.variant_qc.tsv.gz",
                passing_af_writer=combined_writer,
                per_set_af_path=step2_af_dir / f"{label}.fit_pca_alt_freqs_passing.tsv.gz",
            )
            rows.append(summary)
            genome["total_variants"] += int(summary["total_variants"])
            genome["dropped_liftover_missing_or_allele_mismatch"] += int(summary["dropped_liftover_missing_or_allele_mismatch"])
            genome["dropped_af_diff_gt_threshold"] += int(summary["dropped_af_diff_gt_threshold"])
            genome["dropped_fit_pca_maf_lt_threshold"] += int(summary["dropped_fit_pca_maf_lt_threshold"])
            genome["dropped_classified_eur_geno_gt_threshold"] += int(summary["dropped_classified_eur_geno_gt_threshold"])
            genome["passing_variants"] += int(summary["passing_variants"])
            if int(summary["passing_variants"]) > 0:
                genome_min_mac = min(genome_min_mac, float(summary["min_fit_pca_mac_passing"]))
                genome_min_maf = min(genome_min_maf, float(summary["min_fit_pca_maf_passing"]))
                genome_max_missing = max(genome_max_missing, float(summary["max_classified_eur_missing_rate_passing"]))
                genome_max_diff = max(genome_max_diff, float(summary["max_abs_alt_freq_diff_passing"]))

    rows.append(
        {
            "label": "step2_wgs_chr1_22",
            "total_variants": genome["total_variants"],
            "dropped_liftover_missing_or_allele_mismatch": genome["dropped_liftover_missing_or_allele_mismatch"],
            "dropped_af_diff_gt_threshold": genome["dropped_af_diff_gt_threshold"],
            "dropped_fit_pca_maf_lt_threshold": genome["dropped_fit_pca_maf_lt_threshold"],
            "dropped_classified_eur_geno_gt_threshold": genome["dropped_classified_eur_geno_gt_threshold"],
            "passing_variants": genome["passing_variants"],
            "min_fit_pca_mac_passing": "" if genome["passing_variants"] == 0 else genome_min_mac,
            "min_fit_pca_maf_passing": "" if genome["passing_variants"] == 0 else genome_min_maf,
            "max_classified_eur_missing_rate_passing": "" if genome["passing_variants"] == 0 else genome_max_missing,
            "max_abs_alt_freq_diff_passing": "" if genome["passing_variants"] == 0 else genome_max_diff,
        }
    )

    write_summary(args.output_dir / "gwas_genotype_qc.summary.tsv", rows)
    write_params(args.output_dir / "gwas_genotype_qc.params.tsv", args, liftover_counts)
    write_filter_steps(
        step2_dir / "gwas_step2_wgs.filter_steps.tsv",
        "Step 2 WGS pfiles",
        {
            "total": genome["total_variants"],
            "dropped_step_1": genome["dropped_liftover_missing_or_allele_mismatch"],
            "dropped_step_2": genome["dropped_af_diff_gt_threshold"],
            "dropped_step_3": genome["dropped_fit_pca_maf_lt_threshold"],
            "dropped_step_4": genome["dropped_classified_eur_geno_gt_threshold"],
        },
        {"af_diff_max": args.step2_af_diff_max, "maf_min": args.step2_maf_min, "geno_max": args.step2_geno_max},
    )

    print(f"Wrote GWAS genotype filter outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
