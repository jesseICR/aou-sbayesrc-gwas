#!/usr/bin/env python3
"""Build HapMap3 WGS HQ extract lists from existing GWAS WGS QC metrics."""

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


def fmt_float(value: float | None) -> str:
    if value is None or math.isnan(value):
        return ""
    return f"{value:.10g}"


def read_hapmap_rsids(path: Path) -> tuple[list[str], set[str]]:
    rsids: list[str] = []
    seen: set[str] = set()
    duplicates: list[str] = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            rsid = line.strip()
            if not rsid:
                raise ValueError(f"{path}: empty rsid at line {line_no}")
            if rsid in seen:
                duplicates.append(rsid)
                if len(duplicates) >= 5:
                    break
            seen.add(rsid)
            rsids.append(rsid)
    if duplicates:
        raise ValueError(f"{path}: duplicate rsids found, e.g. {duplicates[:5]}")
    return rsids, seen


def read_liftover(path: Path, wanted: set[str]) -> tuple[dict[str, tuple[str, str, float]], Counter]:
    mapping: dict[str, tuple[str, str, float]] = {}
    counts: Counter = Counter()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader.fieldnames, {"ID", "A1_hg38", "A2_hg38", "A1Freq"}, path)
        for row in reader:
            rsid = row["ID"]
            if rsid not in wanted:
                continue
            if rsid in mapping:
                counts["duplicate_target_id"] += 1
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
            counts["loaded_target_rows"] += 1
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


def write_filter_summary(
    path: Path,
    *,
    total_requested: int,
    counts: Counter,
    per_chrom_pass: dict[str, int],
    thresholds: dict[str, float],
    liftover_counts: Counter,
) -> None:
    rows: list[tuple[int, str, int, int, int, str]] = []
    remaining = total_requested
    rows.append((0, "requested_hapmap3_rsids", total_requested, 0, total_requested, "tracked HapMap3 rsid list"))

    steps = [
        ("present_in_wgs_metrics", "drop HapMap3 rsids not found in WGS QC metrics", counts["absent_from_wgs"]),
        (
            "liftover_af_available_and_alleles_match",
            "drop variants without unambiguous SBayesRC/snp.info ALT frequency",
            counts["dropped_liftover_missing_or_allele_mismatch"],
        ),
        (
            f"abs_fit_pca_alt_freq_diff_le_{thresholds['af_diff_max']:g}",
            "drop variants whose fit-PCA ALT frequency differs from SBayesRC/snp.info",
            counts["dropped_af_diff_gt_threshold"],
        ),
        (
            f"fit_pca_maf_ge_{thresholds['maf_min']:g}",
            "drop variants with low MAF in fit-PCA Europeans",
            counts["dropped_fit_pca_maf_lt_threshold"],
        ),
        (
            f"classified_eur_missingness_le_{thresholds['geno_max']:g}",
            "drop variants with high missingness in classified Europeans",
            counts["dropped_classified_eur_geno_gt_threshold"],
        ),
    ]
    for step_no, (name, description, dropped) in enumerate(steps, start=1):
        before = remaining
        remaining -= dropped
        rows.append((step_no, name, before, dropped, remaining, description))

    with path.open("w", newline="") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"requested_hapmap3_rsids\t{total_requested}\n")
        handle.write(f"hapmap3_present_in_wgs_metrics\t{counts['present_in_wgs']}\n")
        handle.write(f"hapmap3_absent_from_wgs_metrics\t{counts['absent_from_wgs']}\n")
        handle.write(f"final_hapmap3_hq_snps\t{counts['pass']}\n")
        handle.write(f"af_diff_max\t{thresholds['af_diff_max']}\n")
        handle.write(f"fit_pca_maf_min\t{thresholds['maf_min']}\n")
        handle.write(f"classified_eur_missingness_max\t{thresholds['geno_max']}\n")
        for key in sorted(liftover_counts):
            handle.write(f"liftover_{key}\t{liftover_counts[key]}\n")
        for chrom in sorted(per_chrom_pass, key=lambda x: int(x.removeprefix("chr"))):
            handle.write(f"{chrom}_passing_hapmap3_hq_snps\t{per_chrom_pass[chrom]}\n")
        handle.write("\n")
        handle.write("step\tfilter\tinput_variants\tdropped_this_step\tremaining\tdescription\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hapmap-rsids", type=Path, required=True)
    parser.add_argument("--metrics-dir", type=Path, required=True)
    parser.add_argument("--liftover", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--af-diff-max", type=float, default=0.03)
    parser.add_argument("--maf-min", type=float, default=0.007)
    parser.add_argument("--geno-max", type=float, default=0.01)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    extracts_dir = args.output_dir / "extracts"
    extracts_dir.mkdir(parents=True, exist_ok=True)

    hapmap_order, hapmap_set = read_hapmap_rsids(args.hapmap_rsids)
    print(f"Loaded {len(hapmap_order)} HapMap3 rsids from {args.hapmap_rsids}", flush=True)
    liftover, liftover_counts = read_liftover(args.liftover, hapmap_set)
    print(f"Loaded {liftover_counts['loaded_target_rows']} target liftover rows", flush=True)

    qc_path = args.output_dir / "hapmap3_bfile_hq.variant_qc.tsv.gz"
    absent_path = args.output_dir / "hapmap3_bfile_hq.absent_from_wgs.txt"
    per_chrom_path = args.output_dir / "hapmap3_bfile_hq.per_chrom_summary.tsv"
    summary_path = args.output_dir / "hapmap3_bfile_hq.filter_summary.tsv"

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
        "pass_hapmap3_hq",
    ]

    counts: Counter = Counter()
    per_chrom_rows: list[dict[str, int | str]] = []
    per_chrom_pass: dict[str, int] = {}
    found: set[str] = set()

    with gzip.open(qc_path, "wt", newline="") as qc_handle:
        qc_writer = csv.DictWriter(qc_handle, delimiter="\t", fieldnames=qc_fields)
        qc_writer.writeheader()

        for chrom_num in range(1, 23):
            chrom = f"chr{chrom_num}"
            acount_path = args.metrics_dir / "wgs" / f"{chrom}.fit_pca.acount"
            vmiss_path = args.metrics_dir / "wgs" / f"{chrom}.our_eur.vmiss"
            if not acount_path.exists() or not vmiss_path.exists():
                raise FileNotFoundError(f"missing WGS metric inputs for {chrom}: {acount_path}, {vmiss_path}")

            chrom_counts: Counter = Counter()
            vmiss = read_vmiss(vmiss_path)
            extract_path = extracts_dir / f"{chrom}.extract.txt"
            with acount_path.open(newline="") as acount_handle, extract_path.open("w") as extract_handle:
                reader = csv.DictReader(acount_handle, delimiter="\t")
                require_columns(reader.fieldnames, {"#CHROM", "ID", "REF", "ALT", "ALT_CTS", "OBS_CT"}, acount_path)
                for row in reader:
                    rsid = row["ID"]
                    if rsid not in hapmap_set:
                        continue
                    if rsid in found:
                        raise ValueError(f"duplicate HapMap3 rsid found in WGS metrics: {rsid}")
                    found.add(rsid)
                    counts["present_in_wgs"] += 1
                    chrom_counts["present_in_wgs"] += 1

                    ref = row["REF"].upper()
                    alt = row["ALT"].upper()
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
                    if lift_af is None or math.isnan(fit_af):
                        abs_diff = math.nan
                    else:
                        abs_diff = abs(fit_af - lift_af)

                    if lift_af is None or math.isnan(abs_diff):
                        filter_step = "liftover_missing_or_allele_mismatch"
                        count_key = "dropped_liftover_missing_or_allele_mismatch"
                    elif abs_diff > args.af_diff_max:
                        filter_step = "af_diff_gt_threshold"
                        count_key = "dropped_af_diff_gt_threshold"
                    elif math.isnan(fit_maf) or fit_maf < args.maf_min:
                        filter_step = "fit_pca_maf_lt_threshold"
                        count_key = "dropped_fit_pca_maf_lt_threshold"
                    elif math.isnan(missing) or missing > args.geno_max:
                        filter_step = "classified_eur_geno_gt_threshold"
                        count_key = "dropped_classified_eur_geno_gt_threshold"
                    else:
                        filter_step = "pass"
                        count_key = "pass"
                        extract_handle.write(f"{rsid}\n")

                    counts[count_key] += 1
                    chrom_counts[count_key] += 1
                    qc_writer.writerow(
                        {
                            "chrom": row["#CHROM"],
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
                            "pass_hapmap3_hq": str(filter_step == "pass"),
                        }
                    )

            per_chrom_pass[chrom] = chrom_counts["pass"]
            per_chrom_rows.append(
                {
                    "chrom": chrom,
                    "present_in_wgs": chrom_counts["present_in_wgs"],
                    "dropped_liftover_missing_or_allele_mismatch": chrom_counts[
                        "dropped_liftover_missing_or_allele_mismatch"
                    ],
                    "dropped_af_diff_gt_threshold": chrom_counts["dropped_af_diff_gt_threshold"],
                    "dropped_fit_pca_maf_lt_threshold": chrom_counts["dropped_fit_pca_maf_lt_threshold"],
                    "dropped_classified_eur_geno_gt_threshold": chrom_counts[
                        "dropped_classified_eur_geno_gt_threshold"
                    ],
                    "passing_variants": chrom_counts["pass"],
                    "extract_file": str(extract_path),
                }
            )
            print(f"{chrom}: {chrom_counts['pass']} passing HapMap3 HQ variants", flush=True)

    counts["absent_from_wgs"] = len(hapmap_order) - len(found)
    with absent_path.open("w") as absent_handle:
        for rsid in hapmap_order:
            if rsid not in found:
                absent_handle.write(f"{rsid}\n")

    with per_chrom_path.open("w", newline="") as handle:
        fieldnames = [
            "chrom",
            "present_in_wgs",
            "dropped_liftover_missing_or_allele_mismatch",
            "dropped_af_diff_gt_threshold",
            "dropped_fit_pca_maf_lt_threshold",
            "dropped_classified_eur_geno_gt_threshold",
            "passing_variants",
            "extract_file",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_chrom_rows)

    write_filter_summary(
        summary_path,
        total_requested=len(hapmap_order),
        counts=counts,
        per_chrom_pass=per_chrom_pass,
        thresholds={"af_diff_max": args.af_diff_max, "maf_min": args.maf_min, "geno_max": args.geno_max},
        liftover_counts=liftover_counts,
    )

    print(f"Final HapMap3 HQ variants: {counts['pass']}", flush=True)
    print(f"Wrote {summary_path}", flush=True)
    print(f"Wrote {qc_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
