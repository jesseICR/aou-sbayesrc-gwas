"""Build the high-quality direct-SNP extract list from EUR QC metrics.

This is Step 4's local filter builder. It combines:

  * the original UKBB direct-SNP rsid list
  * plink2 --freq output in AoU EUR samples
  * plink2 --missing variant-only output in AoU EUR samples
  * SBayesRC liftover alleles/frequencies

and writes a sequential filter summary plus a final rsid extract list.
"""

import argparse
from pathlib import Path

import pandas as pd


def require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")


def read_direct_rsids(path: Path) -> pd.DataFrame:
    rsids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(rsids) != len(set(rsids)):
        duplicated = pd.Series(rsids)[pd.Series(rsids).duplicated()].head(5).tolist()
        raise ValueError(f"{path}: duplicate rsids found, e.g. {duplicated}")
    return pd.DataFrame({"rsid": rsids, "direct_list_order": range(1, len(rsids) + 1)})


def read_afreq(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(frame, {"ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"}, path)
    frame = frame.rename(
        columns={
            "#CHROM": "chrom",
            "ID": "rsid",
            "REF": "aou_ref",
            "ALT": "aou_alt",
            "ALT_FREQS": "aou_eur_alt_freq",
            "OBS_CT": "aou_eur_obs_ct",
        }
    )
    frame["aou_eur_alt_freq"] = pd.to_numeric(frame["aou_eur_alt_freq"], errors="coerce")
    frame["aou_eur_obs_ct"] = pd.to_numeric(frame["aou_eur_obs_ct"], errors="coerce").astype("Int64")
    return frame[["rsid", "chrom", "aou_ref", "aou_alt", "aou_eur_alt_freq", "aou_eur_obs_ct"]]


def read_vmiss(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(frame, {"ID", "F_MISS"}, path)
    frame = frame.rename(columns={"ID": "rsid", "F_MISS": "aou_eur_missing_rate"})
    frame["aou_eur_missing_rate"] = pd.to_numeric(frame["aou_eur_missing_rate"], errors="coerce")
    return frame[["rsid", "aou_eur_missing_rate"]]


def read_liftover(path: Path, wanted: set[str]) -> pd.DataFrame:
    usecols = ["ID", "A1_hg38", "A2_hg38", "A1Freq", "status", "strand_flip", "ref_match"]
    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=usecols, dtype=str, chunksize=1_000_000):
        hit = chunk[chunk["ID"].isin(wanted)]
        if not hit.empty:
            parts.append(hit.copy())
    if parts:
        frame = pd.concat(parts, ignore_index=True)
    else:
        frame = pd.DataFrame(columns=usecols)
    duplicate_count = int(frame.duplicated("ID").sum())
    if duplicate_count:
        examples = frame.loc[frame.duplicated("ID"), "ID"].head(5).tolist()
        raise ValueError(f"{path}: duplicate target rsids found, e.g. {examples}")
    frame = frame.rename(
        columns={
            "ID": "rsid",
            "A1_hg38": "liftover_A1_hg38",
            "A2_hg38": "liftover_A2_hg38",
            "A1Freq": "liftover_A1Freq",
        }
    )
    frame["liftover_A1Freq"] = pd.to_numeric(frame["liftover_A1Freq"], errors="coerce")
    return frame


def classify_rows(
    frame: pd.DataFrame,
    af_diff_max: float,
    maf_min: float,
    missing_max: float,
) -> pd.DataFrame:
    for column in ["aou_ref", "aou_alt", "liftover_A1_hg38", "liftover_A2_hg38"]:
        frame[column] = frame[column].str.upper()

    a1_is_alt = (frame["liftover_A1_hg38"] == frame["aou_alt"]) & (
        frame["liftover_A2_hg38"] == frame["aou_ref"]
    )
    a1_is_ref = (frame["liftover_A1_hg38"] == frame["aou_ref"]) & (
        frame["liftover_A2_hg38"] == frame["aou_alt"]
    )

    frame["liftover_freq_alignment"] = "allele_mismatch"
    frame.loc[frame["liftover_A1_hg38"].isna(), "liftover_freq_alignment"] = "missing_from_liftover"
    frame.loc[a1_is_alt, "liftover_freq_alignment"] = "A1_is_alt"
    frame.loc[a1_is_ref, "liftover_freq_alignment"] = "A1_is_ref"

    frame["liftover_alt_freq"] = pd.NA
    frame.loc[a1_is_alt, "liftover_alt_freq"] = frame.loc[a1_is_alt, "liftover_A1Freq"]
    frame.loc[a1_is_ref, "liftover_alt_freq"] = 1.0 - frame.loc[a1_is_ref, "liftover_A1Freq"]
    frame["liftover_alt_freq"] = pd.to_numeric(frame["liftover_alt_freq"], errors="coerce")

    frame["abs_alt_freq_diff"] = (frame["aou_eur_alt_freq"] - frame["liftover_alt_freq"]).abs()
    frame["aou_eur_maf"] = frame["aou_eur_alt_freq"].where(
        frame["aou_eur_alt_freq"] <= 0.5, 1.0 - frame["aou_eur_alt_freq"]
    )

    frame["filter_step"] = "pass"
    frame.loc[frame["aou_eur_alt_freq"].isna(), "filter_step"] = "absent_from_direct_bfile"
    frame.loc[
        (frame["filter_step"] == "pass") & frame["liftover_alt_freq"].isna(),
        "filter_step",
    ] = "liftover_missing_or_allele_mismatch"
    frame.loc[
        (frame["filter_step"] == "pass") & (frame["abs_alt_freq_diff"] > af_diff_max),
        "filter_step",
    ] = "af_diff_gt_threshold"
    frame.loc[
        (frame["filter_step"] == "pass") & (frame["aou_eur_maf"] < maf_min),
        "filter_step",
    ] = "eur_maf_lt_threshold"
    frame.loc[
        (frame["filter_step"] == "pass") & (frame["aou_eur_missing_rate"] > missing_max),
        "filter_step",
    ] = "eur_missingness_gt_threshold"
    frame["pass_hq_direct"] = frame["filter_step"] == "pass"
    return frame


def write_filter_summary(
    path: Path,
    frame: pd.DataFrame,
    af_diff_max: float,
    maf_min: float,
    missing_max: float,
) -> None:
    total = len(frame)
    remaining = total
    rows: list[tuple[int, str, int, int, int, str]] = []

    direct_present = frame["aou_eur_alt_freq"].notna()
    dropped = int((~direct_present).sum())
    remaining -= dropped
    rows.append((1, "present_in_aou_direct_bfile", total, dropped, remaining, "drop direct-list rsids absent from AoU direct bfile"))

    aligned = direct_present & frame["liftover_alt_freq"].notna()
    dropped = int((direct_present & ~aligned).sum())
    before = remaining
    remaining -= dropped
    rows.append((2, "liftover_af_available_and_alleles_match", before, dropped, remaining, "drop rsids without an unambiguous liftover ALT frequency"))

    af_pass = aligned & (frame["abs_alt_freq_diff"] <= af_diff_max)
    dropped = int((aligned & ~af_pass).sum())
    before = remaining
    remaining -= dropped
    rows.append((3, f"abs_alt_freq_diff_le_{af_diff_max:g}", before, dropped, remaining, "drop variants whose AoU EUR ALT frequency differs too much from SBayesRC"))

    maf_pass = af_pass & (frame["aou_eur_maf"] >= maf_min)
    dropped = int((af_pass & ~maf_pass).sum())
    before = remaining
    remaining -= dropped
    rows.append((4, f"aou_eur_maf_ge_{maf_min:g}", before, dropped, remaining, "drop variants with low EUR minor allele frequency"))

    miss_pass = maf_pass & (frame["aou_eur_missing_rate"] <= missing_max)
    dropped = int((maf_pass & ~miss_pass).sum())
    before = remaining
    remaining -= dropped
    rows.append((5, f"aou_eur_missingness_le_{missing_max:g}", before, dropped, remaining, "drop variants with high EUR missingness"))

    with path.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"requested_direct_snps\t{total}\n")
        handle.write(f"af_diff_max\t{af_diff_max}\n")
        handle.write(f"aou_eur_maf_min\t{maf_min}\n")
        handle.write(f"aou_eur_missing_rate_max\t{missing_max}\n")
        handle.write(f"final_hq_direct_snps\t{remaining}\n")
        handle.write("\n")
        handle.write("step\tfilter\tinput_variants\tdropped_this_step\tremaining_from_original_direct_list\tdescription\n")
        handle.write(f"0\trequested_direct_snps\t{total}\t0\t{total}\toriginal UKBB direct-SNP rsid list\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct-snps", type=Path, required=True)
    parser.add_argument("--afreq", type=Path, required=True)
    parser.add_argument("--vmiss", type=Path, required=True)
    parser.add_argument("--liftover", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--af-diff-max", type=float, default=0.04)
    parser.add_argument("--maf-min", type=float, default=0.007)
    parser.add_argument("--missing-max", type=float, default=0.05)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    direct = read_direct_rsids(args.direct_snps)
    afreq = read_afreq(args.afreq)
    vmiss = read_vmiss(args.vmiss)
    liftover = read_liftover(args.liftover, set(afreq["rsid"]))

    frame = direct.merge(afreq, on="rsid", how="left", validate="one_to_one")
    frame = frame.merge(vmiss, on="rsid", how="left", validate="one_to_one")
    frame = frame.merge(liftover, on="rsid", how="left", validate="one_to_one")
    frame = classify_rows(frame, args.af_diff_max, args.maf_min, args.missing_max)

    variant_qc = args.output_dir / "chr1_22_merged_hq.variant_qc.tsv"
    extract = args.output_dir / "chr1_22_merged_hq.extract.txt"
    summary = args.output_dir / "chr1_22_merged_hq.filter_summary.tsv"
    params = args.output_dir / "chr1_22_merged_hq.params.tsv"

    output_columns = [
        "rsid",
        "direct_list_order",
        "chrom",
        "aou_ref",
        "aou_alt",
        "aou_eur_alt_freq",
        "liftover_alt_freq",
        "abs_alt_freq_diff",
        "aou_eur_maf",
        "aou_eur_missing_rate",
        "aou_eur_obs_ct",
        "liftover_A1_hg38",
        "liftover_A2_hg38",
        "liftover_A1Freq",
        "liftover_freq_alignment",
        "status",
        "strand_flip",
        "ref_match",
        "filter_step",
        "pass_hq_direct",
    ]
    frame[output_columns].to_csv(variant_qc, sep="\t", index=False)
    passing = frame.loc[frame["pass_hq_direct"], "rsid"].tolist()
    extract.write_text("".join(f"{rsid}\n" for rsid in passing))
    write_filter_summary(summary, frame, args.af_diff_max, args.maf_min, args.missing_max)

    with params.open("w") as handle:
        handle.write("parameter\tvalue\n")
        handle.write(f"af_diff_max\t{args.af_diff_max}\n")
        handle.write(f"aou_eur_maf_min\t{args.maf_min}\n")
        handle.write(f"aou_eur_missing_rate_max\t{args.missing_max}\n")
        handle.write(f"requested_direct_snps\t{len(direct)}\n")
        handle.write(f"direct_bfile_variants\t{int(frame['aou_eur_alt_freq'].notna().sum())}\n")
        handle.write(f"liftover_file_size\t{args.liftover.stat().st_size}\n")

    print(f"requested direct SNPs: {len(direct)}")
    print(f"present in AoU direct bfile: {int(frame['aou_eur_alt_freq'].notna().sum())}")
    print(f"passing high-quality direct filters: {len(passing)}")
    print(f"wrote {extract}")
    print(f"wrote {summary}")
    print(f"wrote {variant_qc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
