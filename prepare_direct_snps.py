"""Prepare per-chromosome direct-SNP extraction and missing-variant metadata.

Inputs:
  * a one-rsid-per-line direct SNP list
  * the SBayesRC hg38 alignment CSV (chrom,pos,ref,alt,rsid)
  * the extracted AoU SBayesRC pfiles with rsid IDs

Outputs under --output-dir:
  * chr{N}.extract.txt: direct-SNP rsids present in chr{N}.pvar
  * summary.tsv: per-chromosome requested/present/missing counts
  * missing_direct_snps.tsv: direct-SNP rsids absent from the extracted pfiles
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_direct_rsids(path: Path) -> list[str]:
    rsids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    seen: set[str] = set()
    duplicates: list[str] = []
    for rsid in rsids:
        if rsid in seen:
            duplicates.append(rsid)
        seen.add(rsid)
    if duplicates:
        examples = ", ".join(duplicates[:5])
        raise ValueError(f"{path}: duplicate rsids found, e.g. {examples}")
    return rsids


def read_alignment(path: Path, wanted: set[str]) -> dict[str, tuple[int, int, str, str]]:
    found: dict[str, tuple[int, int, str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"chrom", "pos", "ref", "alt", "rsid"}
        missing_cols = required.difference(reader.fieldnames or [])
        if missing_cols:
            raise ValueError(f"{path}: missing required columns: {sorted(missing_cols)}")
        for row in reader:
            rsid = row["rsid"]
            if rsid not in wanted:
                continue
            if rsid in found:
                raise ValueError(f"{path}: duplicate alignment row for {rsid}")
            chrom = int(row["chrom"])
            if chrom < 1 or chrom > 22:
                raise ValueError(f"{path}: direct SNP {rsid} is on non-autosome {chrom}")
            found[rsid] = (chrom, int(row["pos"]), row["ref"], row["alt"])
    return found


def scan_available_rsids(
    pfile_dir: Path,
    direct_set: set[str],
    alignment: dict[str, tuple[int, int, str, str]],
) -> dict[int, set[str]]:
    available: dict[int, set[str]] = defaultdict(set)
    seen: set[str] = set()
    for chrom in range(1, 23):
        pvar = pfile_dir / f"chr{chrom}.pvar"
        if not pvar.exists():
            raise FileNotFoundError(f"Missing extracted pvar: {pvar}")
        with pvar.open() as handle:
            header: list[str] | None = None
            for line in handle:
                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if line.startswith("#"):
                    continue
                if header is None:
                    raise ValueError(f"{pvar}: variant row encountered before #CHROM header")
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                rsid = row["ID"]
                if rsid not in direct_set:
                    continue
                if rsid in seen:
                    raise ValueError(f"Direct SNP {rsid} appears more than once in extracted pvars")
                seen.add(rsid)
                expected_chrom, expected_pos, expected_ref, expected_alt = alignment[rsid]
                expected = (f"chr{expected_chrom}", str(expected_pos), expected_ref, expected_alt)
                observed = (row["#CHROM"], row["POS"], row["REF"], row["ALT"])
                if observed != expected:
                    raise ValueError(
                        f"{rsid}: extracted pvar allele mismatch; expected {expected}, got {observed}"
                    )
                if expected_chrom != chrom:
                    raise ValueError(f"{rsid}: found in chr{chrom}.pvar, expected chr{expected_chrom}")
                available[chrom].add(rsid)
    return available


def write_outputs(
    output_dir: Path,
    direct_rsids: list[str],
    alignment: dict[str, tuple[int, int, str, str]],
    available: dict[int, set[str]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    desired_by_chrom: dict[int, list[str]] = defaultdict(list)
    for rsid in direct_rsids:
        chrom, _, _, _ = alignment[rsid]
        desired_by_chrom[chrom].append(rsid)

    summary_rows: list[tuple[int, int, int, int]] = []
    missing_rows: list[tuple[int, int, str, str, str]] = []

    for chrom in range(1, 23):
        desired = desired_by_chrom.get(chrom, [])
        present = sorted(available.get(chrom, set()), key=lambda r: (alignment[r][1], r))
        missing = [rsid for rsid in desired if rsid not in available.get(chrom, set())]
        missing.sort(key=lambda r: (alignment[r][1], r))

        extract_path = output_dir / f"chr{chrom}.extract.txt"
        extract_path.write_text("".join(f"{rsid}\n" for rsid in present))

        for rsid in missing:
            row_chrom, pos, ref, alt = alignment[rsid]
            missing_rows.append((row_chrom, pos, rsid, ref, alt))

        summary_rows.append((chrom, len(desired), len(present), len(missing)))

    with (output_dir / "summary.tsv").open("w") as handle:
        handle.write("chrom\tdesired\tavailable_in_wgs_pfiles\tmissing_from_wgs_pfiles\n")
        for chrom, desired, present, missing in summary_rows:
            handle.write(f"chr{chrom}\t{desired}\t{present}\t{missing}\n")

    with (output_dir / "missing_direct_snps.tsv").open("w") as handle:
        handle.write("chrom\tpos\trsid\tref\talt\n")
        for chrom, pos, rsid, ref, alt in sorted(missing_rows):
            handle.write(f"chr{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct-snps", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--wgs-pfile-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    direct_rsids = read_direct_rsids(args.direct_snps)
    direct_set = set(direct_rsids)
    alignment = read_alignment(args.alignment, direct_set)
    missing_alignment = [rsid for rsid in direct_rsids if rsid not in alignment]
    if missing_alignment:
        examples = ", ".join(missing_alignment[:10])
        raise ValueError(
            f"{len(missing_alignment)} direct SNP rsids are absent from the alignment CSV; "
            f"examples: {examples}"
        )
    available = scan_available_rsids(args.wgs_pfile_dir, direct_set, alignment)
    write_outputs(args.output_dir, direct_rsids, alignment, available)

    total_available = sum(len(rsids) for rsids in available.values())
    print(f"direct SNPs requested: {len(direct_rsids)}")
    print(f"direct SNPs available in WGS pfiles: {total_available}")
    print(f"direct SNPs absent from WGS pfiles: {len(direct_rsids) - total_available}")
    print(f"wrote per-chromosome direct SNP metadata to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
