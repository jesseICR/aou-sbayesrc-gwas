#!/usr/bin/env python3
"""Align ADMIXTURE reference frequencies to a PLINK .bim allele order.

The ADMIXTURE projection .P file must have one row per retained SNP, in the
same order as the PLINK .bim, with frequencies for the .bim A1 allele.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


POP_COLS = ["European", "East Asian", "American", "African", "South Asian", "Oceanian"]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-tsv", default="admixture_allele_freqs.tsv", type=Path)
    parser.add_argument("--bim", default="aou_admixture_extracted.bim", type=Path)
    parser.add_argument("--p-out", default="ref_aligned.P", type=Path)
    parser.add_argument("--snps-out", default="snps_aligned.txt", type=Path)
    parser.add_argument("--log-out", default="admixture_align_log.txt", type=Path)
    parser.add_argument("--summary-out", default="admixture_align_summary.tsv", type=Path)
    return parser.parse_args()


def load_reference(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"snp_id", "a1", "a2", *POP_COLS}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(sorted(missing))}")
        return {row["snp_id"]: row for row in reader}


def read_bim(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                raise ValueError(f"{path}:{line_number} has fewer than 6 BIM columns")
            rows.append(
                {
                    "chrom": parts[0],
                    "rsid": parts[1],
                    "cm": parts[2],
                    "pos": parts[3],
                    "a1": parts[4],
                    "a2": parts[5],
                }
            )
    return rows


def comp(allele: str) -> str | None:
    return COMPLEMENT.get(allele)


def align_row(bim_a1: str, bim_a2: str, ref_a1: str, ref_a2: str, freqs: list[float]) -> tuple[list[float] | None, str]:
    allele_set = frozenset((bim_a1, bim_a2))
    if allele_set in AMBIGUOUS:
        return None, "strand_ambiguous"

    if bim_a1 == ref_a1 and bim_a2 == ref_a2:
        return freqs, "same"
    if bim_a1 == ref_a2 and bim_a2 == ref_a1:
        return [1.0 - value for value in freqs], "swapped"
    if bim_a1 == comp(ref_a1) and bim_a2 == comp(ref_a2):
        return freqs, "strand_flip"
    if bim_a1 == comp(ref_a2) and bim_a2 == comp(ref_a1):
        return [1.0 - value for value in freqs], "swapped_strand_flip"
    return None, "allele_mismatch"


def main() -> None:
    args = parse_args()
    reference = load_reference(args.ref_tsv)
    bim_rows = read_bim(args.bim)

    counts = {
        "bim_variants": len(bim_rows),
        "reference_variants": len(reference),
        "same": 0,
        "swapped": 0,
        "strand_flip": 0,
        "swapped_strand_flip": 0,
        "strand_ambiguous": 0,
        "allele_mismatch": 0,
        "missing_from_reference": 0,
    }

    p_rows: list[list[float]] = []
    kept_snps: list[str] = []

    for row in bim_rows:
        rsid = row["rsid"]
        ref = reference.get(rsid)
        if ref is None:
            counts["missing_from_reference"] += 1
            continue

        freqs = [float(ref[col]) for col in POP_COLS]
        aligned_freqs, status = align_row(row["a1"], row["a2"], ref["a1"], ref["a2"], freqs)
        counts[status] += 1
        if aligned_freqs is None:
            continue

        p_rows.append(aligned_freqs)
        kept_snps.append(rsid)

    counts["retained"] = len(kept_snps)
    counts["excluded"] = len(bim_rows) - len(kept_snps)

    with args.snps_out.open("w") as handle:
        for rsid in kept_snps:
            handle.write(f"{rsid}\n")

    with args.p_out.open("w") as handle:
        for freqs in p_rows:
            handle.write(" ".join(f"{value:.10g}" for value in freqs) + "\n")

    summary_lines = [
        f"Allele alignment retained {counts['retained']} of {counts['bim_variants']} extracted SNPs",
        f"same={counts['same']}",
        f"swapped={counts['swapped']}",
        f"strand_flip={counts['strand_flip']}",
        f"swapped_strand_flip={counts['swapped_strand_flip']}",
        f"strand_ambiguous={counts['strand_ambiguous']}",
        f"allele_mismatch={counts['allele_mismatch']}",
        f"missing_from_reference={counts['missing_from_reference']}",
        f"populations={','.join(POP_COLS)}",
    ]
    log_text = "\n".join(summary_lines) + "\n"
    print(log_text, end="")
    args.log_out.write_text(log_text)

    with args.summary_out.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key, value in counts.items():
            handle.write(f"{key}\t{value}\n")
        handle.write("populations\t" + ",".join(POP_COLS) + "\n")


if __name__ == "__main__":
    main()
