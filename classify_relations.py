"""Classify close relationships from a plink2 KING .kin0 table.

This AoU pipeline uses fixed kinship/IBS0 thresholds and deliberately does not
apply an age-gap sibling filter because no portable AoU birth-year/month
dependency is part of this genotype pipeline.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path


def classify_relationship(
    kinship: float,
    ibs0: float,
    close_lower: float,
    first_degree_upper: float,
    ibs0_cutoff: float,
) -> str | None:
    if kinship >= first_degree_upper and ibs0 < ibs0_cutoff:
        return "identical"
    if close_lower <= kinship < first_degree_upper:
        return "sibling" if ibs0 >= ibs0_cutoff else "parent_child"
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kin0", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--close-lower", type=float, default=0.1767)
    parser.add_argument("--first-degree-upper", type=float, default=0.3535)
    parser.add_argument("--ibs0-cutoff", type=float, default=0.0012)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = args.output_dir / "close_relations.csv"
    out_summary = args.output_dir / "close_relations.summary.tsv"

    counts: Counter[str] = Counter()
    total_kin0_pairs = 0
    total_close = 0

    with args.kin0.open() as fin, out_csv.open("w", newline="") as fout:
        header = fin.readline().strip().split()
        cols = {name.lstrip("#"): idx for idx, name in enumerate(header)}
        required = {"IID1", "IID2", "HETHET", "IBS0", "KINSHIP"}
        missing = required.difference(cols)
        if missing:
            raise ValueError(f"{args.kin0}: missing required columns {sorted(missing)}")

        writer = csv.writer(fout)
        writer.writerow(["eid1", "eid2", "kinship", "ibs0", "hethet", "relationship"])
        for line in fin:
            if not line.strip():
                continue
            total_kin0_pairs += 1
            parts = line.strip().split()
            iid1 = parts[cols["IID1"]]
            iid2 = parts[cols["IID2"]]
            hethet = float(parts[cols["HETHET"]])
            ibs0 = float(parts[cols["IBS0"]])
            kinship = float(parts[cols["KINSHIP"]])
            relationship = classify_relationship(
                kinship,
                ibs0,
                args.close_lower,
                args.first_degree_upper,
                args.ibs0_cutoff,
            )
            if relationship is None:
                continue
            total_close += 1
            counts[relationship] += 1
            writer.writerow([iid1, iid2, f"{kinship:.8g}", f"{ibs0:.8g}", f"{hethet:.8g}", relationship])

    with out_summary.open("w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_king_pairs\t{total_kin0_pairs}\n")
        out.write(f"total_close_relationships\t{total_close}\n")
        out.write(f"sibling\t{counts['sibling']}\n")
        out.write(f"parent_child\t{counts['parent_child']}\n")
        out.write(f"identical\t{counts['identical']}\n")
        out.write(f"kinship_close_lower\t{args.close_lower}\n")
        out.write(f"kinship_first_degree_upper\t{args.first_degree_upper}\n")
        out.write(f"kinship_ibs0_cutoff\t{args.ibs0_cutoff}\n")
        out.write("age_gap_filter\tnot_applied\n")

    print(f"total KING pairs: {total_kin0_pairs}")
    print(f"close relationships: {total_close}")
    print(f"sibling: {counts['sibling']}")
    print(f"parent_child: {counts['parent_child']}")
    print(f"identical: {counts['identical']}")
    print(f"wrote {out_csv}")
    print(f"wrote {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
