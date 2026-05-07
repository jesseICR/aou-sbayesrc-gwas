"""Generate per-chromosome SBayesRC variant ID + idmap files.

Reads the combined alignment CSV (sbayesrc_hg38.csv) containing all chromosomes
and writes two files per chromosome (chr1..chr22) into data/sbayesrc_ids/:

  * chr{N}.extract.txt — one variant ID per line in the form
        chr{N}:{pos}:{REF}:{ALT}
    matching the AoU acaf_threshold/plink_bed bim ID convention. Used as input
    to `plink2 --extract`.

  * chr{N}.idmap.txt — two whitespace-separated columns
        chr{N}:{pos}:{REF}:{ALT}    {rsid}
    for `plink2 --update-name` to rewrite the variant ID column to SBayesRC
    rsids during extraction.

Idempotent: skips a chromosome if both output files already exist non-empty.
"""

import argparse
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR / "data" / "support" / "sbayesrc_hg38.csv"
OUTPUT_DIR = SCRIPT_DIR / "data" / "sbayesrc_ids"


def make_aou_id(row: pd.Series) -> str:
    return f"chr{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-file", type=Path, default=INPUT_FILE,
        help="Alignment CSV (cols: chrom,pos,ref,alt,rsid)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=OUTPUT_DIR,
        help="Directory to write per-chromosome ID + idmap files",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    full_df = pd.read_csv(args.input_file)

    for chrom in range(1, 23):
        extract_path = args.output_dir / f"chr{chrom}.extract.txt"
        idmap_path = args.output_dir / f"chr{chrom}.idmap.txt"

        if (
            extract_path.exists() and extract_path.stat().st_size > 0
            and idmap_path.exists() and idmap_path.stat().st_size > 0
        ):
            print(f"chr{chrom}: skipping — outputs already exist")
            continue

        df = full_df.loc[full_df["chrom"] == chrom].copy()
        df["aou_id"] = df.apply(make_aou_id, axis=1)

        # Sanity: every row must have an rsid; aou_id must be unique within chrom.
        assert df["rsid"].notna().all(), f"chr{chrom}: missing rsid values"
        assert df["aou_id"].is_unique, f"chr{chrom}: duplicate chrom/pos/ref/alt"

        extract_path.write_text("\n".join(df["aou_id"]) + "\n")
        idmap_lines = (df["aou_id"] + "\t" + df["rsid"].astype(str)).tolist()
        idmap_path.write_text("\n".join(idmap_lines) + "\n")

        print(f"chr{chrom}: wrote {len(df)} variant IDs to {extract_path.name} + {idmap_path.name}")


if __name__ == "__main__":
    main()
