#!/usr/bin/env python3
"""Create compact per-chromosome REGENIE association result files."""

import argparse
import gzip
import os
from pathlib import Path


SOURCE_TO_OUTPUT = [
    ("ID", "rsid"),
    ("ALLELE1", "allele1"),
    ("A1FREQ", "a1freq"),
    ("N", "n"),
    ("BETA.Y1", "beta"),
    ("SE.Y1", "se"),
    ("LOG10P.Y1", "log10p"),
]


def expand_chroms(spec):
    out = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start_s, end_s = token.split("-", 1)
            start, end = int(start_s), int(end_s)
            if start > end:
                raise ValueError(f"invalid chromosome range: {token}")
            out.extend(range(start, end + 1))
        else:
            out.append(int(token))
    seen = set()
    chroms = []
    for chrom in out:
        if chrom < 1 or chrom > 22:
            raise ValueError(f"chromosome outside 1..22: {chrom}")
        if chrom not in seen:
            chroms.append(chrom)
            seen.add(chrom)
    if not chroms:
        raise ValueError(f"no chromosomes selected by {spec!r}")
    return chroms


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def find_source(output_dir, chrom, result_prefix):
    chrom_label = f"chr{chrom}"
    chrom_dir = output_dir / "step2" / chrom_label
    for suffix in (".regenie.gz", ".regenie"):
        candidate = chrom_dir / f"{chrom_label}_{result_prefix}{suffix}"
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise FileNotFoundError(
        f"missing current-style REGENIE result for {chrom_label}: "
        f"{chrom_dir}/{chrom_label}_{result_prefix}.regenie[.gz]"
    )


def count_rows_gzip(path):
    with gzip.open(path, "rt") as handle:
        n_lines = sum(1 for _ in handle)
    return max(0, n_lines - 1)


def export_chrom(source, target):
    with open_text(source) as src:
        header = src.readline().strip().split()
        if not header:
            raise ValueError(f"empty source file: {source}")
        col_idx = {name: i for i, name in enumerate(header)}
        missing = [src_col for src_col, _ in SOURCE_TO_OUTPUT if src_col not in col_idx]
        if missing:
            raise ValueError(f"{source} missing required columns: {', '.join(missing)}")

        indices = [col_idx[src_col] for src_col, _ in SOURCE_TO_OUTPUT]
        max_idx = max(indices)

        tmp = target.with_suffix(target.suffix + f".tmp.{os.getpid()}")
        source_rows = 0
        with gzip.open(tmp, "wt") as out:
            out.write("\t".join(out_col for _, out_col in SOURCE_TO_OUTPUT) + "\n")
            for line in src:
                fields = line.strip().split()
                if not fields:
                    continue
                if len(fields) <= max_idx:
                    raise ValueError(f"short row in {source}: {line[:120]!r}")
                out.write("\t".join(fields[i] for i in indices) + "\n")
                source_rows += 1
        os.replace(tmp, target)
    output_rows = count_rows_gzip(target)
    if output_rows != source_rows:
        raise RuntimeError(f"row-count mismatch for {target}: source={source_rows}, output={output_rows}")
    return source_rows, output_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--result-prefix", required=True)
    parser.add_argument("--chroms", required=True)
    args = parser.parse_args()

    output_dir = args.output_dir
    light_dir = output_dir / "lightweight"
    light_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    total_rows = 0
    columns = ",".join(out_col for _, out_col in SOURCE_TO_OUTPUT)
    for chrom in expand_chroms(args.chroms):
        chrom_label = f"chr{chrom}"
        source = find_source(output_dir, chrom, args.result_prefix)
        target = light_dir / f"{chrom_label}.{args.result_prefix}.regenie_lite.tsv.gz"
        source_rows, output_rows = export_chrom(source, target)
        total_rows += output_rows
        summary_rows.append((chrom_label, source, target, source_rows, output_rows, columns))
        print(f"  {chrom_label}: {output_rows} rows -> {target}")

    summary_path = light_dir / "regenie_lite.summary.tsv"
    tmp_summary = summary_path.with_suffix(summary_path.suffix + f".tmp.{os.getpid()}")
    with tmp_summary.open("w") as out:
        out.write("chrom\tsource_file\toutput_file\tsource_rows\toutput_rows\tcolumns\n")
        for row in summary_rows:
            out.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n"
            )
    os.replace(tmp_summary, summary_path)
    print(f"  Summary: {summary_path}")
    print(f"  Total lightweight rows: {total_rows}")


if __name__ == "__main__":
    main()
