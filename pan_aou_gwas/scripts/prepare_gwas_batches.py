#!/usr/bin/env python3
"""Prepare idempotent batched PLINK2 phenotype files for pan-AoU GWAS."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from pan_aou_gwas import write_batch_pheno


def read_manifest(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def dedupe_manifest_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Keep one row per PLINK phenotype name.

    Some AoU survey concepts are repeated under multiple live
    question_concept_ids while sharing the same stable item/answer identifier.
    The phenotype builder writes those to the same durable paths. PLINK2
    rejects duplicated --pheno-name values, so batch prep must collapse exact
    generated phenotype duplicates before constructing combined phenotype TSVs.
    """
    unique: list[dict[str, str]] = []
    duplicates: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        pheno_name = row.get("pheno_name") or pheno_name_from_file(Path(row["pheno_path"]))
        key = pheno_name
        if key in seen:
            dup = dict(row)
            dup["duplicate_key"] = key
            duplicates.append(dup)
            continue
        seen.add(key)
        unique.append(row)
    return unique, duplicates


def pheno_name_from_file(path: Path) -> str:
    with open(path, newline="") as f:
        header = f.readline().rstrip("\n").split("\t")
    if len(header) < 3:
        raise ValueError(f"phenotype file has fewer than 3 columns: {path}")
    return header[2]


def output_complete(row: dict[str, str]) -> bool:
    glm = Path(row.get("glm", ""))
    sumstats = Path(row.get("sumstats", ""))
    return glm.exists() and glm.stat().st_size > 0 and sumstats.exists() and sumstats.stat().st_size > 0


def pending_jobs(rows: list[dict[str, str]], force: bool) -> list[dict[str, object]]:
    jobs: list[dict[str, object]] = []
    for row in rows:
        pheno_id = row["pheno_id"]
        pheno_path = Path(row["pheno_path"])
        if not pheno_path.exists() or pheno_path.stat().st_size == 0:
            raise FileNotFoundError(f"missing phenotype file for {pheno_id}: {pheno_path}")
        pheno_name = row.get("pheno_name") or pheno_name_from_file(pheno_path)
        glm = Path(row["glm"])
        sumstats = Path(row["sumstats"])
        if not force and output_complete(row):
            continue
        jobs.append({
            "pheno_id": pheno_id,
            "pheno_name": pheno_name,
            "pheno_path": pheno_path,
            "glm": glm,
            "sumstats": sumstats,
        })
    return jobs


def write_batch_manifest(path: Path, batch_jobs: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["pheno_id", "pheno_name", "glm", "sumstats"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for job in batch_jobs:
            writer.writerow({field: str(job[field]) for field in fields})


def write_duplicate_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        if path.exists():
            path.unlink()
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["duplicate_key"]
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.batch_size < 1:
        raise SystemExit("--batch-size must be >= 1")

    manifest_rows = read_manifest(args.manifest)
    rows, duplicate_rows = dedupe_manifest_rows(manifest_rows)
    jobs = pending_jobs(rows, args.force)
    args.workdir.mkdir(parents=True, exist_ok=True)
    write_duplicate_manifest(args.workdir / "duplicate_manifest_rows.tsv", duplicate_rows)

    for old in args.workdir.glob("batch_*.pheno.tsv"):
        old.unlink()
    for old in args.workdir.glob("batch_*.keep.tsv"):
        old.unlink()
    for old in args.workdir.glob("batch_*.manifest.tsv"):
        old.unlink()

    plan_path = args.workdir / "batch_plan.tsv"
    fields = ["batch_index", "pheno_tsv", "keep_tsv", "manifest_tsv", "n_phenotypes"]
    with open(plan_path, "w", newline="") as plan:
        writer = csv.DictWriter(plan, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i in range(0, len(jobs), args.batch_size):
            batch_index = i // args.batch_size + 1
            batch = jobs[i:i + args.batch_size]
            stem = args.workdir / f"batch_{batch_index:05d}"
            pheno_tsv = stem.with_suffix(".pheno.tsv")
            keep_tsv = stem.with_suffix(".keep.tsv")
            manifest_tsv = stem.with_suffix(".manifest.tsv")
            write_batch_pheno(batch, pheno_tsv, keep_tsv)
            write_batch_manifest(manifest_tsv, batch)
            writer.writerow({
                "batch_index": batch_index,
                "pheno_tsv": pheno_tsv,
                "keep_tsv": keep_tsv,
                "manifest_tsv": manifest_tsv,
                "n_phenotypes": len(batch),
            })

    summary_path = args.workdir / "batch_plan.summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["manifest_rows", len(manifest_rows)])
        writer.writerow(["duplicate_manifest_rows", len(duplicate_rows)])
        writer.writerow(["manifest_phenotypes", len(rows)])
        writer.writerow(["pending_phenotypes", len(jobs)])
        writer.writerow(["batch_size", args.batch_size])
        writer.writerow(["batches", math.ceil(len(jobs) / args.batch_size) if jobs else 0])

    print(
        f"manifest_rows={len(manifest_rows)} duplicate_manifest_rows={len(duplicate_rows)} "
        f"manifest_phenotypes={len(rows)} pending_phenotypes={len(jobs)} "
        f"batches={math.ceil(len(jobs) / args.batch_size) if jobs else 0}"
    )


if __name__ == "__main__":
    main()
