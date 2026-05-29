"""Select unrelated European IIDs for fitting PCA.

This AoU version mirrors the UKBB PCA-training selector, but replaces the
UKBB White British seed set with AoU samples that this pipeline classified as
European and that are involved in sibling or identical/twin/duplicate pairs.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


def sort_key(iid: str) -> tuple[int, int | str]:
    try:
        return (0, int(iid))
    except ValueError:
        return (1, iid)


def parse_seed_relationships(text: str) -> set[str]:
    values = {part.strip().lower() for part in text.split(",") if part.strip()}
    if not values:
        raise ValueError("--seed-relationships cannot be empty")
    return values


def read_iids(path: Path) -> tuple[set[str], dict[str, str]]:
    iids: set[str] = set()
    fid_by_iid: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                fid, iid = parts[0], parts[1]
            else:
                fid = iid = parts[0]
            iids.add(iid)
            fid_by_iid.setdefault(iid, fid)
    return iids, fid_by_iid


def write_id_file(path: Path, iids: set[str], fid_by_iid: dict[str, str]) -> None:
    with path.open("w") as out:
        for iid in sorted(iids, key=sort_key):
            out.write(f"{fid_by_iid.get(iid, '0')}\t{iid}\n")


def write_psam(path: Path, iids: set[str], fid_by_iid: dict[str, str]) -> None:
    with path.open("w") as out:
        out.write("#FID\tIID\n")
        for iid in sorted(iids, key=sort_key):
            out.write(f"{fid_by_iid.get(iid, '0')}\t{iid}\n")


def read_close_relation_seeds(
    path: Path,
    european_iids: set[str],
    seed_relationships: set[str],
) -> tuple[set[str], set[str], int]:
    all_seed_relation_samples: set[str] = set()
    european_seed_samples: set[str] = set()
    seed_rows = 0

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"eid1", "eid2", "relationship"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")

        for row in reader:
            relationship = row["relationship"].strip().lower()
            if relationship not in seed_relationships:
                continue
            seed_rows += 1
            pair = {row["eid1"].strip(), row["eid2"].strip()}
            all_seed_relation_samples.update(pair)
            european_seed_samples.update(pair & european_iids)

    return all_seed_relation_samples, european_seed_samples, seed_rows


def read_kinship_edges(
    path: Path,
    threshold: float,
    fid_by_iid: dict[str, str],
) -> tuple[list[tuple[str, str, float]], dict[str, set[str]]]:
    edges: list[tuple[str, str, float]] = []
    graph: dict[str, set[str]] = defaultdict(set)

    with path.open() as handle:
        header = handle.readline().strip().split()
        cols = {name.lstrip("#"): idx for idx, name in enumerate(header)}
        required = {"FID1", "IID1", "FID2", "IID2", "KINSHIP"}
        missing = required.difference(cols)
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")

        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            iid1 = parts[cols["IID1"]]
            iid2 = parts[cols["IID2"]]
            fid_by_iid.setdefault(iid1, parts[cols["FID1"]])
            fid_by_iid.setdefault(iid2, parts[cols["FID2"]])
            kinship = float(parts[cols["KINSHIP"]])
            if kinship >= threshold:
                edges.append((iid1, iid2, kinship))
                graph[iid1].add(iid2)
                graph[iid2].add(iid1)

    return edges, graph


def read_plink_in_ids(path: Path) -> set[str]:
    iids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                iids.add(parts[1])
            else:
                iids.add(parts[0])
    return iids


def write_summary(path: Path, rows: list[tuple[str, str | int | float]]) -> None:
    with path.open("w") as out:
        out.write("metric\tvalue\n")
        for key, value in rows:
            out.write(f"{key}\t{value}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--europeans", type=Path, required=True)
    parser.add_argument("--close-relations", type=Path, required=True)
    parser.add_argument("--kin0", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plink2", type=Path, required=True)
    parser.add_argument("--kinship-threshold", type=float, default=0.0441941)
    parser.add_argument("--seed-relationships", default="sibling,identical")
    args = parser.parse_args()

    seed_relationships = parse_seed_relationships(args.seed_relationships)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scrap_dir = args.output_dir / "scrap"
    scrap_dir.mkdir(parents=True, exist_ok=True)

    log_lines: list[str] = []

    def log(message: str) -> None:
        print(message, flush=True)
        log_lines.append(message)

    log("=== Select PCA European IIDs ===")
    log(f"European keep-list: {args.europeans}")
    log(f"Close relationships: {args.close_relations}")
    log(f"KING table: {args.kin0}")
    log(f"KINSHIP threshold: {args.kinship_threshold}")
    log(f"Seed relationships: {','.join(sorted(seed_relationships))}")

    european_iids, fid_by_iid = read_iids(args.europeans)
    log(f"Classified Europeans: {len(european_iids)}")
    if not european_iids:
        raise ValueError("European keep-list is empty")

    all_seed_relation_samples, european_seed_samples, seed_rows = read_close_relation_seeds(
        args.close_relations,
        european_iids,
        seed_relationships,
    )
    log(f"Seed relationship rows: {seed_rows}")
    log(f"Samples in seed relationship rows: {len(all_seed_relation_samples)}")
    log(f"European seed samples: {len(european_seed_samples)}")

    edges, kinship_graph = read_kinship_edges(args.kin0, args.kinship_threshold, fid_by_iid)
    log(
        "Kinship graph: "
        f"{len(kinship_graph)} nodes, {len(edges)} edges "
        f"(KINSHIP >= {args.kinship_threshold})"
    )

    expanded_exclusions = set(european_seed_samples)
    for iid in european_seed_samples:
        expanded_exclusions.update(kinship_graph.get(iid, set()))
    expanded_exclusions_in_eur = expanded_exclusions & european_iids
    candidate_pca_europeans = european_iids - expanded_exclusions
    if not candidate_pca_europeans:
        raise ValueError("No candidate PCA Europeans remain after seed-relative exclusion")

    log("")
    log("=== Sample filtering ===")
    log(f"Expanded exclusions: {len(expanded_exclusions)}")
    log(f"Expanded exclusions in Europeans: {len(expanded_exclusions_in_eur)}")
    log(f"Candidate PCA Europeans before KING cutoff: {len(candidate_pca_europeans)}")

    seed_path = args.output_dir / "seed_sibling_identical_iids.txt"
    expanded_path = args.output_dir / "expanded_sibling_identical_relatives_iids.txt"
    candidate_path = args.output_dir / "candidate_pca_europeans_iids.txt"
    candidate_psam = scrap_dir / "candidate_pca_europeans.psam"
    write_id_file(seed_path, european_seed_samples, fid_by_iid)
    write_id_file(expanded_path, expanded_exclusions, fid_by_iid)
    write_id_file(candidate_path, candidate_pca_europeans, fid_by_iid)
    write_psam(candidate_psam, candidate_pca_europeans, fid_by_iid)

    log("")
    log("=== KING cutoff ===")
    out_prefix = scrap_dir / "fit_pca"
    cmd = [
        str(args.plink2),
        "--psam",
        str(candidate_psam),
        "--king-cutoff-table",
        str(args.kin0),
        str(args.kinship_threshold),
        "--out",
        str(out_prefix),
    ]
    log("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)

    plink_in = out_prefix.with_suffix(".king.cutoff.in.id")
    plink_out = out_prefix.with_suffix(".king.cutoff.out.id")
    if not plink_in.exists():
        raise FileNotFoundError(plink_in)
    fit_pca_iids = read_plink_in_ids(plink_in)
    removed_by_king = candidate_pca_europeans - fit_pca_iids

    final_path = args.output_dir / "fit_pca_iids.txt"
    write_id_file(final_path, fit_pca_iids, fid_by_iid)

    log("")
    log("=== Verification ===")
    not_european = fit_pca_iids - european_iids
    in_expanded = fit_pca_iids & expanded_exclusions
    related_edges_in_final = [
        (iid1, iid2)
        for iid1, iid2, _kinship in edges
        if iid1 in fit_pca_iids and iid2 in fit_pca_iids
    ]
    log(f"Final fit_pca_iids: {len(fit_pca_iids)}")
    log(f"Removed by KING cutoff: {len(removed_by_king)}")
    log(f"Verification not European: {len(not_european)}")
    log(f"Verification in expanded exclusions: {len(in_expanded)}")
    log(f"Verification related pairs in final set: {len(related_edges_in_final)}")

    if not_european:
        raise ValueError("Verification failed: final fit_pca_iids contains non-European IIDs")
    if in_expanded:
        raise ValueError("Verification failed: final fit_pca_iids contains expanded exclusions")
    if related_edges_in_final:
        raise ValueError(
            "Verification failed: final fit_pca_iids contains pairs with "
            f"KINSHIP >= {args.kinship_threshold}"
        )
    if len(fit_pca_iids) != sum(1 for _line in final_path.open()):
        raise ValueError("Verification failed: final output line count mismatch")

    summary_rows: list[tuple[str, str | int | float]] = [
        ("classified_europeans", len(european_iids)),
        ("seed_relationships", ",".join(sorted(seed_relationships))),
        ("seed_relationship_rows", seed_rows),
        ("seed_relationship_samples_total", len(all_seed_relation_samples)),
        ("seed_relationship_europeans", len(european_seed_samples)),
        ("kinship_threshold", args.kinship_threshold),
        ("kinship_edges_ge_threshold", len(edges)),
        ("kinship_graph_nodes", len(kinship_graph)),
        ("expanded_exclusions_total", len(expanded_exclusions)),
        ("expanded_exclusions_in_europeans", len(expanded_exclusions_in_eur)),
        ("expanded_exclusions_not_in_europeans", len(expanded_exclusions - european_iids)),
        ("candidate_pca_europeans_before_king_cutoff", len(candidate_pca_europeans)),
        ("removed_by_king_cutoff", len(removed_by_king)),
        ("fit_pca_iids", len(fit_pca_iids)),
        ("verification_not_european", len(not_european)),
        ("verification_in_expanded_exclusions", len(in_expanded)),
        ("verification_related_pairs_ge_threshold", len(related_edges_in_final)),
        ("plink_in_file", str(plink_in)),
        ("plink_out_file", str(plink_out)),
    ]
    summary_path = args.output_dir / "select_pca_europeans.summary.tsv"
    write_summary(summary_path, summary_rows)

    log("")
    log("=== Outputs ===")
    log(f"Wrote {final_path}")
    log(f"Wrote {summary_path}")
    log(f"Wrote {seed_path}")
    log(f"Wrote {expanded_path}")
    log(f"Wrote {candidate_path}")

    log_path = args.output_dir / "pca_eur_log.txt"
    with log_path.open("w") as out:
        out.write("\n".join(log_lines) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
