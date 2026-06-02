#!/usr/bin/env python3
"""Build sample-QC exclusions from identical-genotype relationship components."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def sort_iids(iids):
    return sorted(iids, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))


def read_fam(path):
    fid_by_iid = {}
    with path.open() as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                fid_by_iid[parts[1]] = parts[0]
    return fid_by_iid


def connected_components(graph):
    seen = set()
    out = []
    for node in sort_iids(graph):
        if node in seen:
            continue
        stack = [node]
        seen.add(node)
        comp = []
        while stack:
            current = stack.pop()
            comp.append(current)
            for neighbor in graph[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        out.append(sort_iids(comp))
    return out


def edge_count_for_component(graph, component):
    nodes = set(component)
    return sum(1 for node in nodes for neighbor in graph[node] if neighbor in nodes) // 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--close-relations", type=Path, required=True)
    parser.add_argument("--fam", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--exclude-min-component-size", type=int, default=3)
    args = parser.parse_args()

    if args.exclude_min_component_size < 2:
        raise ValueError("--exclude-min-component-size must be >= 2")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fid_by_iid = read_fam(args.fam)

    graph = defaultdict(set)
    total_close_rows = 0
    total_identical_pairs = 0
    with args.close_relations.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"eid1", "eid2", "relationship"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{args.close_relations} missing columns: {sorted(missing)}")
        for row in reader:
            total_close_rows += 1
            if row["relationship"].strip() != "identical":
                continue
            iid1 = row["eid1"].strip()
            iid2 = row["eid2"].strip()
            if not iid1 or not iid2 or iid1 == iid2:
                continue
            graph[iid1].add(iid2)
            graph[iid2].add(iid1)
            total_identical_pairs += 1

    components = connected_components(graph)
    component_records = []
    excluded_iids = set()
    size_counts = Counter()
    complete_counts = Counter()

    for idx, component in enumerate(components, start=1):
        component_size = len(component)
        n_edges = edge_count_for_component(graph, component)
        expected_edges = component_size * (component_size - 1) // 2
        is_complete_clique = n_edges == expected_edges
        excluded = component_size >= args.exclude_min_component_size
        size_counts[component_size] += 1
        complete_counts[is_complete_clique] += 1
        if excluded:
            excluded_iids.update(component)
        component_records.append(
            {
                "component_id": f"identical_component_{idx:06d}",
                "component_size": component_size,
                "n_identical_pairs": n_edges,
                "expected_complete_pairs": expected_edges,
                "is_complete_clique": is_complete_clique,
                "exclude_component": excluded,
                "members": component,
            }
        )

    missing_fam_iids = sorted(excluded_iids.difference(fid_by_iid))
    if missing_fam_iids:
        preview = ", ".join(missing_fam_iids[:10])
        raise ValueError(
            f"{len(missing_fam_iids)} excluded IIDs are missing from {args.fam}; first: {preview}"
        )

    members_path = args.out_dir / "identical_components.tsv"
    component_summary_path = args.out_dir / "identical_component_summary.tsv"
    exclude_path = (
        args.out_dir
        / f"exclude_identical_component_size_ge{args.exclude_min_component_size}_iids.txt"
    )
    summary_path = args.out_dir / "identical_component_sample_qc.summary.tsv"
    log_path = args.out_dir / "identical_component_sample_qc.log"

    with members_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "component_id",
                "IID",
                "FID",
                "component_size",
                "n_identical_pairs",
                "expected_complete_pairs",
                "is_complete_clique",
                "exclude_component",
            ]
        )
        for record in component_records:
            for iid in record["members"]:
                writer.writerow(
                    [
                        record["component_id"],
                        iid,
                        fid_by_iid.get(iid, ""),
                        record["component_size"],
                        record["n_identical_pairs"],
                        record["expected_complete_pairs"],
                        int(record["is_complete_clique"]),
                        int(record["exclude_component"]),
                    ]
                )

    with component_summary_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "component_id",
                "component_size",
                "n_identical_pairs",
                "expected_complete_pairs",
                "is_complete_clique",
                "exclude_component",
            ]
        )
        for record in component_records:
            writer.writerow(
                [
                    record["component_id"],
                    record["component_size"],
                    record["n_identical_pairs"],
                    record["expected_complete_pairs"],
                    int(record["is_complete_clique"]),
                    int(record["exclude_component"]),
                ]
            )

    with exclude_path.open("w") as f:
        for iid in sort_iids(excluded_iids):
            f.write(f"{fid_by_iid[iid]} {iid}\n")

    summary_rows = [
        ("total_close_relationship_rows", total_close_rows),
        ("total_identical_pairs", total_identical_pairs),
        ("unique_iids_in_identical_graph", len(graph)),
        ("identical_components", len(components)),
        ("exclude_min_component_size", args.exclude_min_component_size),
        (
            f"components_size_ge_{args.exclude_min_component_size}",
            sum(1 for r in component_records if r["exclude_component"]),
        ),
        (
            f"iids_in_components_size_ge_{args.exclude_min_component_size}",
            len(excluded_iids),
        ),
        ("complete_clique_components", complete_counts[True]),
        ("noncomplete_components", complete_counts[False]),
        ("missing_fam_iids_among_excluded", len(missing_fam_iids)),
    ]
    for size in sorted(size_counts):
        summary_rows.append((f"component_size_{size}_count", size_counts[size]))

    with summary_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary_rows)

    log_lines = [
        "=== Identical-genotype component sample QC ===",
        f"close_relations: {args.close_relations}",
        f"fam: {args.fam}",
        f"exclude_min_component_size: {args.exclude_min_component_size}",
        "",
        "=== Summary ===",
    ]
    log_lines.extend(f"{k}: {v}" for k, v in summary_rows)
    log_lines.extend(
        [
            "",
            f"wrote {members_path}",
            f"wrote {component_summary_path}",
            f"wrote {exclude_path}",
            f"wrote {summary_path}",
        ]
    )
    log_path.write_text("\n".join(log_lines) + "\n")

    for line in log_lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
