#!/usr/bin/env python3
"""Build AoU height GWAS phenotype/covariate files."""

import argparse
import csv
import os
import statistics
import sys


def log(lines, msg):
    print(msg)
    lines.append(msg)


def sort_iids(iids):
    return sorted(iids, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))


def read_keep_iids(path):
    iids = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                iids.add(parts[-1])
    return iids


def load_height_rows(path):
    rows = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"IID", "height", "age_at_height", "n_height_records"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["IID"].strip()
            if not iid:
                continue
            height = float(row["height"])
            age = float(row["age_at_height"])
            n_height_records = int(row["n_height_records"])
            rows[iid] = {
                "height": height,
                "age_at_height": age,
                "n_height_records": n_height_records,
            }
    return rows


def load_sex_covar(path):
    sex_map = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"FID", "IID", "sex_01"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["IID"].strip()
            if not iid:
                continue
            sex_01 = int(row["sex_01"])
            if sex_01 not in (0, 1):
                raise ValueError(f"Unexpected sex_01={sex_01} for IID {iid}")
            sex_map[iid] = sex_01
    return sex_map


def load_fam_fids(path):
    fid_by_iid = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            fid_by_iid[parts[1]] = parts[0]
    return fid_by_iid


def load_projected_pcs(path, n_pcs):
    pc_data = {}
    pc_headers = [f"PC{i}_AVG" for i in range(1, n_pcs + 1)]
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header and header[0] == "#FID":
            header[0] = "FID"
        missing = [h for h in ["IID"] + pc_headers if h not in header]
        if missing:
            raise ValueError(f"{path} missing columns: {missing}")
        iid_idx = header.index("IID")
        pc_indices = [header.index(h) for h in pc_headers]
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(pc_indices):
                continue
            iid = fields[iid_idx]
            pc_data[iid] = [fields[i] for i in pc_indices]
    return pc_data, pc_headers


def read_iids_from_file(path, has_header, sep):
    out = []
    with open(path) as f:
        if has_header:
            f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            out.append(line.split(sep)[1])
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height-query", required=True)
    parser.add_argument("--europeans", required=True)
    parser.add_argument("--sex-covar", required=True)
    parser.add_argument("--fam", required=True)
    parser.add_argument("--sscore", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--height-min-cm", type=float, default=140.0)
    parser.add_argument("--n-pcs", type=int, default=10)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    log_lines = []
    log(log_lines, "=== AoU height GWAS setup ===")
    log(log_lines, f"height_query: {args.height_query}")
    log(log_lines, f"height_min_cm: {args.height_min_cm}")
    log(log_lines, f"n_pcs: {args.n_pcs}")

    europeans = read_keep_iids(args.europeans)
    sex_map = load_sex_covar(args.sex_covar)
    fid_by_iid = load_fam_fids(args.fam)
    height_rows = load_height_rows(args.height_query)
    pc_data, pc_headers = load_projected_pcs(args.sscore, args.n_pcs)

    log(log_lines, "")
    log(log_lines, "=== Input counts ===")
    log(log_lines, f"classified Europeans: {len(europeans)}")
    log(log_lines, f"sex covariate rows: {len(sex_map)}")
    log(log_lines, f"fam rows: {len(fid_by_iid)}")
    log(log_lines, f"height query rows: {len(height_rows)}")
    log(log_lines, f"projected PC rows: {len(pc_data)}")

    candidates = set(europeans)
    missing_height = candidates - set(height_rows)
    candidates &= set(height_rows)
    missing_sex = candidates - set(sex_map)
    candidates &= set(sex_map)
    missing_fam = candidates - set(fid_by_iid)
    candidates &= set(fid_by_iid)
    missing_pcs = candidates - set(pc_data)
    candidates &= set(pc_data)

    gwas_iids = candidates
    gwas_iids = sort_iids(gwas_iids)
    if not gwas_iids:
        raise RuntimeError("No GWAS samples remain after height/sex/PC filters")

    ages = [height_rows[iid]["age_at_height"] for iid in gwas_iids]
    mean_age = statistics.mean(ages)

    covar_data = {}
    for iid in gwas_iids:
        age_c = height_rows[iid]["age_at_height"] - mean_age
        sex_c = sex_map[iid] - 0.5
        covar_data[iid] = (age_c, sex_c, age_c * sex_c)

    males = sum(1 for iid in gwas_iids if sex_map[iid] == 1)
    females = sum(1 for iid in gwas_iids if sex_map[iid] == 0)
    heights = [height_rows[iid]["height"] for iid in gwas_iids]

    log(log_lines, "")
    log(log_lines, "=== Filtering counts ===")
    log(log_lines, f"Europeans missing height query row: {len(missing_height)}")
    log(log_lines, f"Europeans with height row but missing sex covariate: {len(missing_sex)}")
    log(log_lines, f"Europeans with height+sex but missing fam row: {len(missing_fam)}")
    log(log_lines, f"Europeans with height+sex+fam but missing projected PCs: {len(missing_pcs)}")
    log(log_lines, f"Final GWAS samples: {len(gwas_iids)}")
    log(log_lines, f"Final female/male: {females}/{males}")
    log(log_lines, f"Mean age at height measurement: {mean_age:.6f}")
    log(log_lines, f"Mean height: {statistics.mean(heights):.6f}")
    log(log_lines, f"Median height: {statistics.median(heights):.6f}")

    training_path = os.path.join(args.out_dir, "training_iids.txt")
    phen_path = os.path.join(args.out_dir, "phen.txt")
    base_covar_path = os.path.join(args.out_dir, "base_covar.txt")
    covar_path = os.path.join(args.out_dir, "covar.txt")
    log_path = os.path.join(args.out_dir, "height_gwas_log.txt")
    summary_path = os.path.join(args.out_dir, "height_gwas.summary.tsv")

    with open(training_path, "w") as f:
        for iid in gwas_iids:
            f.write(f"{fid_by_iid[iid]} {iid}\n")

    with open(phen_path, "w") as f:
        f.write("FID\tIID\theight\n")
        for iid in gwas_iids:
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{height_rows[iid]['height']:.10g}\n")

    with open(base_covar_path, "w") as f:
        f.write("FID\tIID\tage_c\tsex_c\tage_c_sex_c_inter\n")
        for iid in gwas_iids:
            age_c, sex_c, inter = covar_data[iid]
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{age_c:.12g}\t{sex_c:.1f}\t{inter:.12g}\n")

    with open(covar_path, "w") as f:
        f.write("FID\tIID\tage_c\tsex_c\tage_c_sex_c_inter\t" + "\t".join(pc_headers) + "\n")
        for iid in gwas_iids:
            age_c, sex_c, inter = covar_data[iid]
            pcs = "\t".join(pc_data[iid])
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{age_c:.12g}\t{sex_c:.1f}\t{inter:.12g}\t{pcs}\n")

    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"classified_europeans\t{len(europeans)}\n")
        f.write(f"sex_covar_rows\t{len(sex_map)}\n")
        f.write(f"fam_rows\t{len(fid_by_iid)}\n")
        f.write(f"height_query_rows\t{len(height_rows)}\n")
        f.write(f"projected_pc_rows\t{len(pc_data)}\n")
        f.write(f"europeans_missing_height_query\t{len(missing_height)}\n")
        f.write(f"height_candidates_missing_sex_covar\t{len(missing_sex)}\n")
        f.write(f"height_sex_candidates_missing_fam\t{len(missing_fam)}\n")
        f.write(f"height_sex_fam_candidates_missing_pcs\t{len(missing_pcs)}\n")
        f.write(f"gwas_samples\t{len(gwas_iids)}\n")
        f.write(f"gwas_female\t{females}\n")
        f.write(f"gwas_male\t{males}\n")
        f.write(f"height_min_cm\t{args.height_min_cm}\n")
        f.write(f"height_mean_cm\t{statistics.mean(heights):.10g}\n")
        f.write(f"height_median_cm\t{statistics.median(heights):.10g}\n")
        f.write(f"age_mean\t{mean_age:.10g}\n")
        f.write(f"n_pcs\t{args.n_pcs}\n")

    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")

    print("\n=== Verification checks ===")
    def check(name, condition, detail=""):
        nonlocal_passed[0] = nonlocal_passed[0] and bool(condition)
        print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))

    nonlocal_passed = [True]
    training_iids = read_iids_from_file(training_path, False, " ")
    phen_iids = read_iids_from_file(phen_path, True, "\t")
    base_iids = read_iids_from_file(base_covar_path, True, "\t")
    covar_iids = read_iids_from_file(covar_path, True, "\t")

    check("same IIDs in all output files", training_iids == phen_iids == base_iids == covar_iids)
    check("all GWAS IIDs are classified European", set(training_iids) <= europeans)
    check("all GWAS IIDs have sex covariate", set(training_iids) <= set(sex_map))
    check("all GWAS IIDs have genotype FID", set(training_iids) <= set(fid_by_iid))
    check("all output FID/IID pairs match genotype fam",
          all(line.split()[0] == fid_by_iid[line.split()[1]]
              for line in open(training_path) if line.strip()))
    check("row count matches summary", len(training_iids) == len(gwas_iids))
    check("all heights >= threshold", all(height_rows[iid]["height"] >= args.height_min_cm for iid in gwas_iids))
    check("mean age_c is approximately zero",
          abs(statistics.mean(covar_data[iid][0] for iid in gwas_iids)) < 1e-8)
    check("sex_c values are {-0.5, 0.5}", {covar_data[iid][1] for iid in gwas_iids} == {-0.5, 0.5})
    check("all GWAS IIDs have requested PCs", set(training_iids) <= set(pc_data))
    if not nonlocal_passed[0]:
        sys.exit(1)

    print("\nAll verification checks passed.")


if __name__ == "__main__":
    main()
