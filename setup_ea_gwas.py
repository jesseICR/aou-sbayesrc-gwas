#!/usr/bin/env python3
"""Build AoU educational-attainment GWAS phenotype/covariate files."""

import argparse
import csv
import os
import statistics
import sys
from collections import Counter


EA_MAPPING = {
    1585941: 1.0,    # Never attended
    1585942: 2.5,    # Grades 1-4 midpoint
    1585943: 6.5,    # Grades 5-8 midpoint
    1585944: 10.0,   # Some high school
    1585945: 13.0,   # Twelve or GED
    1585946: 15.0,   # College one to three
    1585947: 17.0,   # College graduate
    1585948: 20.0,   # Advanced degree
}


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


def load_ea_rows(path):
    rows = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"IID", "ea_years", "yob", "age_at_survey", "answer_concept_id", "answer", "n_ea_records"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["IID"].strip()
            if not iid:
                continue
            answer_concept_id = int(row["answer_concept_id"])
            ea_years = float(row["ea_years"])
            if answer_concept_id not in EA_MAPPING:
                raise ValueError(f"Unexpected EA answer_concept_id={answer_concept_id} for IID {iid}")
            if abs(ea_years - EA_MAPPING[answer_concept_id]) > 1e-8:
                raise ValueError(f"EA mapping mismatch for answer_concept_id={answer_concept_id}")
            rows[iid] = {
                "ea_years": ea_years,
                "yob": float(row["yob"]),
                "age_at_survey": float(row["age_at_survey"]),
                "answer_concept_id": answer_concept_id,
                "answer": row["answer"],
                "n_ea_records": int(row["n_ea_records"]),
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
            pc_data[fields[iid_idx]] = [fields[i] for i in pc_indices]
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


def write_answer_counts(path, ea_rows, gwas_iids):
    all_counts = Counter(row["answer_concept_id"] for row in ea_rows.values())
    final_counts = Counter(ea_rows[iid]["answer_concept_id"] for iid in gwas_iids)
    answer_names = {}
    for row in ea_rows.values():
        answer_names[row["answer_concept_id"]] = row["answer"]

    with open(path, "w") as f:
        f.write("answer_concept_id\tanswer\tea_years\tquery_count\tfinal_gwas_count\n")
        for concept_id in sorted(EA_MAPPING, key=lambda x: EA_MAPPING[x]):
            f.write(
                f"{concept_id}\t{answer_names.get(concept_id, '')}\t"
                f"{EA_MAPPING[concept_id]:.10g}\t{all_counts[concept_id]}\t{final_counts[concept_id]}\n"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ea-query", required=True)
    parser.add_argument("--europeans", required=True)
    parser.add_argument("--sex-covar", required=True)
    parser.add_argument("--exclude-iids", required=True)
    parser.add_argument("--fam", required=True)
    parser.add_argument("--sscore", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-pcs", type=int, default=10)
    parser.add_argument("--min-age-at-survey", type=float, default=26.0)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    log_lines = []
    log(log_lines, "=== AoU educational-attainment GWAS setup ===")
    log(log_lines, f"ea_query: {args.ea_query}")
    log(log_lines, f"exclude_iids: {args.exclude_iids}")
    log(log_lines, f"n_pcs: {args.n_pcs}")
    log(log_lines, f"min_age_at_survey: {args.min_age_at_survey:g}")

    europeans = read_keep_iids(args.europeans)
    sex_map = load_sex_covar(args.sex_covar)
    excluded_iids = read_keep_iids(args.exclude_iids)
    fid_by_iid = load_fam_fids(args.fam)
    ea_rows = load_ea_rows(args.ea_query)
    pc_data, pc_headers = load_projected_pcs(args.sscore, args.n_pcs)

    log(log_lines, "")
    log(log_lines, "=== Input counts ===")
    log(log_lines, f"classified Europeans: {len(europeans)}")
    log(log_lines, f"sex covariate rows: {len(sex_map)}")
    log(log_lines, f"sample-QC exclusion IIDs: {len(excluded_iids)}")
    log(log_lines, f"fam rows: {len(fid_by_iid)}")
    log(log_lines, f"EA query rows: {len(ea_rows)}")
    log(log_lines, f"projected PC rows: {len(pc_data)}")

    duplicate_record_iids = {iid for iid, row in ea_rows.items() if row["n_ea_records"] > 1}
    excluded_europeans = europeans & excluded_iids
    candidates = set(europeans) - excluded_iids
    after_sample_qc = set(candidates)
    missing_ea = candidates - set(ea_rows)
    candidates &= set(ea_rows)
    below_min_age = {iid for iid in candidates if ea_rows[iid]["age_at_survey"] < args.min_age_at_survey}
    candidates -= below_min_age
    missing_sex = candidates - set(sex_map)
    candidates &= set(sex_map)
    missing_fam = candidates - set(fid_by_iid)
    candidates &= set(fid_by_iid)
    missing_pcs = candidates - set(pc_data)
    candidates &= set(pc_data)

    gwas_iids = sort_iids(candidates)
    if not gwas_iids:
        raise RuntimeError("No EA GWAS samples remain after phenotype/sex/PC filters")

    yobs = [ea_rows[iid]["yob"] for iid in gwas_iids]
    ages = [ea_rows[iid]["age_at_survey"] for iid in gwas_iids]
    mean_yob = statistics.mean(yobs)
    covar_data = {}
    for iid in gwas_iids:
        yob_c = ea_rows[iid]["yob"] - mean_yob
        sex_c = sex_map[iid] - 0.5
        covar_data[iid] = (yob_c, sex_c, yob_c * sex_c)

    males = sum(1 for iid in gwas_iids if sex_map[iid] == 1)
    females = sum(1 for iid in gwas_iids if sex_map[iid] == 0)
    ea_values = [ea_rows[iid]["ea_years"] for iid in gwas_iids]

    log(log_lines, "")
    log(log_lines, "=== Filtering counts ===")
    log(log_lines, f"Europeans removed by sample-QC exclusion: {len(excluded_europeans)}")
    log(log_lines, f"Europeans after sample-QC exclusion: {len(after_sample_qc)}")
    log(log_lines, f"Europeans missing codeable EA phenotype: {len(missing_ea)}")
    log(log_lines, f"Europeans with EA phenotype but age_at_survey < {args.min_age_at_survey:g}: {len(below_min_age)}")
    log(log_lines, f"Europeans with EA phenotype but missing sex covariate: {len(missing_sex)}")
    log(log_lines, f"Europeans with EA+sex but missing fam row: {len(missing_fam)}")
    log(log_lines, f"Europeans with EA+sex+fam but missing projected PCs: {len(missing_pcs)}")
    log(log_lines, f"EA query IIDs with multiple codeable records: {len(duplicate_record_iids)}")
    log(log_lines, f"Final GWAS samples: {len(gwas_iids)}")
    log(log_lines, f"Final female/male: {females}/{males}")
    log(log_lines, f"Mean fractional YOB: {mean_yob:.6f}")
    log(log_lines, f"Mean age at EA survey: {statistics.mean(ages):.6f}")
    log(log_lines, f"Minimum age at EA survey: {min(ages):.6f}")
    log(log_lines, f"Mean EA years: {statistics.mean(ea_values):.6f}")
    log(log_lines, f"Median EA years: {statistics.median(ea_values):.6f}")

    training_path = os.path.join(args.out_dir, "training_iids.txt")
    phen_path = os.path.join(args.out_dir, "phen.txt")
    base_covar_path = os.path.join(args.out_dir, "base_covar.txt")
    covar_path = os.path.join(args.out_dir, "covar.txt")
    answer_counts_path = os.path.join(args.out_dir, "ea_answer_counts.tsv")
    log_path = os.path.join(args.out_dir, "ea_gwas_log.txt")
    summary_path = os.path.join(args.out_dir, "ea_gwas.summary.tsv")

    with open(training_path, "w") as f:
        for iid in gwas_iids:
            f.write(f"{fid_by_iid[iid]} {iid}\n")

    with open(phen_path, "w") as f:
        f.write("FID\tIID\tea_years\n")
        for iid in gwas_iids:
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{ea_rows[iid]['ea_years']:.10g}\n")

    with open(base_covar_path, "w") as f:
        f.write("FID\tIID\tyob_c\tsex_c\tyob_c_sex_c_inter\n")
        for iid in gwas_iids:
            yob_c, sex_c, inter = covar_data[iid]
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{yob_c:.12g}\t{sex_c:.1f}\t{inter:.12g}\n")

    with open(covar_path, "w") as f:
        f.write("FID\tIID\tyob_c\tsex_c\tyob_c_sex_c_inter\t" + "\t".join(pc_headers) + "\n")
        for iid in gwas_iids:
            yob_c, sex_c, inter = covar_data[iid]
            pcs = "\t".join(pc_data[iid])
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{yob_c:.12g}\t{sex_c:.1f}\t{inter:.12g}\t{pcs}\n")

    write_answer_counts(answer_counts_path, ea_rows, gwas_iids)

    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"classified_europeans\t{len(europeans)}\n")
        f.write(f"sex_covar_rows\t{len(sex_map)}\n")
        f.write(f"sample_qc_exclusion_iids\t{len(excluded_iids)}\n")
        f.write(f"classified_europeans_removed_by_sample_qc\t{len(excluded_europeans)}\n")
        f.write(f"classified_europeans_after_sample_qc\t{len(after_sample_qc)}\n")
        f.write(f"fam_rows\t{len(fid_by_iid)}\n")
        f.write(f"ea_query_rows\t{len(ea_rows)}\n")
        f.write(f"projected_pc_rows\t{len(pc_data)}\n")
        f.write(f"europeans_missing_codeable_ea\t{len(missing_ea)}\n")
        f.write(f"min_age_at_survey\t{args.min_age_at_survey:.10g}\n")
        f.write(f"ea_candidates_below_min_age_at_survey\t{len(below_min_age)}\n")
        f.write(f"ea_candidates_missing_sex_covar\t{len(missing_sex)}\n")
        f.write(f"ea_sex_candidates_missing_fam\t{len(missing_fam)}\n")
        f.write(f"ea_sex_fam_candidates_missing_pcs\t{len(missing_pcs)}\n")
        f.write(f"ea_query_iids_with_multiple_codeable_records\t{len(duplicate_record_iids)}\n")
        f.write(f"gwas_samples\t{len(gwas_iids)}\n")
        f.write(f"gwas_female\t{females}\n")
        f.write(f"gwas_male\t{males}\n")
        f.write(f"ea_years_mean\t{statistics.mean(ea_values):.10g}\n")
        f.write(f"ea_years_median\t{statistics.median(ea_values):.10g}\n")
        f.write(f"yob_mean\t{mean_yob:.10g}\n")
        f.write(f"age_at_survey_mean\t{statistics.mean(ages):.10g}\n")
        f.write(f"age_at_survey_min\t{min(ages):.10g}\n")
        f.write(f"n_pcs\t{args.n_pcs}\n")
        f.write("covar_cols\tyob_c,sex_c,yob_c_sex_c_inter," + ",".join(pc_headers) + "\n")

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
    check("no GWAS IID is in sample-QC exclusion list", not (set(training_iids) & excluded_iids))
    check("all GWAS IIDs have sex covariate", set(training_iids) <= set(sex_map))
    check("all GWAS IIDs are at least min age at survey",
          all(ea_rows[iid]["age_at_survey"] >= args.min_age_at_survey for iid in training_iids))
    check("all GWAS IIDs have genotype FID", set(training_iids) <= set(fid_by_iid))
    check("all output FID/IID pairs match genotype fam",
          all(line.split()[0] == fid_by_iid[line.split()[1]]
              for line in open(training_path) if line.strip()))
    check("row count matches summary", len(training_iids) == len(gwas_iids))
    check("all EA values are mapped values", set(ea_values) <= set(EA_MAPPING.values()))
    check("mean yob_c is approximately zero",
          abs(statistics.mean(covar_data[iid][0] for iid in gwas_iids)) < 1e-8)
    check("sex_c values are {-0.5, 0.5}", {covar_data[iid][1] for iid in gwas_iids} == {-0.5, 0.5})
    check("all GWAS IIDs have requested PCs", set(training_iids) <= set(pc_data))
    if not nonlocal_passed[0]:
        sys.exit(1)

    print("\nAll verification checks passed.")


if __name__ == "__main__":
    main()
