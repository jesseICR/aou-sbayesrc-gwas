#!/usr/bin/env python3
"""Build AoU household-income GWAS phenotype/covariate files."""

import argparse
import csv
import os
import statistics
import sys
from collections import Counter


INCOME_MAPPING = {
    1585376: 5.0,     # <10k, midpoint of 0-10k
    1585377: 17.5,    # 10k-25k
    1585378: 30.0,    # 25k-35k
    1585379: 42.5,    # 35k-50k
    1585380: 62.5,    # 50k-75k
    1585381: 87.5,    # 75k-100k
    1585382: 125.0,   # 100k-150k
    1585383: 175.0,   # 150k-200k
    1585384: 250.0,   # >200k, top-coded as implicit 200k-300k midpoint
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


def load_income_rows(path):
    rows = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"IID", "income_k", "age_at_survey", "answer_concept_id", "answer", "n_income_records"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["IID"].strip()
            if not iid:
                continue
            answer_concept_id = int(row["answer_concept_id"])
            income_k = float(row["income_k"])
            if answer_concept_id not in INCOME_MAPPING:
                raise ValueError(f"Unexpected income answer_concept_id={answer_concept_id} for IID {iid}")
            if abs(income_k - INCOME_MAPPING[answer_concept_id]) > 1e-8:
                raise ValueError(f"Income mapping mismatch for answer_concept_id={answer_concept_id}")
            rows[iid] = {
                "income_k": income_k,
                "age_at_survey": float(row["age_at_survey"]),
                "answer_concept_id": answer_concept_id,
                "answer": row["answer"],
                "n_income_records": int(row["n_income_records"]),
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


def write_answer_counts(path, income_rows, gwas_iids):
    all_counts = Counter(row["answer_concept_id"] for row in income_rows.values())
    final_counts = Counter(income_rows[iid]["answer_concept_id"] for iid in gwas_iids)
    answer_names = {}
    for row in income_rows.values():
        answer_names[row["answer_concept_id"]] = row["answer"]

    with open(path, "w") as f:
        f.write("answer_concept_id\tanswer\tincome_k\tquery_count\tfinal_gwas_count\n")
        for concept_id in sorted(INCOME_MAPPING, key=lambda x: INCOME_MAPPING[x]):
            f.write(
                f"{concept_id}\t{answer_names.get(concept_id, '')}\t"
                f"{INCOME_MAPPING[concept_id]:.10g}\t{all_counts[concept_id]}\t{final_counts[concept_id]}\n"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--income-query", required=True)
    parser.add_argument("--europeans", required=True)
    parser.add_argument("--sex-covar", required=True)
    parser.add_argument("--exclude-iids", required=True)
    parser.add_argument("--fam", required=True)
    parser.add_argument("--sscore", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-pcs", type=int, default=10)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    log_lines = []
    log(log_lines, "=== AoU household-income GWAS setup ===")
    log(log_lines, f"income_query: {args.income_query}")
    log(log_lines, f"exclude_iids: {args.exclude_iids}")
    log(log_lines, f"n_pcs: {args.n_pcs}")

    europeans = read_keep_iids(args.europeans)
    sex_map = load_sex_covar(args.sex_covar)
    excluded_iids = read_keep_iids(args.exclude_iids)
    fid_by_iid = load_fam_fids(args.fam)
    income_rows = load_income_rows(args.income_query)
    pc_data, pc_headers = load_projected_pcs(args.sscore, args.n_pcs)

    log(log_lines, "")
    log(log_lines, "=== Input counts ===")
    log(log_lines, f"classified Europeans: {len(europeans)}")
    log(log_lines, f"sex covariate rows: {len(sex_map)}")
    log(log_lines, f"sample-QC exclusion IIDs: {len(excluded_iids)}")
    log(log_lines, f"fam rows: {len(fid_by_iid)}")
    log(log_lines, f"income query rows: {len(income_rows)}")
    log(log_lines, f"projected PC rows: {len(pc_data)}")

    duplicate_record_iids = {iid for iid, row in income_rows.items() if row["n_income_records"] > 1}
    excluded_europeans = europeans & excluded_iids
    candidates = set(europeans) - excluded_iids
    after_sample_qc = set(candidates)
    missing_income = candidates - set(income_rows)
    candidates &= set(income_rows)
    missing_sex = candidates - set(sex_map)
    candidates &= set(sex_map)
    missing_fam = candidates - set(fid_by_iid)
    candidates &= set(fid_by_iid)
    missing_pcs = candidates - set(pc_data)
    candidates &= set(pc_data)

    gwas_iids = sort_iids(candidates)
    if not gwas_iids:
        raise RuntimeError("No income GWAS samples remain after phenotype/sex/PC filters")

    ages = [income_rows[iid]["age_at_survey"] for iid in gwas_iids]
    mean_age = statistics.mean(ages)
    covar_data = {}
    for iid in gwas_iids:
        age_c = income_rows[iid]["age_at_survey"] - mean_age
        age_c_sq = age_c * age_c
        sex_c = sex_map[iid] - 0.5
        covar_data[iid] = (age_c, age_c_sq, sex_c, age_c * sex_c)

    males = sum(1 for iid in gwas_iids if sex_map[iid] == 1)
    females = sum(1 for iid in gwas_iids if sex_map[iid] == 0)
    income_values = [income_rows[iid]["income_k"] for iid in gwas_iids]

    log(log_lines, "")
    log(log_lines, "=== Filtering counts ===")
    log(log_lines, f"Europeans removed by sample-QC exclusion: {len(excluded_europeans)}")
    log(log_lines, f"Europeans after sample-QC exclusion: {len(after_sample_qc)}")
    log(log_lines, f"Europeans missing codeable income phenotype: {len(missing_income)}")
    log(log_lines, f"Europeans with income phenotype but missing sex covariate: {len(missing_sex)}")
    log(log_lines, f"Europeans with income+sex but missing fam row: {len(missing_fam)}")
    log(log_lines, f"Europeans with income+sex+fam but missing projected PCs: {len(missing_pcs)}")
    log(log_lines, f"Income query IIDs with multiple codeable records: {len(duplicate_record_iids)}")
    log(log_lines, f"Final GWAS samples: {len(gwas_iids)}")
    log(log_lines, f"Final female/male: {females}/{males}")
    log(log_lines, f"Mean age at income survey: {mean_age:.6f}")
    log(log_lines, f"Mean income_k: {statistics.mean(income_values):.6f}")
    log(log_lines, f"Median income_k: {statistics.median(income_values):.6f}")

    training_path = os.path.join(args.out_dir, "training_iids.txt")
    phen_path = os.path.join(args.out_dir, "phen.txt")
    base_covar_path = os.path.join(args.out_dir, "base_covar.txt")
    covar_path = os.path.join(args.out_dir, "covar.txt")
    answer_counts_path = os.path.join(args.out_dir, "income_answer_counts.tsv")
    log_path = os.path.join(args.out_dir, "income_gwas_log.txt")
    summary_path = os.path.join(args.out_dir, "income_gwas.summary.tsv")

    with open(training_path, "w") as f:
        for iid in gwas_iids:
            f.write(f"{fid_by_iid[iid]} {iid}\n")

    with open(phen_path, "w") as f:
        f.write("FID\tIID\tincome_k\n")
        for iid in gwas_iids:
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{income_rows[iid]['income_k']:.10g}\n")

    with open(base_covar_path, "w") as f:
        f.write("FID\tIID\tage_c\tage_c_sq\tsex_c\tage_c_sex_c_inter\n")
        for iid in gwas_iids:
            age_c, age_c_sq, sex_c, inter = covar_data[iid]
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{age_c:.12g}\t{age_c_sq:.12g}\t{sex_c:.1f}\t{inter:.12g}\n")

    with open(covar_path, "w") as f:
        f.write("FID\tIID\tage_c\tage_c_sq\tsex_c\tage_c_sex_c_inter\t" + "\t".join(pc_headers) + "\n")
        for iid in gwas_iids:
            age_c, age_c_sq, sex_c, inter = covar_data[iid]
            pcs = "\t".join(pc_data[iid])
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{age_c:.12g}\t{age_c_sq:.12g}\t{sex_c:.1f}\t{inter:.12g}\t{pcs}\n")

    write_answer_counts(answer_counts_path, income_rows, gwas_iids)

    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"classified_europeans\t{len(europeans)}\n")
        f.write(f"sex_covar_rows\t{len(sex_map)}\n")
        f.write(f"sample_qc_exclusion_iids\t{len(excluded_iids)}\n")
        f.write(f"classified_europeans_removed_by_sample_qc\t{len(excluded_europeans)}\n")
        f.write(f"classified_europeans_after_sample_qc\t{len(after_sample_qc)}\n")
        f.write(f"fam_rows\t{len(fid_by_iid)}\n")
        f.write(f"income_query_rows\t{len(income_rows)}\n")
        f.write(f"projected_pc_rows\t{len(pc_data)}\n")
        f.write(f"europeans_missing_codeable_income\t{len(missing_income)}\n")
        f.write(f"income_candidates_missing_sex_covar\t{len(missing_sex)}\n")
        f.write(f"income_sex_candidates_missing_fam\t{len(missing_fam)}\n")
        f.write(f"income_sex_fam_candidates_missing_pcs\t{len(missing_pcs)}\n")
        f.write(f"income_query_iids_with_multiple_codeable_records\t{len(duplicate_record_iids)}\n")
        f.write(f"gwas_samples\t{len(gwas_iids)}\n")
        f.write(f"gwas_female\t{females}\n")
        f.write(f"gwas_male\t{males}\n")
        f.write(f"income_k_mean\t{statistics.mean(income_values):.10g}\n")
        f.write(f"income_k_median\t{statistics.median(income_values):.10g}\n")
        f.write(f"age_at_survey_mean\t{mean_age:.10g}\n")
        f.write(f"n_pcs\t{args.n_pcs}\n")
        f.write("covar_cols\tage_c,age_c_sq,sex_c,age_c_sex_c_inter," + ",".join(pc_headers) + "\n")

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
    check("all GWAS IIDs have genotype FID", set(training_iids) <= set(fid_by_iid))
    check("all output FID/IID pairs match genotype fam",
          all(line.split()[0] == fid_by_iid[line.split()[1]]
              for line in open(training_path) if line.strip()))
    check("row count matches summary", len(training_iids) == len(gwas_iids))
    check("all income values are mapped values", set(income_values) <= set(INCOME_MAPPING.values()))
    check("mean age_c is approximately zero",
          abs(statistics.mean(covar_data[iid][0] for iid in gwas_iids)) < 1e-8)
    check("age_c_sq values are nonnegative", all(covar_data[iid][1] >= 0 for iid in gwas_iids))
    check("sex_c values are {-0.5, 0.5}", {covar_data[iid][2] for iid in gwas_iids} == {-0.5, 0.5})
    check("all GWAS IIDs have requested PCs", set(training_iids) <= set(pc_data))
    if not nonlocal_passed[0]:
        sys.exit(1)

    print("\nAll verification checks passed.")


if __name__ == "__main__":
    main()
