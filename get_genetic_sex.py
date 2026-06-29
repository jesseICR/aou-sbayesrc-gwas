#!/usr/bin/env python3
"""Build sex covariate and sex-at-birth/WGS-ploidy QC outputs for AoU."""

import argparse
import csv
import os
from collections import Counter


FEMALE_CONCEPT_ID = "45878463"
MALE_CONCEPT_ID = "45880669"
FEMALE_PLOIDY = "XX"
MALE_PLOIDY = "XY"


def sort_iids(iids):
    return sorted(iids, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))


def load_fam_samples(path):
    iids = []
    fid_by_iid = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                iids.append(parts[1])
                fid_by_iid[parts[1]] = parts[0]
    return iids, fid_by_iid


def load_sex_query(path):
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {
            "IID",
            "sex_at_birth_concept_id",
            "sex_at_birth_concept_name",
            "sex_at_birth_source_concept_id",
            "sex_at_birth_source_concept_name",
            "sex_at_birth_source_value",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["IID"].strip()
            if iid:
                out[iid] = row
    return out


def load_genomic_metrics(path):
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"research_id", "dragen_sex_ploidy"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["research_id"].strip()
            if iid:
                out[iid] = row
    return out


def sex_01_from_concept(concept_id):
    if concept_id == FEMALE_CONCEPT_ID:
        return 0
    if concept_id == MALE_CONCEPT_ID:
        return 1
    return None


def expected_ploidy(sex_01):
    if sex_01 == 0:
        return FEMALE_PLOIDY
    if sex_01 == 1:
        return MALE_PLOIDY
    return None


def status_for(sex_01, ploidy):
    if sex_01 is None:
        return "missing_or_nonbinary_sex_at_birth"
    if not ploidy:
        return "missing_wgs_sex_ploidy"
    exp = expected_ploidy(sex_01)
    if ploidy == exp:
        return "concordant_binary"
    if ploidy in {FEMALE_PLOIDY, MALE_PLOIDY}:
        return "discordant_binary"
    return "noncanonical_wgs_sex_ploidy"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fam", required=True)
    parser.add_argument("--sex-query", required=True)
    parser.add_argument("--genomic-metrics", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--require-ploidy-concordance", type=int, choices=[0, 1], default=1)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    sample_iids, fid_by_iid = load_fam_samples(args.fam)
    sample_set = set(sample_iids)
    sex_rows = load_sex_query(args.sex_query)
    genomic_rows = load_genomic_metrics(args.genomic_metrics)

    summary = Counter()
    crosstab = Counter()
    qc_records = []
    sex_covar = {}

    for iid in sample_iids:
        sex_row = sex_rows.get(iid)
        genomic_row = genomic_rows.get(iid)

        if sex_row is None:
            concept_id = ""
            concept_name = ""
            source_value = ""
            sex_01 = None
            summary["missing_sex_at_birth_row"] += 1
        else:
            concept_id = sex_row.get("sex_at_birth_concept_id", "").strip()
            concept_name = sex_row.get("sex_at_birth_concept_name", "").strip()
            source_value = sex_row.get("sex_at_birth_source_value", "").strip()
            sex_01 = sex_01_from_concept(concept_id)
            if sex_01 is None:
                pass
            else:
                summary["binary_sex_at_birth"] += 1

        if genomic_row is None:
            metrics_sex_at_birth = ""
            ploidy = ""
            summary["missing_genomic_metrics_row"] += 1
        else:
            metrics_sex_at_birth = genomic_row.get("sex_at_birth", "").strip()
            ploidy = genomic_row.get("dragen_sex_ploidy", "").strip()

        status = status_for(sex_01, ploidy)
        summary[status] += 1

        sex_label = "Female" if sex_01 == 0 else "Male" if sex_01 == 1 else (concept_name or source_value or "MISSING_OR_NONBINARY")
        crosstab[(sex_label, "" if sex_01 is None else str(sex_01), ploidy or "MISSING", status)] += 1

        keep = (sex_01 in (0, 1)) and (
            (not args.require_ploidy_concordance) or status == "concordant_binary"
        )
        if keep:
            sex_covar[iid] = sex_01

        qc_records.append({
            "FID": fid_by_iid[iid],
            "IID": iid,
            "sex_01": "" if sex_01 is None else str(sex_01),
            "sex_at_birth_concept_id": concept_id,
            "sex_at_birth_concept_name": concept_name,
            "sex_at_birth_source_value": source_value,
            "genomic_metrics_sex_at_birth": metrics_sex_at_birth,
            "dragen_sex_ploidy": ploidy,
            "expected_dragen_sex_ploidy": expected_ploidy(sex_01) or "",
            "sex_ploidy_status": status,
            "kept_in_sex_covar": "1" if keep else "0",
        })

    summary["sample_universe"] = len(sample_iids)
    summary["sex_query_rows_total"] = len(sex_rows)
    summary["genomic_metrics_rows_total"] = len(genomic_rows)
    summary["confident_sex_samples"] = len(sex_covar)
    summary["confident_female"] = sum(1 for v in sex_covar.values() if v == 0)
    summary["confident_male"] = sum(1 for v in sex_covar.values() if v == 1)
    summary["require_ploidy_concordance"] = args.require_ploidy_concordance

    universe = len(sample_iids)
    confident_pct = 100.0 * len(sex_covar) / universe if universe else 0.0
    summary["confident_sex_percent"] = f"{confident_pct:.6f}"

    sex_covar_path = os.path.join(args.out_dir, "sex_covar.txt")
    with open(sex_covar_path, "w") as f:
        f.write("FID\tIID\tsex_01\n")
        for iid in sort_iids(sex_covar):
            f.write(f"{fid_by_iid[iid]}\t{iid}\t{sex_covar[iid]}\n")

    qc_path = os.path.join(args.out_dir, "sex_ploidy_qc.tsv")
    qc_fields = [
        "FID", "IID", "sex_01", "sex_at_birth_concept_id",
        "sex_at_birth_concept_name", "sex_at_birth_source_value",
        "genomic_metrics_sex_at_birth", "dragen_sex_ploidy",
        "expected_dragen_sex_ploidy", "sex_ploidy_status",
        "kept_in_sex_covar",
    ]
    with open(qc_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=qc_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_records)

    crosstab_path = os.path.join(args.out_dir, "sex_ploidy_crosstab.tsv")
    with open(crosstab_path, "w") as f:
        f.write("sex_at_birth_label\tsex_01\tdragen_sex_ploidy\tsex_ploidy_status\tn\n")
        for (label, sex_01, ploidy, status), n in sorted(crosstab.items(), key=lambda x: (x[0][3], x[0][0], x[0][2])):
            f.write(f"{label}\t{sex_01}\t{ploidy}\t{status}\t{n}\n")

    summary_path = os.path.join(args.out_dir, "genetic_sex_summary.tsv")
    ordered_metrics = [
        "sample_universe",
        "sex_query_rows_total",
        "genomic_metrics_rows_total",
        "binary_sex_at_birth",
        "missing_sex_at_birth_row",
        "missing_or_nonbinary_sex_at_birth",
        "missing_genomic_metrics_row",
        "missing_wgs_sex_ploidy",
        "concordant_binary",
        "discordant_binary",
        "noncanonical_wgs_sex_ploidy",
        "confident_sex_samples",
        "confident_sex_percent",
        "confident_female",
        "confident_male",
        "require_ploidy_concordance",
    ]
    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        for metric in ordered_metrics:
            f.write(f"{metric}\t{summary.get(metric, 0)}\n")

    log_path = os.path.join(args.out_dir, "genetic_sex_log.txt")
    with open(log_path, "w") as f:
        f.write("=== AoU sex-at-birth / WGS sex-ploidy QC ===\n")
        for metric in ordered_metrics:
            f.write(f"{metric}: {summary.get(metric, 0)}\n")
        f.write("\nCrosstab written to sex_ploidy_crosstab.tsv\n")
        f.write("Per-sample QC table written to sex_ploidy_qc.tsv\n")
        f.write("Covariate file written to sex_covar.txt\n")

    print("=== AoU sex-at-birth / WGS sex-ploidy QC ===")
    for metric in ordered_metrics:
        print(f"{metric}: {summary.get(metric, 0)}")
    print("\n=== Verification checks ===")
    all_passed = True

    def check(name, condition, detail=""):
        nonlocal all_passed
        all_passed = all_passed and bool(condition)
        print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))

    check("sample universe has unique IIDs", len(sample_iids) == len(sample_set))
    check("sex_covar rows match confident_sex_samples", len(sex_covar) == int(summary["confident_sex_samples"]))
    check("sex_covar values are binary", set(sex_covar.values()) <= {0, 1})
    if args.require_ploidy_concordance:
        check("all sex_covar samples are concordant",
              all(r["kept_in_sex_covar"] == "0" or r["sex_ploidy_status"] == "concordant_binary"
                  for r in qc_records))

    if not all_passed:
        raise SystemExit(1)
    print("\nAll verification checks passed.")


if __name__ == "__main__":
    main()
