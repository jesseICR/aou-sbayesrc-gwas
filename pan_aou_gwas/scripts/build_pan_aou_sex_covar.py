#!/usr/bin/env python3
"""Build the pan-AoU GWAS sex covariate.

The core genotype pipeline's sex covariate is deliberately strict: it keeps only
samples with binary sex-at-birth metadata concordant with canonical WGS sex
ploidy (XX/XY).  The pan-AoU atlas uses that file as the base covariate, then
adds a small set of pre-specified rows from sex_ploidy_qc.tsv so phenotype-wide
GWASes do not lose otherwise usable unrelated-European samples:

* assigned sex at birth Male with DRAGEN X0/XO is coded male;
* skipped/prefer-not-to-answer sex at birth with DRAGEN XX is coded female;
* skipped/prefer-not-to-answer sex at birth with DRAGEN XY is coded male.

Other nonbinary/missing/discrepant rows remain excluded.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


WITHHELD_SEX_LABELS = {
    "pmi skip",
    "i prefer not to answer",
}


def iid_sort_key(iid: str) -> tuple[int, int | str]:
    return (0, int(iid)) if iid.isdigit() else (1, iid)


def normalize_label(value: str) -> str:
    text = re.sub(r"\s+", " ", str(value or "").replace("_", " ")).strip().lower()
    text = text.replace(":", " ")
    text = re.sub(r"\s+", " ", text)
    return text


def load_strict_sex_covar(path: Path) -> dict[str, tuple[str, int]]:
    out: dict[str, tuple[str, int]] = {}
    with path.open(newline="") as f:
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
                raise ValueError(f"{path}: nonbinary sex_01 for IID {iid}: {sex_01}")
            out[iid] = (row["FID"].strip() or "0", sex_01)
    return out


def imputed_sex_from_qc_row(row: dict[str, str]) -> tuple[int | None, str]:
    assigned_sex_01 = (row.get("sex_01") or "").strip()
    ploidy = (row.get("dragen_sex_ploidy") or "").strip().upper()
    concept_name = (row.get("sex_at_birth_concept_name") or "").strip()
    source_value = (row.get("sex_at_birth_source_value") or "").strip()
    labels = {normalize_label(concept_name), normalize_label(source_value)}

    if assigned_sex_01 == "1" and ploidy in {"X0", "XO"}:
        return 1, "assigned_male_dragen_x0_xo"

    if assigned_sex_01 == "" and labels & WITHHELD_SEX_LABELS:
        if ploidy == "XX":
            return 0, "withheld_sex_at_birth_dragen_xx"
        if ploidy == "XY":
            return 1, "withheld_sex_at_birth_dragen_xy"

    return None, ""


def build_pan_aou_sex_covar(
    strict_sex_covar: Path,
    sex_ploidy_qc: Path,
    out: Path,
    summary: Path,
    audit_out: Path,
) -> None:
    sex = load_strict_sex_covar(strict_sex_covar)
    strict_iids = set(sex)

    summary_counts = Counter()
    summary_counts["strict_input_rows"] = len(sex)
    added_rows: list[dict[str, str]] = []

    with sex_ploidy_qc.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {
            "FID",
            "IID",
            "sex_01",
            "sex_at_birth_concept_name",
            "sex_at_birth_source_value",
            "dragen_sex_ploidy",
            "sex_ploidy_status",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{sex_ploidy_qc} missing columns: {sorted(missing)}")

        for row in reader:
            summary_counts["sex_ploidy_qc_rows"] += 1
            iid = row["IID"].strip()
            if not iid or iid in strict_iids:
                continue
            imputed_sex_01, rule = imputed_sex_from_qc_row(row)
            if imputed_sex_01 is None:
                continue
            fid = row["FID"].strip() or "0"
            sex[iid] = (fid, imputed_sex_01)
            summary_counts[f"added_{rule}"] += 1
            added_rows.append({
                "FID": fid,
                "IID": iid,
                "sex_01": str(imputed_sex_01),
                "imputation_rule": rule,
                "sex_at_birth_concept_name": row.get("sex_at_birth_concept_name", ""),
                "sex_at_birth_source_value": row.get("sex_at_birth_source_value", ""),
                "dragen_sex_ploidy": row.get("dragen_sex_ploidy", ""),
                "sex_ploidy_status": row.get("sex_ploidy_status", ""),
            })

    summary_counts["added_total"] = len(added_rows)
    summary_counts["pan_aou_sex_covar_rows"] = len(sex)
    summary_counts["pan_aou_female"] = sum(1 for _, sex_01 in sex.values() if sex_01 == 0)
    summary_counts["pan_aou_male"] = sum(1 for _, sex_01 in sex.values() if sex_01 == 1)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("FID\tIID\tsex_01\n")
        for iid in sorted(sex, key=iid_sort_key):
            fid, sex_01 = sex[iid]
            f.write(f"{fid}\t{iid}\t{sex_01}\n")

    audit_out.parent.mkdir(parents=True, exist_ok=True)
    audit_fields = [
        "FID",
        "IID",
        "sex_01",
        "imputation_rule",
        "sex_at_birth_concept_name",
        "sex_at_birth_source_value",
        "dragen_sex_ploidy",
        "sex_ploidy_status",
    ]
    with audit_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=audit_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(added_rows, key=lambda r: iid_sort_key(r["IID"])))

    summary.parent.mkdir(parents=True, exist_ok=True)
    ordered = [
        "strict_input_rows",
        "sex_ploidy_qc_rows",
        "added_assigned_male_dragen_x0_xo",
        "added_withheld_sex_at_birth_dragen_xx",
        "added_withheld_sex_at_birth_dragen_xy",
        "added_total",
        "pan_aou_sex_covar_rows",
        "pan_aou_female",
        "pan_aou_male",
    ]
    with summary.open("w") as f:
        f.write("metric\tvalue\n")
        for metric in ordered:
            f.write(f"{metric}\t{summary_counts.get(metric, 0)}\n")

    print(
        "  pan_aou_sex_covar="
        f"{summary_counts['pan_aou_sex_covar_rows']} "
        f"(strict={summary_counts['strict_input_rows']}, "
        f"added={summary_counts['added_total']})"
    )
    print(f"  summary={summary}")
    print(f"  audit={audit_out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict-sex-covar", type=Path, required=True)
    parser.add_argument("--sex-ploidy-qc", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--audit-out", type=Path, required=True)
    args = parser.parse_args()

    build_pan_aou_sex_covar(
        strict_sex_covar=args.strict_sex_covar,
        sex_ploidy_qc=args.sex_ploidy_qc,
        out=args.out,
        summary=args.summary,
        audit_out=args.audit_out,
    )


if __name__ == "__main__":
    main()
