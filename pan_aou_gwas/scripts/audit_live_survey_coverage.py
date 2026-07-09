#!/usr/bin/env python3
"""Audit live AoU survey question coverage against pan-AoU GWAS manifests.

The invariant checked here is intentionally simple:

  every live survey question observed in the GWAS sample must be represented in
  phenotype_manifest.tsv, represented in skipped_phenotypes.tsv, or explicitly
  excluded by metadata. Anything else is a coverage gap.

This audits question-level coverage. The phenotype builder itself remains the
answer-level authority for one-vs-rest binary phenotypes and writes QC failures
to skipped_phenotypes.tsv.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pan_aou_gwas as pan  # noqa: E402


EXPECTED_DISPOSITIONS = {
    "ordinal_and_binary",
    "binary_only",
    "nominal_binary",
    "flagged_review",
    "numeric",
}


def load_keep(path: Path | None) -> set[str] | None:
    if not path:
        return None
    keep = set()
    with open(path) as f:
        for line in f:
            parts = line.split()
            if parts:
                keep.add(parts[-1])
    return keep


def suppressed_n(n: int) -> str:
    if n == 0:
        return "0"
    if 1 <= n <= 20:
        return "<=20"
    return str(n)


def scan_live_questions(paths: list[Path], keep: set[str] | None):
    seen = {}
    sample_sets = defaultdict(set)
    row_counts = Counter()
    answers_by_qid = defaultdict(set)
    for path in paths:
        if not path or not path.exists() or path.stat().st_size == 0:
            continue
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                pid = (row.get("person_id") or row.get("IID") or "").strip()
                if keep is not None and pid not in keep:
                    continue
                qid = (row.get("question_concept_id") or "").strip()
                if not qid:
                    continue
                survey = (row.get("survey") or "").strip()
                question = (row.get("question") or "").strip()
                if qid not in seen:
                    seen[qid] = {"survey": survey, "question": question}
                elif question and not seen[qid].get("question"):
                    seen[qid]["question"] = question
                row_counts[qid] += 1
                answer = (row.get("answer") or "").strip()
                if answer and not pan.is_missing_answer(answer):
                    answers_by_qid[qid].add(answer)
                if pid:
                    sample_sets[qid].add(pid)
    sample_counts = {qid: len(samples) for qid, samples in sample_sets.items()}
    return seen, row_counts, sample_counts, answers_by_qid


def load_manifest_by_qid(path: Path):
    by_qid = defaultdict(list)
    if not path or not path.exists() or path.stat().st_size == 0:
        return by_qid
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            qid = (row.get("question_concept_id") or "").strip()
            if qid:
                by_qid[qid].append(row)
    return by_qid


def load_manifest_by_pheno_id(path: Path):
    by_pid = defaultdict(list)
    if not path or not path.exists() or path.stat().st_size == 0:
        return by_pid
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("pheno_id") or "").strip()
            if pid:
                by_pid[pid].append(row)
    return by_pid


def build_qman(args):
    qman, _qman_by_item, qman_rows = pan.load_question_manifest(args.question_manifest)
    qman_by_qid, _qid_by_item = pan.load_live_question_crosswalk(
        args.aou_question_concepts, qman_rows, qman
    )
    qman.update(qman_by_qid)
    ea_proxy_rows = pan.load_ea_proxy_feature_sources(args.ea_proxy_feature_manifest, qman)
    qman_rows.extend(ea_proxy_rows)
    for row in ea_proxy_rows:
        qid = (row.get("question_concept_id") or "").strip()
        if qid:
            qman.setdefault(qid, row)
        label = (row.get("field_label") or "").strip()
        if label:
            qman.setdefault(pan.norm_q(label), row)
    live_override_rows = pan.load_live_question_overrides(args.aou_question_concepts, qman)
    qman_rows.extend(live_override_rows)
    for row in live_override_rows:
        qid = (row.get("question_concept_id") or "").strip()
        if qid:
            qman.setdefault(qid, row)
        label = (row.get("field_label") or "").strip()
        if label:
            qman.setdefault(pan.norm_q(label), row)
    return qman, qman_rows, qman_by_qid, ea_proxy_rows, live_override_rows


def metadata_for_live_question(qid: str, question: str, qman: dict):
    return qman.get(qid) or qman.get(pan.norm_q(question))


def classify(qid: str, question: str, qman: dict, passed: list[dict], skipped: list[dict]):
    if passed and skipped:
        return "covered_passed_and_qc_skipped"
    if passed:
        return "covered_passed"
    if skipped:
        return "covered_qc_skipped"

    man = metadata_for_live_question(qid, question, qman)
    if man is None:
        return "missing_metadata_no_encoding_rule"
    disp = (man.get("disposition") or "").strip()
    if disp.startswith("excluded"):
        return "explicitly_excluded"
    if disp in EXPECTED_DISPOSITIONS:
        return "missing_expected_no_manifest_or_skip"
    return "metadata_non_gwas_disposition"


def answer_status(pheno_id: str, passed_by_pid: dict, skipped_by_pid: dict):
    passed = passed_by_pid.get(pheno_id, [])
    skipped = skipped_by_pid.get(pheno_id, [])
    if passed and skipped:
        return "covered_passed_and_qc_skipped"
    if passed:
        return "covered_passed"
    if skipped:
        return "covered_qc_skipped"
    return "missing_expected_no_manifest_or_skip"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--survey-csv", type=Path, required=True)
    ap.add_argument("--bhp-csv", type=Path, default=None)
    ap.add_argument("--keep", type=Path, default=None)
    ap.add_argument("--question-manifest", type=Path, required=True)
    ap.add_argument("--aou-question-concepts", type=Path, required=True)
    ap.add_argument("--ea-proxy-feature-manifest", type=Path, default=None)
    ap.add_argument("--phenotype-manifest", type=Path, required=True)
    ap.add_argument("--skipped-phenotypes", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    keep = load_keep(args.keep)
    qman, qman_rows, qman_by_qid, ea_proxy_rows, live_override_rows = build_qman(args)
    live, row_counts, sample_counts, answers_by_qid = scan_live_questions(
        [args.survey_csv, args.bhp_csv], keep
    )
    passed_by_qid = load_manifest_by_qid(args.phenotype_manifest)
    skipped_by_qid = load_manifest_by_qid(args.skipped_phenotypes)
    passed_by_pid = load_manifest_by_pheno_id(args.phenotype_manifest)
    skipped_by_pid = load_manifest_by_pheno_id(args.skipped_phenotypes)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    answer_rows = []
    for qid in sorted(live, key=lambda x: (0, int(x)) if x.isdigit() else (1, x)):
        man = metadata_for_live_question(qid, live[qid].get("question", ""), qman)
        passed = passed_by_qid.get(qid, [])
        skipped = skipped_by_qid.get(qid, [])
        status = classify(qid, live[qid].get("question", ""), qman, passed, skipped)
        skip_reasons = sorted({r.get("skip_reason", "") for r in skipped if r.get("skip_reason")})
        rows.append({
            "question_concept_id": qid,
            "status": status,
            "survey": live[qid].get("survey", ""),
            "question": live[qid].get("question", ""),
            "metadata_disposition": (man or {}).get("disposition", ""),
            "metadata_item_concept": (man or {}).get("item_concept", ""),
            "metadata_ordinal_rule": (man or {}).get("ordinal_rule", ""),
            "n_passed_phenotypes": len(passed),
            "n_skipped_phenotypes": len(skipped),
            "skip_reasons": ",".join(skip_reasons),
            "n_rows_suppressed": suppressed_n(row_counts[qid]),
            "n_samples_suppressed": suppressed_n(sample_counts.get(qid, 0)),
        })
        if man is None:
            continue
        disp = (man.get("disposition") or "").strip()
        if disp.startswith("excluded"):
            continue
        item_id = pan.manifest_item_id(man, qid)
        if disp in ("ordinal_and_binary", "binary_only", "nominal_binary", "flagged_review"):
            for answer in sorted(answers_by_qid.get(qid, set()), key=pan.answer_slug):
                pheno_id = f"bin_{item_id}__{pan.answer_slug(answer)}"
                status = answer_status(pheno_id, passed_by_pid, skipped_by_pid)
                skip_rows = skipped_by_pid.get(pheno_id, [])
                answer_rows.append({
                    "question_concept_id": qid,
                    "answer": answer,
                    "expected_pheno_id": pheno_id,
                    "status": status,
                    "survey": live[qid].get("survey", ""),
                    "question": live[qid].get("question", ""),
                    "metadata_disposition": disp,
                    "metadata_item_concept": man.get("item_concept", ""),
                    "skip_reasons": ",".join(sorted({r.get("skip_reason", "") for r in skip_rows if r.get("skip_reason")})),
                })
        if disp == "ordinal_and_binary":
            pheno_id = f"ord_{item_id}"
            status = answer_status(pheno_id, passed_by_pid, skipped_by_pid)
            skip_rows = skipped_by_pid.get(pheno_id, [])
            answer_rows.append({
                "question_concept_id": qid,
                "answer": "",
                "expected_pheno_id": pheno_id,
                "status": status,
                "survey": live[qid].get("survey", ""),
                "question": live[qid].get("question", ""),
                "metadata_disposition": disp,
                "metadata_item_concept": man.get("item_concept", ""),
                "skip_reasons": ",".join(sorted({r.get("skip_reason", "") for r in skip_rows if r.get("skip_reason")})),
            })
        if disp == "numeric":
            pheno_id = f"num_{item_id}"
            status = answer_status(pheno_id, passed_by_pid, skipped_by_pid)
            skip_rows = skipped_by_pid.get(pheno_id, [])
            answer_rows.append({
                "question_concept_id": qid,
                "answer": "",
                "expected_pheno_id": pheno_id,
                "status": status,
                "survey": live[qid].get("survey", ""),
                "question": live[qid].get("question", ""),
                "metadata_disposition": disp,
                "metadata_item_concept": man.get("item_concept", ""),
                "skip_reasons": ",".join(sorted({r.get("skip_reason", "") for r in skip_rows if r.get("skip_reason")})),
            })

    fieldnames = [
        "question_concept_id",
        "status",
        "survey",
        "question",
        "metadata_disposition",
        "metadata_item_concept",
        "metadata_ordinal_rule",
        "n_passed_phenotypes",
        "n_skipped_phenotypes",
        "skip_reasons",
        "n_rows_suppressed",
        "n_samples_suppressed",
    ]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    answer_out = args.out.with_name(args.out.stem.replace(".question", "") + ".answer_coverage.tsv")
    answer_fields = [
        "question_concept_id",
        "answer",
        "expected_pheno_id",
        "status",
        "survey",
        "question",
        "metadata_disposition",
        "metadata_item_concept",
        "skip_reasons",
    ]
    with open(answer_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=answer_fields, delimiter="\t")
        w.writeheader()
        w.writerows(answer_rows)

    status_counts = Counter(r["status"] for r in rows)
    answer_status_counts = Counter(r["status"] for r in answer_rows)
    print(f"live_questions={len(rows)}")
    print(f"expected_answer_level_phenotypes={len(answer_rows)}")
    print(f"question_manifest_rows={len(qman_rows)}")
    print(f"live_qid_links={len(qman_by_qid)}")
    print(f"ea_proxy_supplemental_rows={len(ea_proxy_rows)}")
    print(f"live_qid_override_rows={len(live_override_rows)}")
    print("question_status_counts")
    for status, n in status_counts.most_common():
        print(f"{status}\t{n}")
    print("answer_status_counts")
    for status, n in answer_status_counts.most_common():
        print(f"{status}\t{n}")
    print(f"wrote={args.out}")
    print(f"wrote={answer_out}")


if __name__ == "__main__":
    main()
