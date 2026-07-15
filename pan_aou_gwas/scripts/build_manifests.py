#!/usr/bin/env python3
"""Turn the survey item inventory into the four phenotype manifests.

Reads `metadata/survey_item_inventory.tsv` (from parse_codebooks.py) and writes:

  metadata/survey_question_manifest.tsv  -- one row per question, its disposition
  metadata/ordinal_answer_templates.tsv  -- the shared ordinal template library
  metadata/ordinal_mapping_manifest.tsv  -- one row per (ordinal question, answer)
  metadata/flagged_questions.tsv         -- sensitive and uncertain-ordinal items

Dispositions
  ordinal_and_binary  single-select mapped to an ordinal scale (+ per-answer binary)
  binary_only         single/multi-select, no defensible ordinal scale
  nominal_binary      single-select that is explicitly nominal (binary one-vs-rest)
  numeric             free-numeric survey entry (linear GWAS on the number)
  excluded_*          not a GWAS phenotype (with reason)
  flagged_review      single-select we could not confidently classify -> review
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import ordinal_rules as R


INCLUDED_SURVEYS = {
    "The Basics",
    "Lifestyle",
    "Overall Health",
    "Healthcare Access and Utilization",
    "Social Determinants of Health",
    "COVID-19 Participant Experience (COPE)",
    "Minute Survey on COVID-19 Vaccines",
    "Life Functioning",
    "Emotional Health History and Well-Being",
    "Behavioral Health and Personality",
}

# Legacy / family-history sheets: not run as broad phenotypes. Only the PFHH
# self-history allowlist (maintained in the specsheet + pipeline, 33 binary
# phenotypes) is used from these.
FAMILY_HISTORY_SURVEYS = {
    "Personal Medical History",
    "Family Health History",
    "Personal and Family Health History",
}

# Item/label patterns that are excluded even inside included surveys.
EXCLUDE_ITEM_RE = re.compile(
    r"^(race_|whatraceethnicity|aian|asian|black|hispanic|mena|nhpi|white_|"
    r"noneofthese|socialsecurity|thebasics_countryborn)",
    re.I,
)
# Note: address/phone/email/ZIP entry fields are free-text and already excluded
# as such, so they are intentionally NOT listed here (listing "address" wrongly
# caught "How many years have you lived at your current address?").
EXCLUDE_LABEL_RE = re.compile(
    r"\b(race|ethnic|hispanic|latino|tribe|social security number|"
    r"please specify|which categories describe you)\b",
    re.I,
)


def parse_options(optstr: str) -> list[tuple[str, str]]:
    out = []
    for tok in optstr.split(" | "):
        if "=" in tok:
            code, label = tok.split("=", 1)
            out.append((code.strip(), label.strip()))
    return out


def core_labels(options: list[tuple[str, str]]) -> list[str]:
    return [lab for _, lab in options if not R.is_missing(lab)]


def sig_of(options) -> frozenset:
    return frozenset(R.norm(l) for l in core_labels(options))


def try_ordinal(item_concept: str, options: list[tuple[str, str]]):
    """Return (rule_name, source, confidence, {norm_label: value}) or None."""
    labels = core_labels(options)
    norm_labels = [R.norm(l) for l in labels]

    # 1. explicit per-item override
    ov = R.ITEM_OVERRIDES.get(item_concept)
    if ov is not None:
        mapping = {}
        ok = True
        for nl in norm_labels:
            if nl in ov["map"]:
                mapping[nl] = ov["map"][nl]
            elif nl in ov["local_missing"]:
                continue
            else:
                ok = False
                break
        if ok and mapping:
            return ov["rule"], "override", ov["confidence"], mapping

    # 2. shared templates
    for name, tpl in R.TEMPLATES.items():
        mapping = {}
        ok = True
        for nl in norm_labels:
            if nl in tpl["map"]:
                mapping[nl] = tpl["map"][nl]
            elif nl in tpl["local_missing"]:
                continue
            else:
                ok = False
                break
        if ok and mapping:
            n_levels = len(set(mapping.values()))
            signed = any(v < 0 for v in mapping.values())
            if n_levels >= (2 if signed else 3):
                return name, "template", tpl["confidence"], mapping
    return None


def is_nominal(item_concept: str, options) -> bool:
    if item_concept in R.NOMINAL_ITEMS:
        return True
    sig = sig_of(options)
    for hint in R.NOMINAL_SIGNATURE_HINTS:
        if hint <= sig:
            return True
    return False


def has_numeric_sibling(item: str, all_items: set) -> bool:
    """True if a dedicated numeric-entry field exists for this gated radio."""
    return any(item + suf in all_items for suf in ("number", "_number", "age", "_age"))


def disposition(row, all_items) -> dict:
    survey = row["survey"]
    item = row["item_concept"]
    label = row["field_label"]
    cls = row["phenotype_class"]
    options = parse_options(row["options"])
    n_core = len(core_labels(options))
    sens = R.sensitive_topics(label, survey)

    base = {
        "disposition": "",
        "ordinal_rule": "",
        "ordinal_source": "",
        "ordinal_confidence": "",
        "n_binary_phenos": 0,
        "n_ordinal_levels": 0,
        "sensitive_topics": ";".join(sens),
        "flag_reason": "",
        "notes": "",
        "_ordinal_map": {},
    }

    # Non-phenotype field types
    if cls in ("descriptive",):
        base["disposition"] = "excluded_descriptive"
        return base
    if cls == "date_text":
        base["disposition"] = "excluded_date"
        return base
    if cls == "free_text":
        base["disposition"] = "excluded_free_text"
        return base

    # Explicit exclusions inside included surveys (race/ethnicity/PII/specify)
    if EXCLUDE_ITEM_RE.match(item) or EXCLUDE_LABEL_RE.search(label):
        base["disposition"] = "excluded_race_pii_admin"
        return base

    # Family-history sheets: excluded except the PFHH self allowlist (separate).
    if survey in FAMILY_HISTORY_SURVEYS:
        base["disposition"] = "excluded_family_history"
        base["notes"] = "PFHH self-history allowlist handled separately (see specsheet 11.3)."
        return base

    if survey not in INCLUDED_SURVEYS:
        base["disposition"] = "excluded_other_survey"
        return base

    # Numeric free-entry
    if cls == "numeric_text":
        base["disposition"] = "numeric"
        base["notes"] = f"range [{row['validation_min']}, {row['validation_max']}]"
        return base

    # Multi-select checkbox: binary one-vs-rest per option, never ordinal.
    if cls == "multi_select":
        base["disposition"] = "binary_only"
        base["n_binary_phenos"] = n_core
        base["notes"] = "checkbox: one selected-vs-not GWAS per option."
        return base

    # Single-select (radio/dropdown/yesno) and slider
    if cls in ("single_select", "slider"):
        base["n_binary_phenos"] = n_core  # one-vs-rest per valid option
        # "Enter <the amount of time / number / age / response>" radios gate a
        # free numeric entry: the real datum is the number (value_as_number), so
        # treat the question as a continuous phenotype, not a degenerate radio.
        core = [lab for _, lab in options if not R.is_missing(lab)]
        enter_numeric = re.compile(
            r"\s*enter\s+(the\s+)?(number|amount|time|age|years?|months?|weeks?|days?|minutes?|hours?)\b",
            re.I,
        )
        if core and all(enter_numeric.match(lab) for lab in core):
            if has_numeric_sibling(item, all_items):
                # A dedicated _number/_age field already captures the value.
                base["disposition"] = "excluded_gated_prompt"
                base["notes"] = "numeric entry captured by its _number/_age sibling field."
                return base
            base["disposition"] = "numeric"
            base["notes"] = "numeric entry gated by a radio prompt; GWAS the entered value."
            return base
        # Degenerate: no usable contrast after removing missing-type answers.
        if n_core < 2:
            base["disposition"] = "excluded_degenerate"
            base["notes"] = "fewer than 2 informative answer options."
            if sens:
                base["flag_reason"] = "sensitive"
            return base
        if is_nominal(item, options):
            base["disposition"] = "nominal_binary"
            if sens:
                base["flag_reason"] = "sensitive"
            return base
        res = try_ordinal(item, options)
        if res is not None:
            rule, source, conf, mapping = res
            base["disposition"] = "ordinal_and_binary"
            base["ordinal_rule"] = rule
            base["ordinal_source"] = source
            base["ordinal_confidence"] = conf
            base["n_ordinal_levels"] = len(set(mapping.values()))
            base["_ordinal_map"] = mapping
            if conf == "medium" or sens:
                base["flag_reason"] = ";".join(
                    (["medium_confidence_ordinal"] if conf == "medium" else [])
                    + (["sensitive"] if sens else [])
                )
            return base
        # A 2-option single-select with no ordinal scale is just a binary
        # contrast (e.g. Yes/No); it is not an "uncertain ordinal" candidate.
        if n_core == 2:
            base["disposition"] = "binary_only"
            base["notes"] = "binary single-select; one-vs-rest GWAS."
            if sens:
                base["flag_reason"] = "sensitive"
            return base
        # >=3 options and no defensible ordinal scale: binary + review flag.
        base["disposition"] = "flagged_review"
        base["flag_reason"] = ";".join(["uncertain_ordinal"] + (["sensitive"] if sens else []))
        return base

    base["disposition"] = "excluded_other"
    return base


def main() -> None:
    here = Path(__file__).resolve().parent
    root = here.parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", type=Path, default=root / "metadata/survey_item_inventory.tsv")
    ap.add_argument("--outdir", type=Path, default=root / "metadata")
    args = ap.parse_args()

    with open(args.inventory) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # -- template library ----------------------------------------------------
    with open(args.outdir / "ordinal_answer_templates.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rule", "answer_label_normalized", "value", "confidence", "description"])
        for name, tpl in R.TEMPLATES.items():
            for lab, val in sorted(tpl["map"].items(), key=lambda kv: kv[1]):
                w.writerow([name, lab, val, tpl["confidence"], tpl["desc"]])

    # -- per-question manifest + ordinal mapping + flags ---------------------
    q_out = open(args.outdir / "survey_question_manifest.tsv", "w", newline="")
    qw = csv.writer(q_out, delimiter="\t")
    qw.writerow(
        [
            "survey",
            "item_concept",
            "question_concept_id",
            "field_type",
            "phenotype_class",
            "disposition",
            "ordinal_rule",
            "ordinal_source",
            "ordinal_confidence",
            "n_options",
            "n_binary_phenos",
            "n_ordinal_levels",
            "sensitive_topics",
            "flag_reason",
            "notes",
            "field_label",
        ]
    )

    o_out = open(args.outdir / "ordinal_mapping_manifest.tsv", "w", newline="")
    ow = csv.writer(o_out, delimiter="\t")
    ow.writerow(
        [
            "survey",
            "item_concept",
            "question_concept_id",
            "ordinal_rule",
            "ordinal_source",
            "confidence",
            "answer_label",
            "answer_label_normalized",
            "ordinal_value",
            "field_label",
        ]
    )

    f_out = open(args.outdir / "flagged_questions.tsv", "w", newline="")
    fw = csv.writer(f_out, delimiter="\t")
    fw.writerow(
        [
            "survey",
            "item_concept",
            "question_concept_id",
            "flag_reason",
            "sensitive_topics",
            "disposition",
            "ordinal_rule",
            "n_options",
            "field_label",
            "options",
        ]
    )

    from collections import Counter

    all_items = {r["item_concept"] for r in rows}
    counts = Counter()
    for row in rows:
        d = disposition(row, all_items)
        counts[d["disposition"]] += 1
        qw.writerow(
            [
                row["survey"],
                row["item_concept"],
                row["question_concept_id"],
                row["field_type"],
                row["phenotype_class"],
                d["disposition"],
                d["ordinal_rule"],
                d["ordinal_source"],
                d["ordinal_confidence"],
                row["n_options"],
                d["n_binary_phenos"],
                d["n_ordinal_levels"],
                d["sensitive_topics"],
                d["flag_reason"],
                d["notes"],
                row["field_label"],
            ]
        )
        if d["_ordinal_map"]:
            options = parse_options(row["options"])
            for _, label in options:
                nl = R.norm(label)
                if nl in d["_ordinal_map"]:
                    ow.writerow(
                        [
                            row["survey"],
                            row["item_concept"],
                            row["question_concept_id"],
                            d["ordinal_rule"],
                            d["ordinal_source"],
                            d["ordinal_confidence"],
                            label,
                            nl,
                            d["_ordinal_map"][nl],
                            row["field_label"],
                        ]
                    )
        if d["flag_reason"]:
            fw.writerow(
                [
                    row["survey"],
                    row["item_concept"],
                    row["question_concept_id"],
                    d["flag_reason"],
                    d["sensitive_topics"],
                    d["disposition"],
                    d["ordinal_rule"],
                    row["n_options"],
                    row["field_label"],
                    row["options"],
                ]
            )

    for h in (q_out, o_out, f_out):
        h.close()

    print("Disposition counts:")
    for k, v in counts.most_common():
        print(f"  {k:28s} {v}")


if __name__ == "__main__":
    main()
