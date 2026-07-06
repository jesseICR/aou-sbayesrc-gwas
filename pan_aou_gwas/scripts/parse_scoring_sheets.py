#!/usr/bin/env python3
"""Parse the codebook 'Scoring' sheets into a validated-composite-score manifest.

The scoring sheets define validated instruments (GAD-7, PHQ-9, PSS, ASRS,
BFI-2-XS, UCLA loneliness, ACE, PROMIS, ...) as, per item: a question and a set
of answer->value rows. Most sheets bake reverse-keying into the per-item values;
BFI-2-XS lists raw 1-5 and needs the standard domain split + reverse (handled in
the composite builder, not here).

Output: metadata/composite_items_manifest.tsv with one row per
(instrument, item, answer) -> value, plus the item's question text (used to
match the OMOP survey response at runtime, the same way ordinal mapping works).
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import openpyxl

SCORING_SHEETS = {
    "COPE Scoring": "COVID-19 Participant Experience (COPE)",
    "Overall Health Scoring": "Overall Health",
    "Lifestyle Scoring": "Lifestyle",
    "Social Determinants of Health S": "Social Determinants of Health",
    "Emotional Health Scoring": "Emotional Health History and Well-Being",
    "Behavioral Health Scoring": "Behavioral Health and Personality",
}

HEADER_TOKENS = {"instrument:", "instrument", "measure", "source", ""}
ITEM_CODE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def _clean(v) -> str:
    return re.sub(r"\s+", " ", str(v).replace("\xa0", " ")).strip() if v is not None else ""


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def parse_sheet(ws, survey):
    rows_out = []
    current_instrument = ""
    current_item = None
    current_q = ""
    for row in ws.iter_rows(values_only=True):
        cells = [_clean(c) for c in row]
        if not any(cells):
            continue
        c0 = cells[0]
        if c0 and c0.lower() not in HEADER_TOKENS:
            current_instrument = c0
        # answer row: a numeric value somewhere in cells[2:]
        val = None
        for c in cells[2:]:
            n = _num(c)
            if n is not None:
                val = n
                break
        label = cells[1]
        if val is not None and label and current_item:
            rows_out.append((current_instrument, survey, current_item, current_q, label, val))
            continue
        # item header: cells[1] is a question, and there is an item-code token
        code = ""
        for c in reversed(cells[1:]):
            if c and ITEM_CODE_RE.match(c) and _num(c) is None:
                code = c
                break
        if label and code and code != label:
            current_item = code
            current_q = label
    return rows_out


def main() -> None:
    here = Path(__file__).resolve().parent
    root = here.parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", type=Path, default=root.parent.parent / "survey_data_codebooks.xlsx")
    ap.add_argument("--out", type=Path, default=root / "metadata/composite_items_manifest.tsv")
    args = ap.parse_args()

    wb = openpyxl.load_workbook(args.xlsx, read_only=True, data_only=True)
    allrows = []
    for sheet, survey in SCORING_SHEETS.items():
        allrows.extend(parse_sheet(wb[sheet], survey))
    wb.close()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["instrument", "survey", "item_code", "question", "answer_label", "value"])
        w.writerows(allrows)
    # summary
    from collections import defaultdict
    items = defaultdict(set)
    for inst, _s, code, _q, _l, _v in allrows:
        items[inst].add(code)
    print(f"Wrote {args.out} ({len(allrows)} answer-score rows)")
    for inst in sorted(items):
        print(f"  {inst:52s} {len(items[inst]):2d} items")


if __name__ == "__main__":
    main()
