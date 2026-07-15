#!/usr/bin/env python3
"""Parse the All of Us survey codebook workbook into a flat question inventory.

The workbook `survey_data_codebooks.xlsx` is a REDCap-style data dictionary.
Two physical layouts appear:

  * standard sheets (Basics, Lifestyle, ...): header on row 0, no Version column.
  * COPE / Minute Survey sheets: a leading "Version" column, and the Minute
    Survey header row is not the first row.

This module normalises both into one schema and parses the inline answer
"Choices" cell (``code, label | code, label | ...``) into structured options.

It emits one row per *real* survey item (radio/dropdown/checkbox/slider/yesno,
plus numeric-validated text) and skips descriptive/section and free-text items,
recording why each item was skipped so nothing is silently dropped.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import openpyxl


# Survey sheets in the workbook that hold real questionnaire items.
SURVEY_SHEETS = [
    "Basics",
    "Lifestyle",
    "Overall Health",
    "Family Health History",
    "Personal Medical History",
    "Personal and Family Health Hist",
    "Healthcare Access and Utilizati",
    "Social Determinants of Health",
    "COPE ",
    "Minute Survey on COVID-19 Vacci",
    "Life Functioning",
    "Emotional Health History",
    "Behavioral Health",
]

# Canonical survey display names (used across all downstream manifests).
SURVEY_DISPLAY = {
    "Basics": "The Basics",
    "Lifestyle": "Lifestyle",
    "Overall Health": "Overall Health",
    "Family Health History": "Family Health History",
    "Personal Medical History": "Personal Medical History",
    "Personal and Family Health Hist": "Personal and Family Health History",
    "Healthcare Access and Utilizati": "Healthcare Access and Utilization",
    "Social Determinants of Health": "Social Determinants of Health",
    "COPE ": "COVID-19 Participant Experience (COPE)",
    "Minute Survey on COVID-19 Vacci": "Minute Survey on COVID-19 Vaccines",
    "Life Functioning": "Life Functioning",
    "Emotional Health History": "Emotional Health History and Well-Being",
    "Behavioral Health": "Behavioral Health and Personality",
}

CANONICAL_COLUMNS = [
    "Item Concept",
    "Form Name",
    "Section Header",
    "Field Type",
    "Field Label",
    "Choices, Calculations, OR Slider Labels",
    "Text Validation Type OR Show Slider Number",
    "Text Validation Min",
    "Text Validation Max",
    "Branching Logic (Show field only if...)",
]


@dataclass
class Item:
    survey: str
    sheet: str
    version: str
    item_concept: str
    form_name: str
    section_header: str
    field_type: str
    field_label: str
    choices_raw: str
    validation_type: str
    validation_min: str
    validation_max: str
    branching_logic: str
    options: list[tuple[str, str]] = field(default_factory=list)


def _clean(value) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value).replace("\xa0", " ")).strip()


def _find_header_row(rows: list[tuple]) -> int:
    """Return the index of the header row (the one containing 'Item Concept')."""
    for idx, row in enumerate(rows):
        cells = [_clean(c).lower() for c in row]
        if "item concept" in cells and "field type" in cells:
            return idx
    raise ValueError("Header row not found")


def _column_index(header: list[str]) -> dict[str, int]:
    idx = {}
    for i, name in enumerate(header):
        idx.setdefault(_clean(name).lower(), i)
    return idx


def parse_choices(choices_raw: str) -> list[tuple[str, str]]:
    """Parse a REDCap choices cell into ``[(code, label), ...]``.

    Choices are pipe-delimited; within each choice the code and label are
    separated by the *first* comma. Codes never contain commas; labels may.
    """
    out: list[tuple[str, str]] = []
    text = _clean(choices_raw)
    if not text:
        return out
    for chunk in text.split("|"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "," in chunk:
            code, label = chunk.split(",", 1)
            out.append((code.strip(), label.strip()))
        else:
            out.append((chunk.strip(), ""))
    return out


def parse_sheet(ws) -> list[Item]:
    rows = list(ws.iter_rows(values_only=True))
    header_idx = _find_header_row(rows)
    header = [_clean(c) for c in rows[header_idx]]
    col = _column_index(header)
    has_version = "version" in col

    def get(row, name):
        i = col.get(name)
        if i is None or i >= len(row):
            return ""
        return _clean(row[i])

    items: list[Item] = []
    display = SURVEY_DISPLAY[ws.title]
    for row in rows[header_idx + 1 :]:
        if not any(c is not None and _clean(c) for c in row):
            continue
        item_concept = get(row, "item concept")
        field_type = get(row, "field type").lower()
        if not item_concept and not field_type:
            continue
        it = Item(
            survey=display,
            sheet=ws.title,
            version=get(row, "version") if has_version else "",
            item_concept=item_concept,
            form_name=get(row, "form name"),
            section_header=get(row, "section header"),
            field_type=field_type,
            field_label=get(row, "field label"),
            choices_raw=get(row, "choices, calculations, or slider labels"),
            validation_type=get(row, "text validation type or show slider number"),
            validation_min=get(row, "text validation min"),
            validation_max=get(row, "text validation max"),
            branching_logic=get(row, "branching logic (show field only if...)"),
        )
        it.options = parse_choices(it.choices_raw)
        items.append(it)
    return items


def load_all(xlsx_path: Path) -> list[Item]:
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    items: list[Item] = []
    for sheet in SURVEY_SHEETS:
        items.extend(parse_sheet(wb[sheet]))
    wb.close()
    return items


def load_question_concept_map(path: Path) -> dict[str, str]:
    """Map normalised question text -> question_concept_id from the AoU manifest."""
    out: dict[str, str] = {}
    if not path or not path.exists():
        return out
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            q = _normalize_text(row.get("question", ""))
            qid = row.get("question_concept_id", "").strip()
            if q and qid:
                out.setdefault(q, qid)
    return out


def _normalize_text(text: str) -> str:
    text = re.sub(r"<[^>]+>", " ", str(text or ""))
    text = re.sub(r"\s+", " ", text).strip().lower()
    return text


# Field types that are never GWAS phenotypes.
SKIP_FIELD_TYPES = {"descriptive", "calc"}


def classify(it: Item) -> str:
    ft = it.field_type
    if ft in SKIP_FIELD_TYPES:
        return "descriptive"
    if ft in {"radio", "dropdown", "yesno"}:
        return "single_select"
    if ft == "checkbox":
        return "multi_select"
    if ft == "slider":
        return "slider"
    if ft == "text":
        vt = it.validation_type.lower()
        if vt in {"number", "integer", "numeric"} or re.search(r"number|integer", vt):
            return "numeric_text"
        if re.search(r"date|datetime|time", vt):
            return "date_text"
        return "free_text"
    return "other:" + ft


def main() -> None:
    here = Path(__file__).resolve().parent
    root = here.parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", type=Path, default=root.parent.parent / "survey_data_codebooks.xlsx")
    ap.add_argument(
        "--concept-map",
        type=Path,
        default=root.parent / "data/aou_metadata/aou_ds_survey_question_concepts.tsv",
    )
    ap.add_argument("--out", type=Path, default=root / "metadata/survey_item_inventory.tsv")
    args = ap.parse_args()

    items = load_all(args.xlsx)
    concept_map = load_question_concept_map(args.concept_map)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "survey",
                "sheet",
                "version",
                "item_concept",
                "field_type",
                "phenotype_class",
                "question_concept_id",
                "n_options",
                "validation_type",
                "validation_min",
                "validation_max",
                "has_branching",
                "field_label",
                "options",
            ]
        )
        for it in items:
            qid = concept_map.get(_normalize_text(it.field_label), "")
            opt_str = " | ".join(f"{c}={l}" for c, l in it.options)
            w.writerow(
                [
                    it.survey,
                    it.sheet,
                    it.version,
                    it.item_concept,
                    it.field_type,
                    classify(it),
                    qid,
                    len(it.options),
                    it.validation_type,
                    it.validation_min,
                    it.validation_max,
                    "yes" if it.branching_logic else "no",
                    it.field_label,
                    opt_str,
                ]
            )
    print(f"Wrote {args.out} ({len(items)} items)")


if __name__ == "__main__":
    main()
