#!/usr/bin/env python3
"""Build residualized phenotypes and run covariate-free PLINK2 GWAS on All of Us.

This is the worker behind run_pan_aou_gwas.sh. It consumes the extracted survey /
measurement CSVs and the phenotype manifests, then for every eligible phenotype:

  * binary one-vs-rest (linear probability model residual)   -> linear GWAS
  * ordinal (mapped -> IRNT -> residual)                     -> linear GWAS
  * numeric (parsed -> IRNT -> residual)                     -> linear GWAS
  * physical measurement (QC'd -> IRNT -> residual)          -> linear GWAS
  * PFHH self-history binary (LPM residual)                  -> linear GWAS

Covariates (age_c, sex_c, age_c:sex_c, PC1..PC10) are regressed out BEFORE PLINK2;
PLINK2 runs `--glm allow-no-covars`. See SPECSHEET.md. The residualize()/rint()
routines are ported verbatim from ~/projects/ukgwas/covariate_experiment so the
method matches the validated experiment.

Question<->manifest matching is by normalized question text, because the codebook
uses REDCap item concepts while OMOP carries question_concept_id + text; the
ordinal answer mappings are answer-text driven for the same reason.
"""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import json
import math
import os
import re
import sys
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ordinal_rules as R  # noqa: E402

PC_COLUMNS = [f"PC{i}_AVG" for i in range(1, 11)]

MIN_CASES = 200
MIN_CONTROLS = 200
MIN_QUANT_N = 500
MIN_ORDINAL_LEVELS = 3

MEASUREMENTS = {
    # measurement phenotype -> (concept_ids, unit_to_canonical, (lo, hi))
    "height_cm": ({903133}, "cm", (100.0, 250.0)),
    "systolic_bp_mmhg": ({903118}, "mmhg", (70.0, 260.0)),
    "diastolic_bp_mmhg": ({903115}, "mmhg", (40.0, 160.0)),
    "heart_rate_bpm": ({903126}, "bpm", (30.0, 220.0)),
    # BMI is derived from height+weight or taken from concept 903124; handled below.
}

ZIP3_SES_TRAITS = {
    "zip3_deprivation_index": (
        "deprivation_index",
        "ZIP3 Area SES: deprivation index",
    ),
    "zip3_median_income": (
        "median_income",
        "ZIP3 Area SES: median income",
    ),
    "zip3_fraction_poverty": (
        "fraction_poverty",
        "ZIP3 Area SES: fraction in poverty",
    ),
    "zip3_fraction_assisted_income": (
        "fraction_assisted_income",
        "ZIP3 Area SES: fraction receiving assisted income",
    ),
    "zip3_fraction_no_health_ins": (
        "fraction_no_health_ins",
        "ZIP3 Area SES: fraction without health insurance",
    ),
    "zip3_fraction_vacant_housing": (
        "fraction_vacant_housing",
        "ZIP3 Area SES: fraction vacant housing",
    ),
    "zip3_fraction_high_school_edu": (
        "fraction_high_school_edu",
        "ZIP3 Area SES: fraction high-school educated",
    ),
}


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def norm_q(text: str) -> str:
    text = re.sub(r"<[^>]+>", " ", str(text or ""))
    text = text.translate(str.maketrans({"‘": "'", "’": "'", "“": '"', "”": '"'}))
    text = re.sub(r"\s+", " ", text).strip().lower().strip(" ?.\"'")
    return text


def compact_key(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text or "").lower())


# --------------------------------------------------------------------------- #
# ported from covariate_experiment
# --------------------------------------------------------------------------- #
def rint(values: np.ndarray) -> np.ndarray:
    from scipy import stats

    if np.isnan(values).any():
        raise RuntimeError("RINT received missing values.")
    ranks = stats.rankdata(values, method="average")
    quantiles = (ranks - 0.5) / len(values)
    return stats.norm.ppf(quantiles)


def residualize(y: np.ndarray, covars: np.ndarray) -> np.ndarray:
    x = np.column_stack([np.ones(len(y)), covars])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


# --------------------------------------------------------------------------- #
# loaders
# --------------------------------------------------------------------------- #
def load_keep(path: Path) -> list[str]:
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if p:
                out.append(p[-1])
    return out


def load_sex(path: Path) -> dict[str, int]:
    out = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            iid = row["IID"].strip()
            try:
                out[iid] = int(row["sex_01"])
            except (KeyError, ValueError):
                pass
    return out


def load_pcs(path: Path) -> dict[str, list[float]]:
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header and header[0] == "#FID":
            header[0] = "FID"
        idx = {h: i for i, h in enumerate(header)}
        for h in ["IID", *PC_COLUMNS]:
            if h not in idx:
                raise ValueError(f"{path} missing column {h}")
        for line in f:
            fld = line.rstrip("\n").split("\t")
            try:
                out[fld[idx["IID"]]] = [float(fld[idx[h]]) for h in PC_COLUMNS]
            except (IndexError, ValueError):
                continue
    return out


def load_fam_fids(path: Path) -> dict[str, str]:
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                out[p[1]] = p[0]
    return out


def load_question_manifest(path: Path) -> tuple[dict[str, dict], dict[str, dict], list[dict]]:
    """Return manifest maps by normalized prompt text/qid and compact item_concept."""
    by_text = {}
    by_item = {}
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
            by_text.setdefault(norm_q(row["field_label"]), row)
            qid = (row.get("question_concept_id") or "").strip()
            if qid:
                by_text.setdefault(qid, row)
            item = compact_key(row.get("item_concept", ""))
            if item:
                by_item.setdefault(item, row)
    return by_text, by_item, rows


def load_live_question_crosswalk(path: Path | None, manifest_rows: list[dict], by_text: dict[str, dict]):
    """Map live AoU question_concept_id and item_concept to manifest rows.

    AoU ds_survey uses short labels such as "Overall Health: General Health",
    while the codebook manifest stores full participant-facing prompts. This
    crosswalk links live question IDs back to codebook item_concept values using
    the repo's metadata export.
    """
    qid_to_manifest = {}
    item_to_qid = {}
    if not path or not path.exists():
        return qid_to_manifest, item_to_qid

    rows_by_survey = defaultdict(list)
    for row in manifest_rows:
        rows_by_survey[row.get("survey", "")].append(row)

    with open(path, newline="") as f:
        for live in csv.DictReader(f, delimiter="\t"):
            qid = (live.get("question_concept_id") or "").strip()
            if not qid:
                continue
            live_text = live.get("question", "")
            live_norm = norm_q(live_text)
            man = by_text.get(live_norm)
            if man is None:
                live_compact = compact_key(live_text)
                candidates = []
                for row in rows_by_survey.get(live.get("survey", ""), []):
                    item_compact = compact_key(row.get("item_concept", ""))
                    if not item_compact or not live_compact:
                        continue
                    if (
                        item_compact == live_compact
                        or item_compact.endswith(live_compact)
                        or live_compact.endswith(item_compact)
                    ):
                        candidates.append(row)
                if len(candidates) == 1:
                    man = candidates[0]
            if man is None:
                continue
            qid_to_manifest[qid] = man
            item = (man.get("item_concept") or "").strip()
            if item:
                item_to_qid[item] = qid
    return qid_to_manifest, item_to_qid


def load_ea_proxy_feature_sources(path: Path | None, existing_qman: dict[str, dict]):
    """Supplement metadata with EA-proxy source questions absent from codebooks.

    The pan-AoU metadata is primarily REDCap/codebook-derived. The SES-EA XGBoost
    feature contract was built directly from live AoU question IDs, and a subset
    of those v9 IDs do not round-trip through the codebook text matcher. This
    supplemental manifest adds only missing live question IDs, leaving codebook
    rows in charge wherever both sources exist.
    """
    rows = []
    if not path or not path.exists():
        return rows
    encoding_to_disposition = {
        "ordinal": "ordinal_and_binary",
        "one_hot": "nominal_binary",
        "allowlisted_one_hot": "nominal_binary",
        "numeric": "numeric",
    }
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (row.get("include_exclude") or "").strip() != "include":
                continue
            qid = (row.get("question_concept_id") or "").strip()
            if not qid.isdigit() or qid in existing_qman:
                continue
            encoding = (row.get("encoding") or "").strip()
            disp = encoding_to_disposition.get(encoding)
            if not disp:
                continue
            rows.append({
                "survey": row.get("survey", ""),
                "item_concept": f"xgb_q{qid}",
                "question_concept_id": qid,
                "field_type": "text",
                "phenotype_class": "numeric_text" if encoding == "numeric" else "single_select",
                "disposition": disp,
                "ordinal_rule": "ea_proxy_ordinal_text" if encoding == "ordinal" else "",
                "ordinal_source": "ea_proxy_feature_manifest" if encoding == "ordinal" else "",
                "ordinal_confidence": "medium" if encoding == "ordinal" else "",
                "n_options": "",
                "n_binary_phenos": "",
                "n_ordinal_levels": "",
                "sensitive_topics": "",
                "flag_reason": "",
                "notes": "Supplemental source question from SES-EA XGBoost feature contract.",
                "field_label": row.get("item_name", ""),
            })
    return rows


def live_question_override_spec(qid: str, survey: str, question: str) -> dict | None:
    """Metadata overrides for live AoU qids absent from codebook-derived rows.

    These are deliberately narrow.  They cover v9 question_concept_id values
    observed in the live survey extracts that failed to round-trip through the
    REDCap/codebook crosswalk, while keeping excluded/admin fields explicit.
    """
    # Repeated COVID-19 vaccine follow-up fields.  The dose-specific "received
    # another dose", dose-type, vaccine-name, and adverse-reaction questions are
    # closed categorical fields; one-vs-rest binary phenotypes are appropriate.
    vaccine_additional_dose = {
        "765937",
        *{str(x) for x in range(765953, 765967)},
    }
    vaccine_adverse_reactions = {str(x) for x in range(765973, 766002, 2)}
    vaccine_type = {str(x) for x in range(766007, 766036, 2)}
    vaccine_name = {str(x) for x in range(766037, 766066, 2)}

    face_mask_qids = {"1310051", "1310052", "1310053", "1310056", "1310060", "1310062"}
    cope_phq_qids = {
        "1333274", "1333275", "1333276", "1333277", "1333278",
        "1333279", "1333280", "1333281", "1333285",
    }
    ehhw_phq_gad_qids = {
        "1703920", "1703977", "1703983", "1703984", "1703987",
        "1703995", "1703996", "1704000", "1704004", "1704024",
        "1704026", "1704028", "1704038", "1704039", "1704041",
        "1704042",
    }
    ehhw_cidi_frequency_qids = {
        "1703979", "1703982", "1704006", "1704032", "1704040",
        "1704043", "1704050", "1704052", "1704053",
    }
    everyday_discrimination_qids = {
        "40192380", "40192395", "40192416", "40192451", "40192466",
        "40192489", "40192490", "40192496", "40192519",
    }
    health_care_discrimination_qids = {
        "40192383", "40192394", "40192423", "40192425",
        "40192497", "40192503", "40192505",
    }

    def spec(disposition: str, ordinal_rule: str = "", phenotype_class: str = "single_select",
             field_type: str = "radio", notes: str = "") -> dict:
        return {
            "survey": survey,
            "item_concept": f"live_q{qid}",
            "question_concept_id": qid,
            "field_type": field_type,
            "phenotype_class": phenotype_class,
            "disposition": disposition,
            "ordinal_rule": ordinal_rule,
            "ordinal_source": "live_qid_override" if ordinal_rule else "",
            "ordinal_confidence": "medium" if ordinal_rule else "",
            "n_options": "",
            "n_binary_phenos": "",
            "n_ordinal_levels": "",
            "sensitive_topics": "",
            "flag_reason": "",
            "notes": notes or "Live v9 qid override for question absent from codebook crosswalk.",
            "field_label": question,
        }

    if qid == "713888":
        return spec("ordinal_and_binary", "symptom_onset_month_ordinal")
    if qid in vaccine_additional_dose:
        return spec("nominal_binary")
    if qid in vaccine_adverse_reactions:
        return spec("nominal_binary", phenotype_class="multi_select", field_type="checkbox")
    if qid in vaccine_type or qid in vaccine_name:
        return spec("nominal_binary")
    if qid == "905040":
        return spec("nominal_binary", phenotype_class="multi_select", field_type="checkbox")
    if qid in face_mask_qids:
        return spec("ordinal_and_binary", "freq_never_always_na_0_3")
    if qid == "1310067":
        return spec("nominal_binary", phenotype_class="multi_select", field_type="checkbox")
    if qid == "1310133":
        return spec("nominal_binary")
    if qid == "1310136":
        return spec("ordinal_and_binary", "remote_childcare_0_3")
    if qid == "1310140":
        return spec("nominal_binary")
    if qid in {"1310144", "1310145"}:
        return spec("nominal_binary", phenotype_class="multi_select", field_type="checkbox")
    if qid == "1332748":
        return spec("ordinal_and_binary", "days_last5_midpoint")
    if qid in cope_phq_qids or qid in ehhw_phq_gad_qids:
        return spec("ordinal_and_binary", "phq_gad_0_3")
    if qid in {"1333286", "1333288", "1333289"}:
        return spec("nominal_binary")
    if qid == "1333287":
        return spec("excluded_descriptive", notes="Branching unit selector for split sitting-time numeric fields.")
    if qid == "1333300":
        return spec("ordinal_and_binary", "hygiene_frequency_1_4")
    if qid == "1333301":
        return spec("ordinal_and_binary", "days_last5_midpoint")
    if qid in {"903641", "903642"}:
        return spec("numeric", notes="Sitting time component; numeric field.")
    if qid in {"1703882", "1703886", "1703895", "1703923"}:
        return spec("binary_only")
    if qid in ehhw_cidi_frequency_qids:
        return spec("ordinal_and_binary", "time_all_none_0_4")
    if qid in {"1704015", "1704046"}:
        return spec("nominal_binary", phenotype_class="multi_select", field_type="checkbox")
    if qid == "1704030":
        return spec("binary_only")
    if qid == "1704135":
        return spec("numeric", notes="range [2.0, 99.0]")
    if qid in everyday_discrimination_qids:
        return spec("ordinal_and_binary", "freq_event_0_5")
    if qid in health_care_discrimination_qids:
        return spec("ordinal_and_binary", "freq_never_always_0_4")
    if qid == "40192428":
        return spec("nominal_binary", phenotype_class="multi_select", field_type="checkbox")
    if qid == "43529899":
        return spec("ordinal_and_binary", "freq_always_none_0_3")
    if qid == "43529901":
        return spec("ordinal_and_binary", "importance_0_3")
    if qid == "43529902":
        return spec("ordinal_and_binary", "freq_always_none_0_3")
    return None


def load_live_question_overrides(path: Path | None, existing_qman: dict[str, dict]):
    """Build manifest rows for selected live qids absent from the codebook map."""
    rows = []
    seen_qids = set()
    if not path or not path.exists():
        pass
    else:
        with open(path, newline="") as f:
            for live in csv.DictReader(f, delimiter="\t"):
                qid = (live.get("question_concept_id") or "").strip()
                if not qid or qid in existing_qman:
                    continue
                row = live_question_override_spec(
                    qid,
                    (live.get("survey") or "").strip(),
                    (live.get("question") or "").strip(),
                )
                if row:
                    rows.append(row)
                    seen_qids.add(qid)
    static_live_questions = {
        "836838": ("Personal and Family Health History", "Who in your family has had a kidney condition? Select all that apply."),
        "1703882": ("Behavioral Health and Personality", "Did you ever talk to a health professional about any of these experiences?"),
        "1703886": ("Behavioral Health and Personality", "Were you ever in your life frightened by agoraphobia-like situations?"),
        "1703895": ("Behavioral Health and Personality", "Have you ever been bothered with recurring intrusive thoughts, images, or urges?"),
        "1703923": ("Behavioral Health and Personality", "Have you ever had a high, excited, or hyper period that others noticed or that got you into trouble?"),
        "1703920": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by becoming easily annoyed or irritable"),
        "1703977": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by thoughts that you would be better off dead or of hurting yourself in some way"),
        "1703979": ("Emotional Health History and Well-Being", "During those 6 months, how often did you worry excessively or too much?"),
        "1703982": ("Emotional Health History and Well-Being", "During those 6 months, how often did you have trouble controlling your worry?"),
        "1703983": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by trouble falling or staying asleep, or sleeping too much"),
        "1703984": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by feeling nervous, anxious, or on edge"),
        "1703987": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by trouble relaxing"),
        "1703995": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by not being able to stop or control worrying"),
        "1703996": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by moving or speaking slowly or being fidgety/restless"),
        "1704000": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by worrying too much about different things"),
        "1704004": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by poor appetite or overeating"),
        "1704006": ("Emotional Health History and Well-Being", "During those 6 months, how often did you worry about a number of different things in your life?"),
        "1704015": ("Emotional Health History and Well-Being", "During feelings of depression or loss of interest did you ever try any of the following for these problems?"),
        "1704024": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by feeling down, depressed, or hopeless"),
        "1704026": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by little interest or pleasure in doing things"),
        "1704028": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by being so restless that it is hard to sit still"),
        "1704030": ("Emotional Health History and Well-Being", "During feelings of depression or loss of interest did you ever tell a professional about these problems?"),
        "1704032": ("Emotional Health History and Well-Being", "During those 6 months, how often did you have difficulty concentrating or your mind going blank?"),
        "1704038": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by trouble concentrating?"),
        "1704039": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by feeling tired or having little energy"),
        "1704040": ("Emotional Health History and Well-Being", "During those 6 months, how often did you have difficulty falling or staying asleep or have restless sleep?"),
        "1704041": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by feeling bad about yourself or that you are a failure"),
        "1704042": ("Emotional Health History and Well-Being", "Over the last 2 weeks, how often have you been bothered by feeling afraid as if something awful might happen"),
        "1704043": ("Emotional Health History and Well-Being", "During those 6 months, how often did you have muscle aches or tension?"),
        "1704046": ("Emotional Health History and Well-Being", "During feelings of depression or loss of interest did you try any of the following medications?"),
        "1704050": ("Emotional Health History and Well-Being", "During those 6 months, how often did you feel restless, keyed up, or on edge?"),
        "1704052": ("Emotional Health History and Well-Being", "During those 6 months, how often did you feel worried and anxious?"),
        "1704053": ("Emotional Health History and Well-Being", "During those 6 months, how often did you feel irritated, annoyed, or grouchy?"),
        "1704135": ("Emotional Health History and Well-Being", "During those 6 months, about how old were you when you first began having problems with anxiety or worrying?"),
    }
    for qid, (survey, question) in static_live_questions.items():
        if qid in existing_qman or qid in seen_qids:
            continue
        if qid == "836838":
            row = {
                "survey": survey,
                "item_concept": f"live_q{qid}",
                "question_concept_id": qid,
                "field_type": "checkbox",
                "phenotype_class": "multi_select",
                "disposition": "excluded_family_history",
                "ordinal_rule": "",
                "ordinal_source": "",
                "ordinal_confidence": "",
                "n_options": "",
                "n_binary_phenos": "",
                "n_ordinal_levels": "",
                "sensitive_topics": "",
                "flag_reason": "",
                "notes": "Handled by PFHH relatedness-burden allowlist.",
                "field_label": question,
            }
        else:
            row = live_question_override_spec(qid, survey, question)
        if row:
            rows.append(row)
    return rows


def load_ordinal_lookup(path: Path) -> dict[tuple[str, str], float]:
    """(question key or item_concept, normalized answer label) -> ordinal value."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ans = R.norm(row["answer_label"])
            value = float(row["ordinal_value"])
            out[(norm_q(row["field_label"]), ans)] = value
            item = (row.get("item_concept") or "").strip()
            if item:
                out[(item, ans)] = value
    return out


# --------------------------------------------------------------------------- #
# survey ingest: latest response per (person, question)
# --------------------------------------------------------------------------- #
def read_survey_rows(
    paths: list[Path],
    keep: set[str],
    allowed_qids: set[str] | None = None,
    allowed_question_texts: set[str] | None = None,
):
    """Yield dict rows for retained samples from one or more survey CSVs."""
    allowed_qids = allowed_qids or set()
    allowed_question_texts = allowed_question_texts or set()
    filter_questions = bool(allowed_qids or allowed_question_texts)
    for path in paths:
        if not path or not path.exists() or path.stat().st_size == 0:
            continue
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                pid = (row.get("person_id") or "").strip()
                if pid not in keep:
                    continue
                if filter_questions:
                    qid = (row.get("question_concept_id") or "").strip()
                    if qid not in allowed_qids and norm_q(row.get("question") or "") not in allowed_question_texts:
                        continue
                yield row


def build_latest_responses(
    survey_paths,
    keep,
    allowed_qids: set[str] | None = None,
    allowed_question_texts: set[str] | None = None,
):
    """Return dict[qid] -> {question, pid -> (age, [(ans_text)...])} using latest datetime."""
    # Per (pid, qid): keep max datetime, collect answers at that datetime.
    # The full v9 survey export is large enough that storing string tuple keys
    # can exceed default AoU Jupyter RAM. Use compact integer ids internally and
    # expand back to strings only after aggregation completes.
    key_mult = 100000
    pid_by_index = list(keep)
    pid_to_index = {pid: i for i, pid in enumerate(pid_by_index)}
    qid_to_index = {}
    qid_by_index = []
    question_text_by_qid_index = []
    answer_to_index = {}
    answer_by_index = []
    latest = {}  # compact_key -> [datetime, age, answer_index | list[answer_index]]

    def qid_index(qid: str, question_text: str) -> int:
        idx = qid_to_index.get(qid)
        if idx is None:
            idx = len(qid_by_index)
            if idx >= key_mult:
                raise RuntimeError(f"Too many survey question IDs for key multiplier {key_mult}")
            qid_to_index[qid] = idx
            qid_by_index.append(qid)
            question_text_by_qid_index.append(question_text)
        elif question_text and not question_text_by_qid_index[idx]:
            question_text_by_qid_index[idx] = question_text
        return idx

    def answer_index(answer: str) -> int:
        idx = answer_to_index.get(answer)
        if idx is None:
            idx = len(answer_by_index)
            answer_to_index[answer] = idx
            answer_by_index.append(answer)
        return idx

    for row in read_survey_rows(survey_paths, keep, allowed_qids, allowed_question_texts):
        pid = row["person_id"].strip()
        pid_idx = pid_to_index.get(pid)
        if pid_idx is None:
            continue
        qid = sys.intern((row.get("question_concept_id") or "").strip())
        if not qid:
            continue
        dt = sys.intern((row.get("survey_datetime") or "").strip())
        ans_idx = answer_index(sys.intern((row.get("answer") or "").strip()))
        try:
            age = float(row.get("age_at_survey") or "nan")
        except ValueError:
            age = float("nan")
        qid_idx = qid_index(qid, row.get("question") or "")
        k = pid_idx * key_mult + qid_idx
        cur = latest.get(k)
        if cur is None or dt > cur[0]:
            latest[k] = [dt, age, ans_idx]
        elif dt == cur[0]:
            answers = cur[2]
            if isinstance(answers, int):
                if ans_idx != answers:
                    cur[2] = [answers, ans_idx]
            elif ans_idx not in answers:
                answers.append(ans_idx)

    questions = defaultdict(lambda: {"question": "", "responses": {}})
    while latest:
        k, (_dt, age, answer_indexes) = latest.popitem()
        pid_idx, qid_idx = divmod(k, key_mult)
        pid = pid_by_index[pid_idx]
        qid = qid_by_index[qid_idx]
        q = questions[qid]
        if not q["question"]:
            q["question"] = question_text_by_qid_index[qid_idx]
        if isinstance(answer_indexes, int):
            answers = (answer_by_index[answer_indexes],)
        else:
            answers = tuple(answer_by_index[i] for i in answer_indexes)
        q["responses"][pid] = (age, answers)
    return questions


def load_tsv_column_values(path: Path | None, column: str) -> set[str]:
    values = set()
    if not path or not path.exists():
        return values
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            value = (row.get(column) or "").strip()
            if value:
                values.add(value)
    return values


# --------------------------------------------------------------------------- #
# phenotype builders -> {"kind": binary|quant, "values": {iid: (y, age)}}
# --------------------------------------------------------------------------- #
def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")[:40] or "x"


def answer_tail(text: str) -> str:
    text = str(text or "").strip()
    if ":" in text:
        _prefix, tail = text.split(":", 1)
        if tail.strip():
            return tail.strip()
    return text


def answer_norm(text: str) -> str:
    return R.norm(answer_tail(text))


def answer_slug(text: str) -> str:
    return slug(answer_tail(text))


def is_missing_answer(text: str) -> bool:
    return R.is_missing(text) or R.is_missing(answer_tail(text))


def manifest_item_id(man: dict, fallback: str) -> str:
    """Stable codebook-backed identifier used in phenotype IDs."""
    item = (man.get("item_concept") or "").strip()
    return slug(item or fallback)


def ordinal_value_from_rule(rule: str, answer: str) -> float | None:
    """Map an answer to a value through a named ordinal template/override."""
    rule = (rule or "").strip()
    if not rule:
        return None
    ans = answer_norm(answer)
    candidates = [ans]
    if ";" in ans:
        candidates.append(ans.split(";", 1)[0].strip())
    template = R.TEMPLATES.get(rule)
    if template:
        if any(a in template.get("local_missing", set()) for a in candidates):
            return None
        for a in candidates:
            if a in template["map"]:
                return float(template["map"][a])
    for override in R.ITEM_OVERRIDES.values():
        if override.get("rule") != rule:
            continue
        if any(a in override.get("local_missing", set()) for a in candidates):
            return None
        for a in candidates:
            if a in override["map"]:
                return float(override["map"][a])
    if rule == "importance_0_3":
        # v9 live HCAU labels use "Not Important", while the codebook override
        # says "Not important at all".
        if ans == "not important":
            return 0.0
    return None


LIVE_ORDINAL_VALUE_BY_ITEM_ANSWER = {
    "educationlevel_highestgrade": {
        "never attended": 9.0,
        "one through four": 9.0,
        "five through eight": 9.0,
        "nine through eleven": 10.0,
        "twelve or ged": 13.0,
        "college one to three": 15.0,
        "college graduate": 18.0,
        "advanced degree": 20.0,
    },
    "income_annualincome": {
        "less 10k": 5.0,
        "10k 25k": 17.5,
        "25k 35k": 30.0,
        "35k 50k": 42.5,
        "50k 75k": 62.5,
        "75k 100k": 87.5,
        "100k 150k": 125.0,
        "150k 200k": 175.0,
        "more 200k": 250.0,
    },
}


def ea_proxy_ordinal_value_from_answer(answer: str):
    """Ordinal parser used by setup_ses_ea_proxy_gwas.py for XGBoost inputs."""
    tail = answer_norm(answer)
    exact = {
        "no": 0.0,
        "yes": 1.0,
        "never": 0.0,
        "rarely": 1.0,
        "sometimes": 2.0,
        "often": 3.0,
        "very often": 4.0,
        "not at all": 0.0,
        "some days": 1.0,
        "every day": 2.0,
        "never in last year": 0.0,
        "less than monthly": 1.0,
        "monthly": 2.0,
        "weekly": 3.0,
        "daily": 4.0,
        "monthly or less": 1.0,
        "2 to 4 per month": 2.0,
        "2 to 3 per week": 3.0,
        "4 or more per week": 4.0,
        "disagree strongly": 1.0,
        "disagree a little": 2.0,
        "neutral; no opinion": 3.0,
        "agree a little": 4.0,
        "agree strongly": 5.0,
        "strongly disagree": 1.0,
        "disagree": 2.0,
        "neither agree nor disagree": 3.0,
        "agree": 4.0,
        "strongly agree": 5.0,
        "poor": 1.0,
        "fair": 2.0,
        "good": 3.0,
        "very good": 4.0,
        "excellent": 5.0,
        "not at all confident": 1.0,
        "a little bit confident": 2.0,
        "somewhat confident": 3.0,
        "quite a bit confident": 4.0,
        "extremely confident": 5.0,
        "unable to do": 0.0,
        "with much difficulty": 1.0,
        "with some difficulty": 2.0,
        "with a little difficulty": 3.0,
        "without any difficulty": 4.0,
    }
    if tail in exact:
        return exact[tail]
    if tail.startswith("never"):
        return 0.0
    if tail.startswith("less than"):
        return 1.0
    if tail.startswith("monthly"):
        return 2.0
    if tail.startswith("weekly"):
        return 3.0
    if tail.startswith("daily"):
        return 4.0
    return None


# Minimum age-at-survey by ordinal rule, matching the repo's dedicated EA/income
# GWAS (setup_ea_gwas.py / setup_income_gwas.py, --min-age-at-survey default 26):
# exclude respondents who may not have completed education / are early-career.
MIN_AGE_BY_RULE = {
    "education_years_ea_proxy": 26.0,
    "income_midpoint_k": 26.0,
}


def build_survey_phenotypes(questions, qman, ord_lookup):
    """Yield (pheno_id, trait_type, kind, {iid: (y, age)}, meta)."""
    for qid, q in questions.items():
        qtext = q["question"]
        man = qman.get(qid) or qman.get(norm_q(qtext))
        if man is None:
            continue  # question not in our included/classified manifest
        disp = man["disposition"]
        if disp.startswith("excluded") or disp == "numeric":
            # numeric handled separately below via value_as_number
            if disp != "numeric":
                continue
        responses = q["responses"]

        # Age-at-survey minimum for EA/income (matches the repo GWAS scripts),
        # applied to every phenotype derived from that question.
        min_age = MIN_AGE_BY_RULE.get(man["ordinal_rule"])
        if min_age is not None:
            responses = {
                pid: (age, ans)
                for pid, (age, ans) in responses.items()
                if age is not None and not math.isnan(age) and age >= min_age
            }

        is_multi = man["phenotype_class"] == "multi_select"
        item_id = manifest_item_id(man, qid)
        # ---- binary one-vs-rest per observed valid answer ------------------
        if disp in ("ordinal_and_binary", "binary_only", "nominal_binary", "flagged_review"):
            # collect the valid (non-missing) answer universe for this question
            valid_answers = set()
            for _, (_, answers) in responses.items():
                for a in answers:
                    if not is_missing_answer(a):
                        valid_answers.add(a)
            for ans in sorted(valid_answers):
                values = {}
                for pid, (age, answers) in responses.items():
                    non_missing = {a for a in answers if not is_missing_answer(a)}
                    if not non_missing:
                        continue
                    if ans in non_missing:
                        values[pid] = (1.0, age)
                    else:
                        # single-select control = answered another valid option;
                        # checkbox control = question shown and option not selected.
                        values[pid] = (0.0, age)
                pid_ = f"bin_{item_id}__{answer_slug(ans)}"
                yield pid_, "binary", "binary", values, {
                    "question_concept_id": qid,
                    "item_concept": man.get("item_concept", ""),
                    "question": qtext,
                    "answer": ans,
                    "ordinal_rule": "",
                }
        # ---- ordinal -------------------------------------------------------
        if disp == "ordinal_and_binary":
            values = {}
            for pid, (age, answers) in responses.items():
                non_missing = [a for a in answers if not is_missing_answer(a)]
                if len(non_missing) != 1:
                    continue
                item = man.get("item_concept", "")
                ans_key = answer_norm(non_missing[0])
                v = LIVE_ORDINAL_VALUE_BY_ITEM_ANSWER.get(item, {}).get(ans_key)
                if v is None:
                    v = ord_lookup.get((item, ans_key))
                if v is None:
                    v = ord_lookup.get((norm_q(qtext), ans_key))
                if v is None:
                    v = ordinal_value_from_rule(man.get("ordinal_rule", ""), non_missing[0])
                if v is None and man.get("ordinal_rule") == "ea_proxy_ordinal_text":
                    v = ea_proxy_ordinal_value_from_answer(non_missing[0])
                if v is not None:
                    values[pid] = (float(v), age)
            yield f"ord_{item_id}", "ordinal", "quant", values, {
                "question_concept_id": qid,
                "item_concept": man.get("item_concept", ""),
                "question": qtext,
                "answer": "",
                "ordinal_rule": man["ordinal_rule"],
            }


def build_numeric_phenotypes(questions, qman):
    for qid, q in questions.items():
        man = qman.get(qid) or qman.get(norm_q(q["question"]))
        if man is None or man["disposition"] != "numeric":
            continue
        item_id = manifest_item_id(man, qid)
        # parse validation range "range [lo, hi]" from notes if present
        lo, hi = None, None
        m = re.search(r"range \[([^,]+),\s*([^\]]+)\]", man.get("notes", ""))
        if m:
            try:
                lo = float(m.group(1))
                hi = float(m.group(2))
            except ValueError:
                pass
        values = {}
        for pid, (age, answers) in q["responses"].items():
            nums = []
            for a in answers:
                try:
                    nums.append(float(a))
                except (TypeError, ValueError):
                    continue
            if len(nums) != 1:
                continue
            v = nums[0]
            if lo is not None and (v < lo or v > hi):
                continue
            values[pid] = (v, age)
        yield f"num_{item_id}", "numeric", "quant", values, {
            "question_concept_id": qid,
            "item_concept": man.get("item_concept", ""),
            "question": q["question"],
            "answer": "",
            "ordinal_rule": "",
        }


# Group -> PFHH category-screen question_concept_id used to recover controls.
PFHH_SCREEN_QID = {
    "Brain and nervous system": "43529272",
    "Mental health or substance use": "43529217",
    "Added skeletal/pain/injury": "702786",
    "Kidney conditions": "43529158",
}

# Family-history burden weights = coefficient of relationship (r) to the
# participant: self = 1, first-degree (parent/sibling/child) = 0.5,
# second-degree (grandparent, half-sibling) = 0.25. The burden phenotype is the
# sum of these weights over the relations the participant selected for a
# condition, so it is a within-family genetic-liability proxy for that condition.
PFHH_RELATION_WEIGHT = {
    "self": 1.0,
    "father": 0.5,
    "mother": 0.5,
    "parent": 0.5,
    "sibling": 0.5,
    "brother": 0.5,
    "sister": 0.5,
    "son": 0.5,
    "daughter": 0.5,
    "child": 0.5,
    "grandparent": 0.25,
    "grandmother": 0.25,
    "grandfather": 0.25,
    "half-sibling": 0.25,
    "half sibling": 0.25,
    "none": 0.0,  # explicit "no one" -> 0, distinct from unknown
}


def pfhh_relation_norm(answer: str) -> str:
    """Normalize PFHH checkbox answers to relation labels.

    Live AoU PFHH answer concepts can arrive as full labels like
    "<question text> - Self" rather than just "Self".  The PFHH burden scorer
    needs the relation suffix to match PFHH_RELATION_WEIGHT.
    """
    raw = str(answer or "").strip()
    if ":" in raw:
        prefix = raw.split(":", 1)[0].strip()
        first = prefix.split()[0].lower() if prefix.split() else ""
        if first in PFHH_RELATION_WEIGHT:
            return R.norm(first)
    text = answer_tail(answer)
    if " - " in text:
        text = text.rsplit(" - ", 1)[1]
    return R.norm(text)


def is_missing_pfhh_answer(answer: str) -> bool:
    return R.is_missing(answer) or R.is_missing(answer_tail(answer)) or R.is_missing(pfhh_relation_norm(answer))


def _pfhh_screen_completers(questions):
    out = {}
    for grp, sqid in PFHH_SCREEN_QID.items():
        pids = {}
        q = questions.get(sqid)
        if q:
            for pid, (age, answers) in q["responses"].items():
                if any(not is_missing_pfhh_answer(a) for a in answers):
                    pids[pid] = age
        out[grp] = pids
    return out


def build_pfhh_phenotypes(questions, allowlist_path):
    """Yield, per allowlisted PFHH condition, TWO phenotypes:

      pfhh_self_has_<cond>   binary  self-only
          case    = participant selected "Self"
          control = completed the relevant category screen, did not self-report
      pfhh_burden_<cond>     quant   genetic-relatedness-weighted family burden
          score   = sum of PFHH_RELATION_WEIGHT over selected relations
                    (self 1, first-degree 0.5, grandparent 0.25); "None"/not
                    shown after completing the screen -> 0; PMI-only -> missing

    Broad family-history phenotypes and non-self relative indicators are never
    run on their own; the burden score is an aggregate liability proxy, not a
    phenotype about any individual relative.
    """
    if not allowlist_path or not Path(allowlist_path).exists():
        return
    with open(allowlist_path, newline="") as f:
        allow = list(csv.DictReader(f, delimiter="\t"))

    screen_completers = _pfhh_screen_completers(questions)

    for row in allow:
        qid = row["question_concept_id"].strip()
        grp = row["pfhh_group"].strip()
        pheno_id = row["phenotype_id"].strip()
        cond = pheno_id.replace("self_has_", "")
        q = questions.get(qid)
        if q is None:
            continue

        # --- binary self_has ------------------------------------------------
        cases = {}
        for pid, (age, answers) in q["responses"].items():
            if any(pfhh_relation_norm(a) == "self" for a in answers):
                cases[pid] = age
        bin_values = {pid: (1.0, age) for pid, age in cases.items()}
        for pid, age in screen_completers.get(grp, {}).items():
            if pid not in cases:
                bin_values[pid] = (0.0, age)
        yield f"pfhh_{pheno_id}", "binary", "binary", bin_values, {
            "question_concept_id": qid,
            "question": row.get("question", ""),
            "answer": "Self",
            "ordinal_rule": "",
        }

        # --- quantitative family-burden sumscore ---------------------------
        burden = {pid: (0.0, age) for pid, age in screen_completers.get(grp, {}).items()}
        for pid, (age, answers) in q["responses"].items():
            non_missing = [a for a in answers if not is_missing_pfhh_answer(a)]
            if not non_missing:
                burden.pop(pid, None)  # answered only PMI -> drop from denominator
                continue
            score = sum(PFHH_RELATION_WEIGHT.get(pfhh_relation_norm(a), 0.0) for a in non_missing)
            burden[pid] = (score, age)
        yield f"pfhh_burden_{cond}", "pfhh_sumscore", "quant", burden, {
            "question_concept_id": qid,
            "question": row.get("question", ""),
            "answer": "relatedness-weighted burden (self=1, 1st-deg=0.5, grandparent=0.25)",
            "ordinal_rule": "pfhh_relatedness_burden",
        }


def build_composite_phenotypes(
    questions, manifest_path, ordinal_manifest_path=None, ord_lookup=None, qid_by_item=None
):
    """Yield validated composite scores (GAD-7, PHQ-9, PSS, BFI-2 Big Five, ...).

    Each composite is a prorated sum over its items (matched to survey responses
    by question text, reverse-keyed per composite_rules), requiring >= 80% of
    items answered. Residualized on the full covariate set (survey age is known).

    Scoring-sheet instruments come from `manifest_path`. The mixed-valence
    neighborhood/walkability/hunger composites are built from EXPLICIT_COMPOSITES
    (item lists), reusing the ordinal scores (`ord_lookup` + `ordinal_manifest`).
    """
    import composite_rules as CR

    rows = []
    if manifest_path and Path(manifest_path).exists():
        with open(manifest_path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

    qid_by_item = qid_by_item or {}
    qtext_to_qid = {}
    for qid, q in questions.items():
        qtext_to_qid.setdefault(norm_q(q["question"]), qid)

    # per instrument slug: {norm_question: {norm_answer: value}} merged across
    # administrations; plus per item_code answer maps for BFI-2.
    inst_items = defaultdict(dict)
    code_answer = defaultdict(dict)
    code_q = {}
    for r in rows:
        try:
            val = float(r["value"])
        except (ValueError, TypeError):
            continue
        ans = R.norm(r["answer_label"])
        code_answer[r["item_code"]][ans] = val
        code_q[r["item_code"]] = r["question"]
        slug = CR.SUM_INSTRUMENTS.get(r["instrument"])
        if slug:
            inst_items[slug].setdefault(norm_q(r["question"]), {})[ans] = val

    def score_items(item_list):
        """item_list: [(qid, ans_map, reverse, lo, hi)]. Yield {pid:(score,age)}."""
        n_items = len(item_list)
        need = math.ceil(CR.MIN_ITEM_FRACTION * n_items)
        pids = set()
        for qid, *_ in item_list:
            pids |= set(questions[qid]["responses"].keys())
        out = {}
        for pid in pids:
            got, age = [], None
            for qid, ans_map, rev, lo, hi in item_list:
                resp = questions[qid]["responses"].get(pid)
                if not resp:
                    continue
                answers = [a for a in resp[1] if not is_missing_answer(a)]
                if len(answers) != 1:
                    continue
                v = ans_map.get(answer_norm(answers[0]))
                if v is None:
                    continue
                if rev:
                    v = lo + hi - v
                got.append(v)
                age = resp[0]
            if len(got) >= need and age is not None and not math.isnan(age):
                out[pid] = (sum(got) / len(got) * n_items, age)  # prorated sum
        return out, n_items

    # merged-by-text sum instruments
    for slug, items in inst_items.items():
        item_list = []
        for nq, ans_map in items.items():
            qid = qtext_to_qid.get(nq)
            if qid is None or not ans_map:
                continue
            rev = any(frag in nq for frag in CR.REVERSE_TEXT_FRAGMENTS)
            item_list.append((qid, ans_map, rev, min(ans_map.values()), max(ans_map.values())))
        if len(item_list) < 2:
            continue
        vals, n_items = score_items(item_list)
        yield f"comp_{slug}", "composite", "quant", vals, {
            "question_concept_id": "", "question": slug,
            "answer": f"{n_items}-item prorated sum", "ordinal_rule": "composite",
            "covar_mode": "full",
        }

    # BFI-2-XS Big Five domains
    for dom, spec in CR.BFI2_DOMAINS.items():
        item_list = []
        for code, rev in spec:
            q, amap = code_q.get(code), code_answer.get(code)
            if not q or not amap:
                continue
            qid = qid_by_item.get(code) or qtext_to_qid.get(norm_q(q))
            if qid is None:
                continue
            item_list.append((qid, amap, rev, min(amap.values()), max(amap.values())))
        if len(item_list) < 2:
            continue
        vals, n_items = score_items(item_list)
        yield f"comp_{dom}", "composite", "quant", vals, {
            "question_concept_id": "", "question": dom,
            "answer": f"{n_items}-item BFI-2-XS domain", "ordinal_rule": "composite",
            "covar_mode": "full",
        }

    # explicit mixed-valence composites, scored from the ordinal manifest
    if ord_lookup and ordinal_manifest_path and Path(ordinal_manifest_path).exists():
        code_qtext, code_vals = {}, defaultdict(list)
        with open(ordinal_manifest_path, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                code_qtext[r["item_concept"]] = r["field_label"]
                try:
                    code_vals[r["item_concept"]].append(float(r["ordinal_value"]))
                except (ValueError, TypeError):
                    pass
        for slug, (desc, items) in CR.EXPLICIT_COMPOSITES.items():
            item_list = []
            for code, rev in items:
                qt = code_qtext.get(code)
                if not qt:
                    continue
                nq = norm_q(qt)
                qid = qid_by_item.get(code) or qtext_to_qid.get(nq)
                vals_c = code_vals.get(code)
                if qid is None or not vals_c:
                    continue
                ans_map = {a: v for (q_, a), v in ord_lookup.items() if q_ in (code, nq)}
                if not ans_map:
                    continue
                item_list.append((qid, ans_map, rev, min(vals_c), max(vals_c)))
            if len(item_list) < 2:
                continue
            vals, n_items = score_items(item_list)
            yield f"comp_{slug}", "composite", "quant", vals, {
                "question_concept_id": "", "question": slug,
                "answer": desc, "ordinal_rule": "composite", "covar_mode": "full",
            }


FITBIT_MIN_DAYS = 10


def _fitbit_person_means(csv_path, keep, cols):
    """Return streaming per-person sums/counts for retained Fitbit rows."""
    agg = defaultdict(lambda: defaultdict(float))
    if not csv_path or not Path(csv_path).exists() or Path(csv_path).stat().st_size == 0:
        return agg
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            pid = (row.get("person_id") or "").strip()
            if pid not in keep:
                continue
            for c in cols:
                try:
                    v = float(row[c])
                except (KeyError, ValueError, TypeError):
                    continue
                if not math.isfinite(v):
                    continue
                agg[pid][f"{c}_sum"] += v
                agg[pid][f"{c}_n"] += 1
            try:
                age = float(row["age"])
            except (KeyError, ValueError, TypeError):
                continue
            if math.isfinite(age):
                agg[pid]["age_sum"] += age
                agg[pid]["age_n"] += 1
    return agg


def build_fitbit_phenotypes(activity_csv, sleep_csv, keep):
    """Per-person Fitbit averages (>= FITBIT_MIN_DAYS valid days) as quant traits.

    Fitbit-worn behavioural/physiological traits: mean daily steps, sedentary and
    active minutes, sleep duration and efficiency. Residualized on the full
    covariate set (age at wear is known). Requires the AoU Fitbit tables; the
    orchestrator extracts them and skips gracefully if absent.
    """
    act = _fitbit_person_means(activity_csv, keep, ["steps", "sedentary_minutes", "active_minutes"])
    slp = _fitbit_person_means(sleep_csv, keep, ["minute_asleep", "sleep_efficiency"])

    specs = [
        ("fitbit_mean_daily_steps", act, "steps"),
        ("fitbit_sedentary_minutes", act, "sedentary_minutes"),
        ("fitbit_active_minutes", act, "active_minutes"),
        ("fitbit_sleep_minutes", slp, "minute_asleep"),
        ("fitbit_sleep_efficiency", slp, "sleep_efficiency"),
    ]
    for pheno_id, agg, col in specs:
        values = {}
        for pid, d in agg.items():
            n = int(d.get(f"{col}_n", 0))
            age_n = int(d.get("age_n", 0))
            if n < FITBIT_MIN_DAYS or age_n == 0:
                continue
            values[pid] = (d[f"{col}_sum"] / n, d["age_sum"] / age_n)
        yield pheno_id, "fitbit", "quant", values, {
            "question_concept_id": "", "question": pheno_id, "answer": f">= {FITBIT_MIN_DAYS} valid days",
            "ordinal_rule": "", "covar_mode": "full",
        }


def build_chronotype_phenotype(chrono_csv, keep):
    """Chronotype proxy = mean main-sleep onset clock hour (morningness-eveningness).

    Onset hours before noon are wrapped to [24,36) so a 1am onset (=25) sorts after
    a 10pm onset (=22); higher = later/more evening. >= FITBIT_MIN_DAYS nights.
    """
    if not chrono_csv or not Path(chrono_csv).exists() or Path(chrono_csv).stat().st_size == 0:
        return
    agg = defaultdict(lambda: defaultdict(float))
    with open(chrono_csv, newline="") as f:
        for row in csv.DictReader(f):
            pid = (row.get("person_id") or "").strip()
            if pid not in keep:
                continue
            try:
                h = float(row["onset_hour"])
                age = float(row["age"])
            except (KeyError, ValueError, TypeError):
                continue
            if not (math.isfinite(h) and math.isfinite(age)):
                continue
            agg[pid]["onset_sum"] += h + 24.0 if h < 12.0 else h
            agg[pid]["onset_n"] += 1
            agg[pid]["age_sum"] += age
            agg[pid]["age_n"] += 1
    values = {}
    for pid, d in agg.items():
        onset_n = int(d.get("onset_n", 0))
        age_n = int(d.get("age_n", 0))
        if onset_n < FITBIT_MIN_DAYS or age_n == 0:
            continue
        values[pid] = (d["onset_sum"] / onset_n, d["age_sum"] / age_n)
    yield "fitbit_chronotype_sleep_onset", "fitbit", "quant", values, {
        "question_concept_id": "", "question": "fitbit_chronotype_sleep_onset",
        "answer": "mean main-sleep onset clock hour; higher = later/evening", "ordinal_rule": "",
        "covar_mode": "full",
    }


def load_item_labels(inventory_path):
    """item_concept -> field label, from the survey item inventory."""
    out = {}
    if not inventory_path or not Path(inventory_path).exists():
        return out
    with open(inventory_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out.setdefault(r["item_concept"], r["field_label"])
    return out


def build_derived_psych_phenotypes(questions, item_labels, qid_by_item=None):
    """Algorithmic psychiatric phenotypes from the UKB-MHQ / CIDI-SF / PCL items.

    Screening-level derivations (documented in SPECSHEET 11d). All are sensitive
    (mental health / suicidality) and should carry the sensitive release tier.
    """
    qid_by_item = qid_by_item or {}
    qtext_to_qid = {}
    for qid, q in questions.items():
        qtext_to_qid.setdefault(norm_q(q["question"]), qid)

    def resp(code):
        lab = item_labels.get(code)
        qid = qid_by_item.get(code) or (qtext_to_qid.get(norm_q(lab)) if lab else None)
        return questions.get(qid, {}).get("responses", {}) if qid else {}

    def nonmiss(pid, r):
        v = r.get(pid)
        return [a for a in v[1] if not is_missing_answer(a)] if v else []

    def age_of(pid, *rs):
        for r in rs:
            v = r.get(pid)
            if v and v[0] is not None and not math.isnan(v[0]):
                return v[0]
        return None

    def yes(ans_list):
        return any(answer_norm(a) == "yes" or answer_norm(a).startswith("yes,") for a in ans_list)

    def binary_from(case_test, denom_codes, pheno_id, desc):
        """case_test(pid)->True/False/None(exclude); denom = answered any denom item."""
        rs = [resp(c) for c in denom_codes]
        pids = set().union(*[set(r) for r in rs]) if rs else set()
        values = {}
        for pid in pids:
            if not any(nonmiss(pid, r) for r in rs):
                continue
            y = case_test(pid)
            if y is None:
                continue
            a = age_of(pid, *rs)
            if a is not None:
                values[pid] = (1.0 if y else 0.0, a)
        return (pheno_id, "binary", "binary", values,
                {"question_concept_id": "", "question": pheno_id, "answer": desc,
                 "ordinal_rule": "derived_psych", "covar_mode": "full"})

    # ---- psychotic experiences (any of voices / thought-insertion / paranoia) --
    r21, r22, r23 = resp("cidi5_21"), resp("cidi5_22"), resp("cidi5_23")

    def psychosis(pid):
        a = nonmiss(pid, r21) + nonmiss(pid, r22) + nonmiss(pid, r23)
        return yes(a) if a else None

    yield binary_from(psychosis, ["cidi5_21", "cidi5_22", "cidi5_23"],
                      "psych_psychotic_experiences_any",
                      "any lifetime psychotic experience (voices/thought-insertion/paranoia)")

    # ---- suicidality / self-harm (each its own binary) ------------------------
    for code, pid_name, desc in [
        ("ss_1", "psych_self_harm_ideation_lifetime", "lifetime thoughts of purposely hurting yourself"),
        ("ss_2", "psych_suicidal_ideation_lifetime", "lifetime thoughts of killing yourself"),
        ("ss_3", "psych_suicide_attempt_lifetime", "lifetime suicide attempt"),
    ]:
        r = resp(code)
        yield binary_from(lambda pid, r=r: (yes(nonmiss(pid, r)) if nonmiss(pid, r) else None),
                          [code], pid_name, desc)

    # ---- mania screen / probable bipolar (UKB Smith 2013 style) ---------------
    r43, r44, r45, r46, r47 = (resp("mhqukb_43"), resp("mhqukb_44"), resp("mhqukb_45"),
                               resp("mhqukb_46"), resp("mhqukb_47"))

    def mania_core(pid):
        return yes(nonmiss(pid, r43)) or yes(nonmiss(pid, r44))

    def n_manic_symptoms(pid):
        return len(nonmiss(pid, r45))

    def mania_screen(pid):
        if not (nonmiss(pid, r43) or nonmiss(pid, r44)):
            return None
        return mania_core(pid) and n_manic_symptoms(pid) >= 3

    yield binary_from(mania_screen, ["mhqukb_43", "mhqukb_44"], "psych_mania_episode_screen",
                      "manic/irritable episode + >=3 symptoms (screen)")

    def probable_bipolar(pid):
        base = mania_screen(pid)
        if base is None:
            return None
        if not base:
            return False
        dur = {answer_norm(a) for a in nonmiss(pid, r46)}
        long_enough = any("four days" in d or "week or more" in d for d in dur)
        prob = any("needed treatment or caused problems" in answer_norm(a) for a in nonmiss(pid, r47))
        return bool(long_enough and prob)

    yield binary_from(probable_bipolar, ["mhqukb_43", "mhqukb_44"], "psych_probable_bipolar",
                      "probable bipolar: mania screen + >=4 day duration + impairment (UKB-style)")

    # ---- depression: lifetime episode + probable recurrent --------------------
    r5, r6, r24 = resp("mhqukb_5"), resp("mhqukb_6"), resp("mhqukb_24")

    def depressed_episode(pid):
        a = nonmiss(pid, r5) + nonmiss(pid, r6)
        return yes(a) if a else None

    yield binary_from(depressed_episode, ["mhqukb_5", "mhqukb_6"], "psych_lifetime_depressed_episode",
                      "ever a >=2-week period of low mood / anhedonia (screen)")

    def recurrent_dep(pid):
        core = depressed_episode(pid)
        if core is None:
            return None
        if not core:
            return False
        return any("several" in answer_norm(a) for a in nonmiss(pid, r24))

    yield binary_from(recurrent_dep, ["mhqukb_5", "mhqukb_6"], "psych_probable_recurrent_depression",
                      "lifetime depressed episode + several episodes (recurrent-MDD proxy)")


def build_acculturation_phenotype(questions, item_labels, qid_by_item=None):
    """Cultural-assimilation index: US-born + English-at-home + English proficiency.

    Higher = more acculturated. Components normalised to 0-1 and summed; English
    proficiency is imputed maximal for participants who speak English at home.
    """
    qid_by_item = qid_by_item or {}
    qtext_to_qid = {}
    for qid, q in questions.items():
        qtext_to_qid.setdefault(norm_q(q["question"]), qid)

    def resp(code):
        lab = item_labels.get(code)
        qid = qid_by_item.get(code) or (qtext_to_qid.get(norm_q(lab)) if lab else None)
        return questions.get(qid, {}).get("responses", {}) if qid else {}

    born, lang, prof = resp("thebasics_birthplace"), resp("chis_1"), resp("chis_1_xx")
    prof_map = {"not at all": 0.0, "not well": 1.0 / 3, "well": 2.0 / 3, "very well": 1.0}

    def one(pid, r):
        v = r.get(pid)
        nm = [a for a in v[1] if not is_missing_answer(a)] if v else []
        return nm[0] if len(nm) == 1 else None

    pids = set(born) | set(lang) | set(prof)
    values = {}
    for pid in pids:
        b = one(pid, born)
        if b is None:
            continue  # birthplace is the anchor component
        us_born = 1.0 if answer_norm(b) == "usa" else 0.0
        lg = one(pid, lang)
        english_home = 1.0 if (lg and answer_norm(lg) == "no") else (0.0 if lg else None)
        if english_home == 1.0:
            eng = 1.0
        else:
            pv = one(pid, prof)
            eng = prof_map.get(answer_norm(pv)) if pv else None
        comps = [us_born] + ([english_home] if english_home is not None else []) + ([eng] if eng is not None else [])
        age = None
        for r in (born, lang, prof):
            v = r.get(pid)
            if v and v[0] is not None and not math.isnan(v[0]):
                age = v[0]
                break
        if len(comps) >= 1 and age is not None:
            values[pid] = (sum(comps), age)  # 0..3 acculturation index
    yield ("accult_index", "acculturation", "quant", values,
           {"question_concept_id": "", "question": "acculturation_index",
            "answer": "US-born + English-at-home + English proficiency", "ordinal_rule": "",
            "covar_mode": "full"})


def build_state_cluster_phenotypes(questions, item_labels, clusters_path, state_csv, keep, qid_by_item=None):
    """One binary GWAS per state cluster: membership (state in cluster) vs not.

    State source: a person->state CSV via `state_csv` (person_id,state[,age]) if
    given (e.g. residence state from ZIP3); otherwise the survey work-address
    state item. NOTE these phenotypes largely capture residual geographic genetic
    structure (PCs absorb much of it) -- interpret as geography/structure, not trait
    biology.
    """
    if not clusters_path or not Path(clusters_path).exists():
        return
    clusters = {}
    qid_by_item = qid_by_item or {}
    with open(clusters_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            clusters[r["cluster"]] = {s.strip().lower() for s in r["states"].split("|")}

    state_by_person = {}  # pid -> (state_lower, age)
    if state_csv and Path(state_csv).exists() and Path(state_csv).stat().st_size > 0:
        with open(state_csv, newline="") as f:
            for row in csv.DictReader(f):
                pid = (row.get("person_id") or "").strip()
                st = (row.get("state") or "").strip()
                if pid in keep and st:
                    try:
                        age = float(row.get("age", "nan"))
                    except (ValueError, TypeError):
                        age = float("nan")
                    state_by_person[pid] = (st.lower(), age)
        src = f"residence-state CSV ({Path(state_csv).name})"
    else:
        qtext_to_qid = {norm_q(q["question"]): qid for qid, q in questions.items()}
        lab = item_labels.get("employmentworkaddress_state")
        qid = qid_by_item.get("employmentworkaddress_state") or (qtext_to_qid.get(norm_q(lab)) if lab else None)
        if qid:
            for pid, (age, answers) in questions[qid]["responses"].items():
                nm = [a for a in answers if not is_missing_answer(a)]
                if len(nm) == 1:
                    state_by_person[pid] = (answer_norm(nm[0]), age)
        src = "survey work-address state"
    if not state_by_person:
        return
    log(f"  state-cluster source: {src} ({len(state_by_person)} with a state)")

    for cluster, states in clusters.items():
        values = {}
        for pid, (st, age) in state_by_person.items():
            if age is None or math.isnan(age):
                continue
            values[pid] = (1.0 if st in states else 0.0, age)
        yield f"geo_{cluster}", "geographic", "binary", values, {
            "question_concept_id": "", "question": f"state_cluster_{cluster}",
            "answer": "member vs non-member", "ordinal_rule": "", "covar_mode": "full",
        }


def build_external_score_phenotypes(config_path, keep):
    """Yield pre-computed continuous scores (ETM cognitive tasks + EA/SES proxies).

    These come from the ea_proxy workflow (ETM task scoring + the SES-EA/g-EA
    proxy models). They are already age/sex-normalized upstream, so they are
    residualized on sex_c + PC1..PC10 only (covar_mode=sexpc), matching the
    repo's final g-EA proxy GWAS covariate choice. The registry lives in
    metadata/external_scores.tsv; missing score files are skipped with a note.
    """
    if not config_path or not Path(config_path).exists():
        return
    with open(config_path, newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if not r["phenotype_id"].startswith("#")]
    for r in rows:
        path = Path(os.path.expandvars(r["score_file"].strip()))
        if not path.exists() or path.stat().st_size == 0:
            log(f"  external score missing, skipping {r['phenotype_id']}: {path}")
            continue
        iid_col, val_col = r["iid_col"].strip(), r["value_col"].strip()
        delim = "\t" if path.suffix in (".tsv", ".txt") else ","
        vals = {}
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter=delim):
                iid = (row.get(iid_col) or "").strip()
                if iid not in keep:
                    continue
                try:
                    vals[iid] = (float(row[val_col]), float("nan"))
                except (KeyError, ValueError, TypeError):
                    continue
        if not vals:
            log(f"  WARN: {r['phenotype_id']} matched 0 samples in {path} "
                f"(check iid_col={iid_col!r}, value_col={val_col!r}).")
            continue
        yield f"cog_{r['phenotype_id']}", "external_score", "quant", vals, {
            "question_concept_id": "",
            "question": r.get("description", ""),
            "answer": "",
            "ordinal_rule": "",
            "covar_mode": "sexpc",
        }


def external_score_pheno_ids(config_path: Path | None) -> set[str]:
    ids: set[str] = set()
    if not config_path or not Path(config_path).exists():
        return ids
    with open(config_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("phenotype_id") or "").strip()
            if pid and not pid.startswith("#"):
                ids.add(f"cog_{pid}")
    return ids


def build_zip3_ses_phenotypes(zip3_ses_csv: Path | None, keep: set[str]):
    """Yield ZIP3-level socioeconomic context phenotypes.

    The orchestrator extracts the latest row per person from AoU
    ds_zip_code_socioeconomic. Raw ZIP3 and ACS vintage are retained only in the
    extract for auditability; the GWAS phenotypes are the seven numeric SES
    fields. These are contextual/geographic traits, not individual survey
    responses.
    """
    if not zip3_ses_csv or not Path(zip3_ses_csv).exists() or Path(zip3_ses_csv).stat().st_size == 0:
        return
    vals = {pid: {} for pid in ZIP3_SES_TRAITS}
    with open(zip3_ses_csv, newline="") as f:
        for row in csv.DictReader(f):
            iid = (row.get("person_id") or row.get("IID") or "").strip()
            if iid not in keep:
                continue
            try:
                age = float(row.get("age_at_observation") or row.get("age") or "nan")
            except (TypeError, ValueError):
                age = float("nan")
            for pid, (col, _desc) in ZIP3_SES_TRAITS.items():
                try:
                    y = float(row[col])
                except (KeyError, TypeError, ValueError):
                    continue
                if math.isnan(y):
                    continue
                vals[pid][iid] = (y, age)
    for pid, (_col, desc) in ZIP3_SES_TRAITS.items():
        if not vals[pid]:
            continue
        yield pid, "zip3_ses", "quant", vals[pid], {
            "question_concept_id": "",
            "question": desc,
            "answer": "",
            "ordinal_rule": "",
            "covar_mode": "full",
        }


def wants_phenotype_source(only: set[str], prefixes=(), exact=()) -> bool:
    if not only:
        return True
    exact_set = set(exact)
    if exact_set & only:
        return True
    return any(any(pid.startswith(prefix) for prefix in prefixes) for pid in only)


def build_measurement_phenotypes(meas_csv: Path, keep: set[str]):
    if not meas_csv or not meas_csv.exists() or meas_csv.stat().st_size == 0:
        return
    # per person: earliest (closest to baseline) valid value per phenotype
    height = defaultdict(list)
    weight = defaultdict(list)
    per = {name: defaultdict(list) for name in MEASUREMENTS}
    with open(meas_csv, newline="") as f:
        for row in csv.DictReader(f):
            pid = (row.get("person_id") or "").strip()
            if pid not in keep:
                continue
            try:
                cid = int(float(row.get("measurement_concept_id") or "nan"))
                val = float(row.get("value_as_number") or "nan")
                age = float(row.get("age_at_measurement") or "nan")
            except ValueError:
                continue
            dt = (row.get("measurement_datetime") or "").strip()
            if math.isnan(val):
                continue
            if cid == 903133:
                height[pid].append((dt, val, age))
            elif cid == 903121:
                weight[pid].append((dt, val, age))
            for name, (cids, _u, _r) in MEASUREMENTS.items():
                if cid in cids:
                    per[name][pid].append((dt, val, age))

    def earliest(recs, lo, hi):
        recs = [(d, v, a) for d, v, a in recs if lo <= v <= hi]
        if not recs:
            return None
        recs.sort(key=lambda t: t[0])
        return recs[0][1], recs[0][2]

    for name, (_cids, _u, (lo, hi)) in MEASUREMENTS.items():
        values = {}
        for pid, recs in per[name].items():
            e = earliest(recs, lo, hi)
            if e:
                values[pid] = (e[0], e[1])
        yield name, "measurement", "quant", values, {"question_concept_id": "", "question": name, "answer": "", "ordinal_rule": ""}

    # BMI = weight_kg / (height_m^2), from earliest same-era height+weight
    bmi = {}
    for pid in set(height) & set(weight):
        h = earliest(height[pid], 100.0, 250.0)
        w = earliest(weight[pid], 20.0, 400.0)
        if h and w:
            b = w[0] / ((h[0] / 100.0) ** 2)
            if 12.0 <= b <= 80.0:
                bmi[pid] = (b, w[1])
    yield "bmi_kg_m2", "measurement", "quant", bmi, {"question_concept_id": "", "question": "bmi_kg_m2", "answer": "", "ordinal_rule": ""}

    # Pulse pressure and MAP require a same-date SBP+DBP pair (spec 10.3).
    sbp = per["systolic_bp_mmhg"]
    dbp = per["diastolic_bp_mmhg"]
    pp, mapv = {}, {}
    for pid in set(sbp) & set(dbp):
        dmap = {}
        for d, v, a in sbp[pid]:
            if 70.0 <= v <= 260.0:
                dmap.setdefault(d, {})["s"] = (v, a)
        for d, v, a in dbp[pid]:
            if 40.0 <= v <= 160.0 and d in dmap:
                dmap[d]["d"] = (v, a)
        pairs = [
            (d, val["s"][0], val["d"][0], val["s"][1])
            for d, val in dmap.items()
            if "s" in val and "d" in val and val["s"][0] > val["d"][0]
        ]
        if not pairs:
            continue
        pairs.sort(key=lambda t: t[0])  # earliest visit -> closest to baseline
        _, s, dd, age = pairs[0]
        pulse = s - dd
        m = dd + pulse / 3.0
        if 15.0 <= pulse <= 150.0:
            pp[pid] = (pulse, age)
        if 50.0 <= m <= 180.0:
            mapv[pid] = (m, age)
    yield "pulse_pressure_mmhg", "measurement", "quant", pp, {"question_concept_id": "", "question": "pulse_pressure_mmhg", "answer": "", "ordinal_rule": ""}
    yield "mean_arterial_pressure_mmhg", "measurement", "quant", mapv, {"question_concept_id": "", "question": "mean_arterial_pressure_mmhg", "answer": "", "ordinal_rule": ""}


# --------------------------------------------------------------------------- #
# residualize + write + PLINK2
# --------------------------------------------------------------------------- #
def prepare_and_write(pheno_id, kind, values, sex, pcs, fid_by_iid, outdir, covar_mode="full"):
    """Return a prep dict, or a dict with skip_reason if the phenotype fails QC.

    covar_mode: "full"  -> age_c, sex_c, age_c:sex_c, PC1..PC10 (survey/measurement)
                "sexpc" -> sex_c, PC1..PC10 only (pre-age-normalized external scores,
                           matching the repo's final g-EA proxy GWAS covariates).
    """
    need_age = covar_mode == "full"
    rows = []
    for iid, (y, age) in values.items():
        if iid not in sex or iid not in pcs or math.isnan(y):
            continue
        if need_age and (age is None or math.isnan(age)):
            continue
        rows.append((iid, y, age, sex[iid], pcs[iid]))

    if kind == "binary":
        ncase = sum(1 for _, y, *_ in rows if y == 1.0)
        nctrl = sum(1 for _, y, *_ in rows if y == 0.0)
        if ncase < MIN_CASES or nctrl < MIN_CONTROLS:
            return {
                "skip_reason": "too_few_cases_or_controls",
                "n": len(rows),
                "n_cases": ncase,
                "n_controls": nctrl,
            }
    else:
        if len(rows) < MIN_QUANT_N:
            return {
                "skip_reason": "too_few_nonmissing",
                "n": len(rows),
                "n_cases": 0,
                "n_controls": 0,
            }
        levels = len({round(y, 6) for _, y, *_ in rows})
        if levels < MIN_ORDINAL_LEVELS:
            return {
                "skip_reason": "too_few_observed_levels",
                "n": len(rows),
                "n_cases": 0,
                "n_controls": 0,
                "n_levels": levels,
            }
        ncase = nctrl = 0

    iids = [r[0] for r in rows]
    y = np.array([r[1] for r in rows], dtype=float)
    sex_c = np.array([r[3] for r in rows], dtype=float) - 0.5
    pc = np.array([r[4] for r in rows], dtype=float)
    if covar_mode == "sexpc":
        covars = np.column_stack([sex_c, pc])
    else:
        age = np.array([r[2] for r in rows], dtype=float)
        age_c = age - age.mean()
        covars = np.column_stack([age_c, sex_c, age_c * sex_c, pc])

    if kind == "binary":
        pheno_vec = residualize(y, covars)          # LPM residual, no INT
    else:
        pheno_vec = residualize(rint(y), covars)    # INT then residualize

    pdir = outdir / "phenotypes"
    pdir.mkdir(parents=True, exist_ok=True)
    name = f"{pheno_id}_resid"
    raw_path = pdir / f"{pheno_id}.raw.pheno.tsv"
    pheno_path = pdir / f"{pheno_id}.resid.pheno.tsv"
    keep_path = pdir / f"{pheno_id}.keep.tsv"
    with open(raw_path, "w") as f:
        f.write(f"FID\tIID\t{pheno_id}_raw\n")
        for iid, v in zip(iids, y):
            f.write(f"{fid_by_iid.get(iid, iid)}\t{iid}\t{v:.17g}\n")
    with open(pheno_path, "w") as f:
        f.write(f"FID\tIID\t{name}\n")
        for iid, v in zip(iids, pheno_vec):
            f.write(f"{fid_by_iid.get(iid, iid)}\t{iid}\t{v:.17g}\n")
    with open(keep_path, "w") as f:
        for iid in iids:
            f.write(f"{fid_by_iid.get(iid, iid)}\t{iid}\n")
    return {
        "raw_path": raw_path,
        "pheno_path": pheno_path,
        "keep_path": keep_path,
        "pheno_name": name,
        "n": len(rows),
        "n_cases": ncase,
        "n_controls": nctrl,
    }


def final_gwas_paths(outdir: Path, pheno_id: str, pheno_name: str) -> tuple[Path, Path, Path]:
    pdir = outdir / "gwas" / pheno_id
    prefix = pdir / pheno_id
    glm = pdir / f"{pheno_id}.{pheno_name}.glm.linear"
    lite = pdir / f"{pheno_id}.sumstats.tsv.gz"
    return prefix, glm, lite


def read_pheno_values(path: Path, pheno_name: str) -> dict[tuple[str, str], str]:
    vals = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            vals[(row["FID"], row["IID"])] = row[pheno_name]
    return vals


def write_batch_pheno(batch_jobs, pheno_path: Path, keep_path: Path) -> None:
    pheno_path.parent.mkdir(parents=True, exist_ok=True)
    by_name = []
    sample_keys = set()
    for job in batch_jobs:
        vals = read_pheno_values(job["pheno_path"], job["pheno_name"])
        by_name.append((job["pheno_name"], vals))
        sample_keys.update(vals.keys())
    ordered = sorted(sample_keys, key=lambda x: (0, int(x[1])) if x[1].isdigit() else (1, x[1]))
    with open(pheno_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["FID", "IID", *[name for name, _ in by_name]])
        for key in ordered:
            writer.writerow([key[0], key[1], *[vals.get(key, "NA") for _, vals in by_name]])
    with open(keep_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for fid, iid in ordered:
            writer.writerow([fid, iid])


def run_plink2_batch(plink2, bfile, batch_jobs, workdir: Path, batch_index: int):
    workdir.mkdir(parents=True, exist_ok=True)
    batch_prefix = workdir / f"batch_{batch_index:05d}"
    batch_pheno = workdir / f"batch_{batch_index:05d}.pheno.tsv"
    batch_keep = workdir / f"batch_{batch_index:05d}.keep.tsv"
    batch_log = workdir / f"batch_{batch_index:05d}.plink2.log"
    for old in workdir.glob(f"batch_{batch_index:05d}*.glm.linear"):
        old.unlink()
    for old in (workdir / f"batch_{batch_index:05d}.log", batch_log):
        if old.exists():
            old.unlink()
    write_batch_pheno(batch_jobs, batch_pheno, batch_keep)
    cmd = [
        plink2,
        "--bfile", str(bfile),
        "--keep", str(batch_keep),
        "--pheno", str(batch_pheno),
        "--pheno-name", ",".join(job["pheno_name"] for job in batch_jobs),
        "--glm", "allow-no-covars", "cols=chrom,pos,a1freq,nobs,beta,se,p",
        "--no-input-missing-phenotype",
        "--out", str(batch_prefix),
    ]
    start = time.time()
    res = subprocess.run(cmd, text=True, capture_output=True)
    elapsed = time.time() - start
    batch_log.write_text(res.stderr + "\n" + res.stdout)
    if res.returncode != 0:
        raise RuntimeError(f"PLINK2 failed for batch {batch_index}; see {batch_log}.")

    outputs = {}
    for job in batch_jobs:
        local_glm = workdir / f"batch_{batch_index:05d}.{job['pheno_name']}.glm.linear"
        if not local_glm.exists() or local_glm.stat().st_size == 0:
            raise RuntimeError(f"PLINK2 did not write expected output: {local_glm}")
        final_glm = job["glm"]
        final_lite = job["sumstats"]
        final_glm.parent.mkdir(parents=True, exist_ok=True)
        write_lightweight_sumstats(local_glm, final_lite)
        shutil.copy2(local_glm, final_glm)
        shutil.copy2(batch_log, final_glm.parent / f"{job['pheno_id']}.plink2.log")
        outputs[job["pheno_id"]] = (final_glm, final_lite, elapsed)
    return outputs


def write_lightweight_sumstats(glm_path: Path, out_path: Path) -> Path:
    """Write compact association columns for downstream scans."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(glm_path, newline="") as src, gzip.open(out_path, "wt", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fields = ["rsid", "chrom", "pos", "allele1", "a1freq", "n", "beta", "se", "p", "log10p"]
        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            p_raw = row.get("P", "")
            log10p = ""
            try:
                p = float(p_raw)
                if p > 0:
                    log10p = f"{-math.log10(p):.8g}"
            except (TypeError, ValueError):
                pass
            writer.writerow({
                "rsid": row.get("ID", ""),
                "chrom": row.get("#CHROM", row.get("CHROM", "")),
                "pos": row.get("POS", ""),
                "allele1": row.get("A1", ""),
                "a1freq": row.get("A1_FREQ", row.get("A1FREQ", "")),
                "n": row.get("OBS_CT", ""),
                "beta": row.get("BETA", ""),
                "se": row.get("SE", ""),
                "p": p_raw,
                "log10p": log10p,
            })
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bfile", type=Path, required=True)
    ap.add_argument("--keep", type=Path, required=True)
    ap.add_argument("--sex", type=Path, required=True)
    ap.add_argument("--pcs", type=Path, required=True)
    ap.add_argument("--survey-csv", type=Path, required=True)
    ap.add_argument("--bhp-csv", type=Path, default=None)
    ap.add_argument("--measurements-csv", type=Path, default=None)
    ap.add_argument("--zip3-ses-csv", type=Path, default=None,
                    help="latest per-person AoU ds_zip_code_socioeconomic extract.")
    ap.add_argument("--fitbit-activity-csv", type=Path, default=None)
    ap.add_argument("--fitbit-sleep-csv", type=Path, default=None)
    ap.add_argument("--fitbit-chronotype-csv", type=Path, default=None)
    ap.add_argument("--question-manifest", type=Path, required=True)
    ap.add_argument("--aou-question-concepts", type=Path, default=None,
                    help="AoU live ds_survey question_concept_id metadata crosswalk.")
    ap.add_argument("--ea-proxy-feature-manifest", type=Path, default=None,
                    help="Supplemental SES-EA XGBoost source-question manifest.")
    ap.add_argument("--ordinal-manifest", type=Path, required=True)
    ap.add_argument("--item-inventory", type=Path, default=None,
                    help="survey_item_inventory.tsv (for derived-psych/acculturation item labels).")
    ap.add_argument("--state-clusters", type=Path, default=None,
                    help="state_clusters.tsv (geographic/political cluster membership GWAS).")
    ap.add_argument("--state-csv", type=Path, default=None,
                    help="optional person_id,state,age CSV (residence state); else survey work state.")
    ap.add_argument("--pfhh-allowlist", type=Path, default=None)
    ap.add_argument("--composite-manifest", type=Path, default=None,
                    help="composite_items_manifest.tsv (validated sum/domain scores).")
    ap.add_argument("--external-scores", type=Path, default=None,
                    help="registry TSV of pre-computed cognitive/EA-proxy scores to GWAS.")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--gwas-workdir", type=Path, default=None,
                    help="Local working directory for temporary batched PLINK2 output.")
    ap.add_argument("--gwas-batch-size", type=int, default=64,
                    help="Number of residualized phenotypes per PLINK2 genotype scan.")
    ap.add_argument("--phenotypes", default="", help="comma-separated pheno_id filter (smoke test).")
    ap.add_argument("--plink2-bin", default=shutil.which("plink2") or "plink2")
    ap.add_argument("--skip-gwas", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.gwas_batch_size < 1:
        raise SystemExit("--gwas-batch-size must be >= 1")
    if args.gwas_workdir is None:
        args.gwas_workdir = args.outdir / "work" / "gwas"

    keep = set(load_keep(args.keep))
    sex = load_sex(args.sex)
    pcs = load_pcs(args.pcs)
    fid_by_iid = load_fam_fids(Path(f"{args.bfile}.fam"))
    log(f"keep={len(keep)}  sex={len(sex)}  pcs={len(pcs)}  fam={len(fid_by_iid)}")

    qman, qman_by_item, qman_rows = load_question_manifest(args.question_manifest)
    qman_by_qid, qid_by_item = load_live_question_crosswalk(
        args.aou_question_concepts, qman_rows, qman
    )
    qman.update(qman_by_qid)
    ea_proxy_rows = load_ea_proxy_feature_sources(args.ea_proxy_feature_manifest, qman)
    qman_rows.extend(ea_proxy_rows)
    for row in ea_proxy_rows:
        qid = (row.get("question_concept_id") or "").strip()
        if qid:
            qman.setdefault(qid, row)
        label = (row.get("field_label") or "").strip()
        if label:
            qman.setdefault(norm_q(label), row)
    live_override_rows = load_live_question_overrides(args.aou_question_concepts, qman)
    qman_rows.extend(live_override_rows)
    for row in live_override_rows:
        qid = (row.get("question_concept_id") or "").strip()
        if qid:
            qman.setdefault(qid, row)
        label = (row.get("field_label") or "").strip()
        if label:
            qman.setdefault(norm_q(label), row)
    ord_lookup = load_ordinal_lookup(args.ordinal_manifest)
    item_labels = load_item_labels(args.item_inventory)
    log(
        f"manifest questions={len(qman_rows)}  live qid links={len(qman_by_qid)}  "
        f"ea_proxy_supplemental={len(ea_proxy_rows)}  "
        f"live_qid_overrides={len(live_override_rows)}  "
        f"ordinal answer maps={len(ord_lookup)}  item labels={len(item_labels)}"
    )

    survey_paths = [args.survey_csv]
    if args.bhp_csv:
        survey_paths.append(args.bhp_csv)
    allowed_qids = set(qman_by_qid)
    allowed_qids.update(qid for qid in qid_by_item.values() if qid)
    allowed_qids.update(
        (row.get("question_concept_id") or "").strip()
        for row in qman_rows
        if (row.get("question_concept_id") or "").strip()
    )
    allowed_qids.update(load_tsv_column_values(args.pfhh_allowlist, "question_concept_id"))
    allowed_qids.update(PFHH_SCREEN_QID.values())
    allowed_question_texts = {
        norm_q(row.get("field_label") or "")
        for row in qman_rows
        if (row.get("field_label") or "").strip()
        and not (row.get("disposition") or "").startswith("excluded")
    }
    log("Building latest-response table ...")
    log(
        f"survey row filter: qids={len(allowed_qids)} "
        f"question_texts={len(allowed_question_texts)}"
    )
    questions = build_latest_responses(survey_paths, keep, allowed_qids, allowed_question_texts)
    log(f"questions with responses: {len(questions)}")

    only = {p.strip() for p in args.phenotypes.split(",") if p.strip()}

    builders = []
    if wants_phenotype_source(only, prefixes=("bin_", "ord_")):
        builders.append(build_survey_phenotypes(questions, qman, ord_lookup))
    if wants_phenotype_source(only, prefixes=("num_",)):
        builders.append(build_numeric_phenotypes(questions, qman))
    if wants_phenotype_source(only, prefixes=("pfhh_",)):
        builders.append(build_pfhh_phenotypes(questions, args.pfhh_allowlist))
    if wants_phenotype_source(only, prefixes=("comp_",)):
        builders.append(build_composite_phenotypes(
            questions, args.composite_manifest, args.ordinal_manifest, ord_lookup, qid_by_item
        ))
    if wants_phenotype_source(only, prefixes=("psych_",)):
        builders.append(build_derived_psych_phenotypes(questions, item_labels, qid_by_item))
    if wants_phenotype_source(only, exact=("accult_index",)):
        builders.append(build_acculturation_phenotype(questions, item_labels, qid_by_item))
    if wants_phenotype_source(only, prefixes=("geo_",)):
        builders.append(build_state_cluster_phenotypes(
            questions, item_labels, args.state_clusters, args.state_csv, keep, qid_by_item
        ))
    measurement_ids = set(MEASUREMENTS) | {"bmi_kg_m2", "pulse_pressure_mmhg", "mean_arterial_pressure_mmhg"}
    if wants_phenotype_source(only, exact=measurement_ids):
        builders.append(build_measurement_phenotypes(args.measurements_csv, keep))
    if wants_phenotype_source(only, exact=set(ZIP3_SES_TRAITS)):
        builders.append(build_zip3_ses_phenotypes(args.zip3_ses_csv, keep))
    fitbit_ids = {
        "fitbit_mean_daily_steps",
        "fitbit_sedentary_minutes",
        "fitbit_active_minutes",
        "fitbit_sleep_minutes",
        "fitbit_sleep_efficiency",
    }
    if wants_phenotype_source(only, exact=fitbit_ids):
        builders.append(build_fitbit_phenotypes(args.fitbit_activity_csv, args.fitbit_sleep_csv, keep))
    if wants_phenotype_source(only, exact=("fitbit_chronotype_sleep_onset",)):
        builders.append(build_chronotype_phenotype(args.fitbit_chronotype_csv, keep))
    external_ids = external_score_pheno_ids(args.external_scores)
    if wants_phenotype_source(only, exact=external_ids):
        builders.append(build_external_score_phenotypes(args.external_scores, keep))

    manifest_rows = []
    skipped_rows = []
    gwas_jobs = []
    seen_pheno_ids = set()
    metadir = args.outdir / "metadata"
    metadir.mkdir(parents=True, exist_ok=True)

    for gen in builders:
        for pheno_id, trait_type, kind, values, meta in gen:
            if only and pheno_id not in only:
                continue
            if pheno_id in seen_pheno_ids:
                skipped_rows.append({
                    "pheno_id": pheno_id,
                    "trait_type": trait_type,
                    "kind": kind,
                    "skip_reason": "duplicate_pheno_id",
                    "n": len(values),
                    "n_cases": "",
                    "n_controls": "",
                    "n_levels": "",
                    "question_concept_id": meta.get("question_concept_id", ""),
                    "item_concept": meta.get("item_concept", ""),
                    "question": meta.get("question", ""),
                    "answer": meta.get("answer", ""),
                })
                continue
            prep = prepare_and_write(
                pheno_id, kind, values, sex, pcs, fid_by_iid, args.outdir, meta.get("covar_mode", "full")
            )
            if "skip_reason" in prep:
                skipped_rows.append({
                    "pheno_id": pheno_id,
                    "trait_type": trait_type,
                    "kind": kind,
                    "skip_reason": prep["skip_reason"],
                    "n": prep.get("n", 0),
                    "n_cases": prep.get("n_cases", 0),
                    "n_controls": prep.get("n_controls", 0),
                    "n_levels": prep.get("n_levels", ""),
                    "question_concept_id": meta.get("question_concept_id", ""),
                    "item_concept": meta.get("item_concept", ""),
                    "question": meta.get("question", ""),
                    "answer": meta.get("answer", ""),
                })
                continue
            seen_pheno_ids.add(pheno_id)
            row = {
                "pheno_id": pheno_id,
                "trait_type": trait_type,
                "kind": kind,
                "n": prep["n"],
                "n_cases": prep["n_cases"],
                "n_controls": prep["n_controls"],
                "pheno_name": prep["pheno_name"],
                "ordinal_rule": meta.get("ordinal_rule", ""),
                "question_concept_id": meta.get("question_concept_id", ""),
                "item_concept": meta.get("item_concept", ""),
                "question": meta.get("question", ""),
                "answer": meta.get("answer", ""),
                "raw_pheno_path": str(prep["raw_path"]),
                "pheno_path": str(prep["pheno_path"]),
            }
            _, glm, lite = final_gwas_paths(args.outdir, pheno_id, prep["pheno_name"])
            row["glm"] = str(glm)
            row["sumstats"] = str(lite)
            if not args.skip_gwas:
                if glm.exists() and glm.stat().st_size > 0 and lite.exists() and lite.stat().st_size > 0 and not args.force:
                    row["gwas_seconds"] = 0.0
                else:
                    row["gwas_seconds"] = ""
                    gwas_jobs.append({
                        "pheno_id": pheno_id,
                        "pheno_name": prep["pheno_name"],
                        "pheno_path": prep["pheno_path"],
                        "glm": glm,
                        "sumstats": lite,
                        "row": row,
                    })
            manifest_rows.append(row)
            if len(manifest_rows) % 100 == 0:
                log(f"  {len(manifest_rows)} phenotypes done")

    if gwas_jobs:
        log("Releasing phenotype-construction state before GWAS ...")
        builders = None
        questions = None
        qman = None
        qman_by_item = None
        qman_by_qid = None
        qman_rows = None
        qid_by_item = None
        ord_lookup = None
        item_labels = None
        allowed_qids = None
        allowed_question_texts = None
        sex = None
        pcs = None
        fid_by_iid = None
        keep = None
        gc.collect()
        log(
            f"Running {len(gwas_jobs)} GWAS phenotypes in batches of "
            f"{args.gwas_batch_size} under {args.gwas_workdir} ..."
        )
        for i in range(0, len(gwas_jobs), args.gwas_batch_size):
            batch_index = i // args.gwas_batch_size + 1
            batch = gwas_jobs[i:i + args.gwas_batch_size]
            log(f"  PLINK2 batch {batch_index}: {len(batch)} phenotypes")
            outputs = run_plink2_batch(args.plink2_bin, args.bfile, batch, args.gwas_workdir, batch_index)
            for job in batch:
                _, _, elapsed = outputs[job["pheno_id"]]
                job["row"]["gwas_seconds"] = round(elapsed, 1)

    man_path = metadir / "phenotype_manifest.tsv"
    if manifest_rows:
        cols = list(manifest_rows[0].keys())
        with open(man_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(manifest_rows)
    skip_path = metadir / "skipped_phenotypes.tsv"
    if skipped_rows:
        cols = list(skipped_rows[0].keys())
        with open(skip_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(skipped_rows)
    else:
        skip_path.write_text("pheno_id\tskip_reason\n")
    (metadir / "run_manifest.json").write_text(
        json.dumps(
            {
                "created_at": datetime.now().isoformat(),
                "bfile": str(args.bfile),
                "keep": str(args.keep),
                "covariates": ["age_c", "sex_c", "age_c_sex_c_inter", *PC_COLUMNS],
                "min_cases": MIN_CASES,
                "min_controls": MIN_CONTROLS,
                "min_quant_n": MIN_QUANT_N,
                "n_phenotypes_passing_qc": len(manifest_rows),
                "n_phenotypes_skipped_qc": len(skipped_rows),
                "skip_gwas": bool(args.skip_gwas),
            },
            indent=2,
        )
        + "\n"
    )
    log(f"Wrote {man_path} ({len(manifest_rows)} phenotypes passed QC; {len(skipped_rows)} skipped)")


if __name__ == "__main__":
    main()
