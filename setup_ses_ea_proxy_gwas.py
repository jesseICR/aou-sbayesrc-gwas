#!/usr/bin/env python3
"""Build AoU SES-EA proxy scores and GWAS input files.

This helper consumes workspace-local extracts produced by
setup_ses_ea_proxy_gwas.sh. It trains the primary-only SES/behavior proxy with
5-fold out-of-fold scoring over fit_pca_iids, then trains a sixth model on
kinship-clean eligible fit_pca_iids and applies it to the remaining eligible
classified-EUR samples.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import statistics
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd


EA_MAPPING = {
    1585941: 9.0,
    1585942: 9.0,
    1585943: 9.0,
    1585944: 10.0,
    1585945: 13.0,
    1585946: 15.0,
    1585947: 18.0,
    1585948: 20.0,
}

PMI_SKIP_ANSWER_ID = 903096
PMI_PREFER_NOT_ANSWER_ID = 903079
PMI_DONT_KNOW_ANSWER_ID = 903087
PMI_MISSING_ANSWER_IDS = {PMI_SKIP_ANSWER_ID, PMI_PREFER_NOT_ANSWER_ID, PMI_DONT_KNOW_ANSWER_ID}
PMI_NONRESPONSE_ANSWER_IDS = {PMI_PREFER_NOT_ANSWER_ID, PMI_DONT_KNOW_ANSWER_ID}
PFHH_FAMILY_CONDITION_QID = 43529217
PFHH_ADHD_ANY_ANSWER_ID = 43528365
PFHH_FAMILY_CONDITION_ANSWER_IDS = {
    903095,  # None
    836758,  # Alcohol use disorder
    836760,  # Drug use disorder
    836759,  # Autism spectrum disorder
    PFHH_ADHD_ANY_ANSWER_ID,
}
PFHH_ALCOHOL_RELATIVE_QID = 836827
PFHH_ALCOHOL_RELATIVE_ANSWER_IDS = {
    43528372,  # Father
    43528375,  # Sibling
    43528373,  # Grandparent
    43528374,  # Mother
    1384600,  # Self
    43528376,  # Son
    43528371,  # Daughter
}
PFHH_DRUG_RELATIVE_QID = 836851
PFHH_DRUG_RELATIVE_ANSWER_IDS = {
    43528657,  # Sibling
    1384413,  # Self
    43528654,  # Father
    43528656,  # Mother
    43528658,  # Son
    43528655,  # Grandparent
    43528653,  # Daughter
}
PFHH_ALLOWLIST_ANSWER_IDS = {
    PFHH_FAMILY_CONDITION_QID: PFHH_FAMILY_CONDITION_ANSWER_IDS,
    PFHH_ALCOHOL_RELATIVE_QID: PFHH_ALCOHOL_RELATIVE_ANSWER_IDS,
    PFHH_DRUG_RELATIVE_QID: PFHH_DRUG_RELATIVE_ANSWER_IDS,
}
AREA_SES_COLS = [
    "deprivation_index",
    "median_income",
    "fraction_poverty",
    "fraction_assisted_income",
    "fraction_no_health_ins",
    "fraction_vacant_housing",
]

SURVEY_SLUGS = {
    "The Basics": "basics",
    "Lifestyle": "lifestyle",
    "Overall Health": "overall_health",
    "Healthcare Access & Utilization": "hcau",
    "Personal and Family Health History": "pfhh",
    "Social Determinants of Health": "sdoh",
    "Behavioral Health and Personality": "bhp",
}

BHP_CODES = [f"bfi2xs_{i}" for i in range(1, 16)] + [f"asrs_{i}" for i in range(1, 7)]

MAIN_PRIMARY_IDS = {
    # The Basics.
    1585852,  # veteran / active duty
    903574, 903573, 903575, 903577, 903578, 903576,  # disability items
    1585952,  # employment
    1585389, 43528428, 1585386,  # insurance
    1585370, 1585375,  # home ownership, income
    1585402, 1585879, 1585889, 1585890, 1585886,  # living situation
    1585892, 1586135, 1585899, 1585357,  # marital, birthplace, orientation + closer description
    # Lifestyle: all codebook items.
    1586213, 1586198, 1586207, 1585870, 1586174, 1586177, 1586169, 1586166,
    1586185, 1586182, 1585656, 1585686, 1585674, 1585650, 1585704, 1585668,
    1585698, 903058, 1585680, 1585692, 1585636, 1586193, 1586190, 1585857,
    1586162, 1586159, 1585864, 1585873, 1585867, 1585860, 1586201,
    # Overall Health: PROMIS/global, BHLS, travel.
    1585748, 1585760, 1585741, 1585711, 1585729, 1585723, 1585717, 1585754,
    1585735, 1585766, 1585772, 1585778, 1585815,
    # Social Determinants of Health: all non-excluded SDOH IDs.
    40192400, 40192463, 40192499, 40192411, 40192417, 40192386, 40192469,
    40192500, 40192493, 40192420, 40192476, 40192457, 40192412, 40192404,
    40192456, 40192522, 40192384, 40192443, 40192471, 40192498, 40192401,
    40192501, 40192507, 40192398, 40192494, 40192397, 40192504, 40192415,
    40192390, 40192516, 40192475, 40192470, 40192439, 40192442, 40192511,
    40192446, 40192388, 40192480, 40192528, 40192399, 40192441, 40192449,
    40192396, 40192452, 40192419, 40192462, 40192491, 40192525, 40192445,
    40192381, 40192506, 40192440, 40192436, 40192410, 40192492, 40192414,
    40192431, 40192437, 40192402, 40192458, 40192426, 40192517,
    # HCAU: all non-race/religion-provider-concordance items.
    43528666, 43528665, 43530415, 43528662, 43528663, 43530408, 43530409,
    43528664, 43530413, 43530410, 43530411, 43530416, 43530412, 43530417,
    43530557, 43530583, 43529903, 43530585, 43529904, 43530584, 43530594,
    43530268, 43529905, 43529906, 43530437, 43530589, 43529974, 43530438,
    43530591, 43530588, 43529976, 43529977, 43529973, 43529975, 43530592,
    43530562, 43530590, 43530439, 43530399, 43530400, 43530403, 43528660,
    43528661, 43530402, 43530404, 43530401, 43530405, 43530406, 43530595,
    43530407, 43529978, 43530593, 43530559, 43530418,
    # PFHH mental health/substance-use family-history items.
    1740660, PFHH_FAMILY_CONDITION_QID, PFHH_ALCOHOL_RELATIVE_QID, PFHH_DRUG_RELATIVE_QID,
}

NUMERIC_IDS = {
    1585879, 1585889, 1585890,
    1586207, 1585870, 1586162, 1586159, 1585864, 1585873,
    40192441,
}

NONRESPONSE_INDICATOR_QIDS = {
    1585375,  # income
    1585636, 1585650, 1585656, 1585668, 1585674, 1585680, 1585686, 1585692,
    1585698, 1585704, 903058,  # substance use
    1585857, 1585860, 1586159, 1586162,  # smoking
    1586198, 1586201, 1586207, 1586213,  # alcohol
    1585886, 40192402, 40192426, 40192517,  # housing / food insecurity
}

LIFESTYLE_BRANCH_RECODE_ORDINAL = {
    (1585857, 1585859): {1585860: 0.0},  # never smoked 100 cigarettes -> smoke frequency not at all
    (1586166, 1586168): {1586169: 0.0},  # never used e-cigarettes -> current frequency not at all
    (1586174, 1586176): {1586177: 0.0},  # never smoked cigars -> current frequency not at all
    (1586182, 1586184): {1586185: 0.0},  # never smoked hookah -> current frequency not at all
    (1586190, 1586192): {1586193: 0.0},  # never used smokeless tobacco -> frequency not at all
    (1586198, 1586200): {1586201: 0.0, 1586213: 0.0},  # no alcohol -> never in frequency fields
    (1586201, 1586202): {1586213: 0.0},  # never drink -> never 6+ drinks
    (1585636, 1585648): {
        1585650: 0.0, 1585656: 0.0, 1585668: 0.0, 1585674: 0.0, 1585680: 0.0,
        1585686: 0.0, 1585692: 0.0, 1585698: 0.0, 1585704: 0.0, 903058: 0.0,
    },  # no listed recreational drugs -> no past-3-month use of each
}

LIFESTYLE_BRANCH_RECODE_NUMERIC = {
    (1585857, 1585859): {1586159: 0.0, 1586162: 0.0, 1585873: 0.0},
    (1585860, 1585863): {1586159: 0.0},
    (1586198, 1586200): {1586207: 0.0},
    (1586201, 1586202): {1586207: 0.0},
}

LEAKAGE_DENYLIST_RE = re.compile(r"education|school|grade|degree|ged|race|ethnic|discrimination", re.I)
LEAKAGE_WHITELIST_IDS = {
    1740660, PFHH_FAMILY_CONDITION_QID, PFHH_ALCOHOL_RELATIVE_QID, PFHH_DRUG_RELATIVE_QID,
    1586135, 1585766, 1585772, 1585778,
}


def log(lines, msg):
    print(msg, flush=True)
    lines.append(msg)


def sort_iids(iids):
    return sorted(iids, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))


def sanitize(value):
    value = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_").lower()
    return value[:80] or "feature"


def parse_int(value):
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def parse_float_answer(answer):
    if answer is None:
        return None
    text = str(answer).strip()
    if not text:
        return None
    if text.lower().startswith("pmi:"):
        return None
    if re.fullmatch(r"[-+]?\d+(\.\d+)?", text):
        return float(text)
    return None


def answer_tail(answer):
    text = re.sub(r"\s+", " ", str(answer or "").strip())
    if ":" in text:
        text = text.split(":")[-1].strip()
    return text.lower()


def ordinal_value_from_answer(answer):
    tail = answer_tail(answer)
    tail = tail.replace("dont know", "don't know")
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


def read_keep_iids(path):
    out = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                out.add(parts[-1])
    return out


def load_sex_covar(path):
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            iid = row["IID"].strip()
            if iid:
                out[iid] = int(row["sex_01"])
    return out


def load_fam_fids(path):
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                out[parts[1]] = parts[0]
    return out


def load_projected_pcs(path, n_pcs):
    pc_headers = [f"PC{i}_AVG" for i in range(1, n_pcs + 1)]
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header and header[0] == "#FID":
            header[0] = "FID"
        missing = [h for h in ["IID"] + pc_headers if h not in header]
        if missing:
            raise ValueError(f"{path} missing columns: {missing}")
        iid_idx = header.index("IID")
        pc_idx = [header.index(h) for h in pc_headers]
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(pc_idx):
                continue
            out[fields[iid_idx]] = [float(fields[i]) for i in pc_idx]
    return out, pc_headers


def load_ea_rows(path):
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"IID", "ea_years", "yob", "age_at_basics", "answer_concept_id", "answer"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            iid = row["IID"].strip()
            answer_id = int(row["answer_concept_id"])
            if answer_id not in EA_MAPPING:
                raise ValueError(f"Unexpected EA answer_concept_id={answer_id}")
            out[iid] = {
                "ea_years": EA_MAPPING[answer_id],
                "yob": float(row["yob"]),
                "age_at_basics": float(row["age_at_basics"]),
                "answer_concept_id": answer_id,
                "answer": row["answer"],
            }
    return out


def load_area_ses(path):
    out = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        lowered = {h.lower(): h for h in reader.fieldnames or []}
        iid_col = lowered.get("iid") or lowered.get("person_id")
        if not iid_col:
            raise ValueError(f"{path} missing IID/person_id column")
        for row in reader:
            iid = row[iid_col].strip()
            if not iid:
                continue
            vals = {}
            for col in AREA_SES_COLS:
                src = lowered.get(col)
                val = row[src].strip() if src else ""
                vals[col] = float(val) if val not in ("", "NA", "nan", "None") else math.nan
            out[iid] = vals
    return out


def load_local_metadata(path):
    meta = {}
    if not path or not os.path.exists(path):
        return meta
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            qid = parse_int(row.get("question_concept_id"))
            if qid is None:
                continue
            meta[qid] = {
                "survey": row.get("survey", ""),
                "question": row.get("question", ""),
            }
    return meta


def add_survey_feature_rows(path, feature_state, metadata, bhp=False):
    rows = 0
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows += 1
            iid = row["IID"].strip()
            if not iid:
                continue
            qid = int(row["question_concept_id"])
            survey = "Behavioral Health and Personality" if bhp else row["survey"]
            slug = SURVEY_SLUGS.get(survey, sanitize(survey))
            question = row.get("question") or row.get("question_name") or metadata.get(qid, {}).get("question", "")
            qcode = row.get("question_code", "")
            feature_state["survey_seen"][slug].add(iid)
            age = row.get("age_at_survey", "")
            try:
                age_val = float(age)
                prev = feature_state["survey_age"][slug].get(iid)
                if prev is None or age_val > prev:
                    feature_state["survey_age"][slug][iid] = age_val
            except ValueError:
                pass
            feature_state["question_meta"][qid] = {
                "survey": survey,
                "question": question,
                "question_code": qcode,
            }

            answer_id = parse_int(row.get("answer_concept_id"))
            answer = row.get("answer", "")
            if qid in PFHH_ALLOWLIST_ANSWER_IDS:
                if answer_id in PMI_MISSING_ANSWER_IDS:
                    feature_state["nonresponse"][(qid, answer_id)] += 1
                    feature_state["pmi_missing_iids"][(qid, answer_id)].add(iid)
                    continue
                feature_state["answered_qids"][qid].add(iid)
                if answer_id in PFHH_ALLOWLIST_ANSWER_IDS[qid]:
                    feature_state["selected_answers"][qid][iid].add(answer_id)
                    feature_state["answer_names"][(qid, answer_id)] = answer
                continue
            if answer_id in PMI_MISSING_ANSWER_IDS:
                feature_state["nonresponse"][(qid, answer_id)] += 1
                feature_state["pmi_missing_iids"][(qid, answer_id)].add(iid)
                continue
            numeric_value = parse_float_answer(answer)
            if qid in NUMERIC_IDS and numeric_value is not None:
                feature_state["numeric_values"][qid][iid] = numeric_value
                feature_state["answered_qids"][qid].add(iid)
                continue
            if answer_id is None:
                if numeric_value is not None:
                    feature_state["numeric_values"][qid][iid] = numeric_value
                    feature_state["answered_qids"][qid].add(iid)
                continue
            feature_state["answered_qids"][qid].add(iid)
            feature_state["selected_answers"][qid][iid].add(answer_id)
            feature_state["answer_names"][(qid, answer_id)] = answer
            ordinal_value = ordinal_value_from_answer(answer)
            if ordinal_value is not None:
                feature_state["ordinal_answer_values"][(qid, answer_id)] = ordinal_value
    return rows


def has_valid_answer(feature_state, iid, qid, answer_id):
    return answer_id in feature_state["selected_answers"].get(qid, {}).get(iid, set())


def has_any_valid_or_pmi(feature_state, iid, qid):
    if iid in feature_state["answered_qids"].get(qid, set()):
        return True
    for answer_id in PMI_MISSING_ANSWER_IDS:
        if iid in feature_state["pmi_missing_iids"].get((qid, answer_id), set()):
            return True
    return False


def apply_lifestyle_branch_recodes(feature_state):
    counts = Counter()
    for (parent_qid, parent_answer_id), children in LIFESTYLE_BRANCH_RECODE_ORDINAL.items():
        candidate_iids = feature_state["selected_answers"].get(parent_qid, {})
        for iid in candidate_iids:
            if not has_valid_answer(feature_state, iid, parent_qid, parent_answer_id):
                continue
            for child_qid, value in children.items():
                if has_any_valid_or_pmi(feature_state, iid, child_qid):
                    continue
                feature_state["branch_ordinal_values"][child_qid][iid] = value
                counts[(parent_qid, parent_answer_id, child_qid, "ordinal", value)] += 1
    for (parent_qid, parent_answer_id), children in LIFESTYLE_BRANCH_RECODE_NUMERIC.items():
        candidate_iids = feature_state["selected_answers"].get(parent_qid, {})
        for iid in candidate_iids:
            if not has_valid_answer(feature_state, iid, parent_qid, parent_answer_id):
                continue
            for child_qid, value in children.items():
                if has_any_valid_or_pmi(feature_state, iid, child_qid):
                    continue
                feature_state["branch_numeric_values"][child_qid][iid] = value
                counts[(parent_qid, parent_answer_id, child_qid, "numeric", value)] += 1
    rows = []
    for (parent_qid, parent_answer_id, child_qid, encoding, value), n in sorted(counts.items()):
        rows.append({
            "parent_question_concept_id": parent_qid,
            "parent_answer_concept_id": parent_answer_id,
            "child_question_concept_id": child_qid,
            "encoding": encoding,
            "value": value,
            "recoded_samples": n,
        })
    return rows


def qid_is_ordinal(feature_state, qid):
    if qid in PFHH_ALLOWLIST_ANSWER_IDS or qid in NUMERIC_IDS:
        return False
    selected = feature_state["selected_answers"].get(qid, {})
    for aids in selected.values():
        if len(aids) > 1:
            return False
    answer_ids = {aid for aids in selected.values() for aid in aids}
    if not answer_ids and qid not in feature_state["branch_ordinal_values"]:
        return False
    return all((qid, aid) in feature_state["ordinal_answer_values"] for aid in answer_ids)


def build_feature_matrix(iids, sex_map, area_ses, feature_state):
    iids = list(iids)
    iid_index = {iid: idx for idx, iid in enumerate(iids)}
    columns = []
    arrays = []

    def add_col(name, values):
        columns.append(name)
        arrays.append(np.asarray(values, dtype=np.float32))

    add_col("genetic_sex_01", [sex_map[iid] for iid in iids])

    all_slugs = ["basics", "lifestyle", "overall_health", "sdoh", "hcau", "pfhh", "bhp"]
    for slug in all_slugs:
        seen = feature_state["survey_seen"].get(slug, set())
        add_col(f"took_{slug}", [1.0 if iid in seen else 0.0 for iid in iids])
        ages = feature_state["survey_age"].get(slug, {})
        add_col(f"age_at_{slug}", [ages.get(iid, math.nan) for iid in iids])

    for col in AREA_SES_COLS:
        add_col(f"zip3_{col}", [area_ses.get(iid, {}).get(col, math.nan) for iid in iids])

    for qid in sorted(set(feature_state["numeric_values"]) | set(feature_state["branch_numeric_values"])):
        meta = feature_state["question_meta"].get(qid, {})
        qslug = sanitize(meta.get("question") or qid)
        vals = feature_state["numeric_values"][qid]
        branch_vals = feature_state["branch_numeric_values"].get(qid, {})
        feature_state["question_encoding"][qid] = "numeric"
        add_col(f"q{qid}_{qslug}_num", [vals.get(iid, branch_vals.get(iid, math.nan)) for iid in iids])

    for qid in sorted(set(feature_state["selected_answers"]) | set(feature_state["branch_ordinal_values"])):
        if not qid_is_ordinal(feature_state, qid):
            continue
        meta = feature_state["question_meta"].get(qid, {})
        qslug = sanitize(meta.get("question") or qid)
        selected_by_iid = feature_state["selected_answers"].get(qid, {})
        branch_vals = feature_state["branch_ordinal_values"].get(qid, {})
        values = []
        for iid in iids:
            if iid in selected_by_iid and selected_by_iid[iid]:
                answer_id = next(iter(selected_by_iid[iid]))
                values.append(feature_state["ordinal_answer_values"].get((qid, answer_id), math.nan))
            else:
                values.append(branch_vals.get(iid, math.nan))
        feature_state["question_encoding"][qid] = "ordinal"
        add_col(f"q{qid}_{qslug}_ord", values)

    for qid in sorted(NONRESPONSE_INDICATOR_QIDS & set(feature_state["question_meta"])):
        meta = feature_state["question_meta"].get(qid, {})
        qslug = sanitize(meta.get("question") or qid)
        valid = feature_state["answered_qids"].get(qid, set())
        pna = feature_state["pmi_missing_iids"].get((qid, PMI_PREFER_NOT_ANSWER_ID), set())
        dk = feature_state["pmi_missing_iids"].get((qid, PMI_DONT_KNOW_ANSWER_ID), set())
        values = []
        for iid in iids:
            if iid in pna or iid in dk:
                values.append(1.0)
            elif iid in valid:
                values.append(0.0)
            else:
                values.append(math.nan)
        add_col(f"q{qid}_{qslug}_nonresponse", values)

    for qid in sorted(feature_state["selected_answers"]):
        if qid_is_ordinal(feature_state, qid):
            continue
        selected_by_iid = feature_state["selected_answers"][qid]
        answered = feature_state["answered_qids"].get(qid, set())
        answer_ids = sorted({aid for aids in selected_by_iid.values() for aid in aids})
        if qid in PFHH_ALLOWLIST_ANSWER_IDS:
            answer_ids = sorted(PFHH_ALLOWLIST_ANSWER_IDS[qid])
        feature_state["question_encoding"][qid] = "one_hot"
        for answer_id in answer_ids:
            arr = np.full(len(iids), np.nan, dtype=np.float32)
            for iid in answered:
                idx = iid_index.get(iid)
                if idx is not None:
                    arr[idx] = 0.0
            for iid, aids in selected_by_iid.items():
                if answer_id in aids:
                    idx = iid_index.get(iid)
                    if idx is not None:
                        arr[idx] = 1.0
            answer_name = feature_state["answer_names"].get((qid, answer_id), "")
            add_col(f"q{qid}_a{answer_id}_{sanitize(answer_name)}", arr)

    if not arrays:
        raise RuntimeError("No features were built")
    matrix = np.vstack(arrays).T.astype(np.float32, copy=False)
    return matrix, columns


def residualize_and_z(ea, yob_c, sex_c, train_mask, apply_mask):
    inter = yob_c * sex_c
    design = np.column_stack([
        np.ones(len(ea), dtype=np.float64),
        yob_c.astype(np.float64),
        sex_c.astype(np.float64),
        inter.astype(np.float64),
    ])
    coef, *_ = np.linalg.lstsq(design[train_mask], ea[train_mask].astype(np.float64), rcond=None)
    resid_train = ea[train_mask] - design[train_mask].dot(coef)
    mean = float(np.mean(resid_train))
    sd = float(np.std(resid_train, ddof=1))
    if not math.isfinite(sd) or sd <= 0:
        raise RuntimeError("Invalid residual SD while building teacher label")
    out = np.full(len(ea), np.nan, dtype=np.float64)
    out[apply_mask] = (ea[apply_mask] - design[apply_mask].dot(coef) - mean) / sd
    return out, coef, mean, sd


def pearson(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 3:
        return math.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 3:
        return math.nan
    xr = pd.Series(x[mask]).rank(method="average").to_numpy()
    yr = pd.Series(y[mask]).rank(method="average").to_numpy()
    return float(np.corrcoef(xr, yr)[0, 1])


def write_tsv(path, rows, header):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compute_metric_row(label, pred, teacher_z, ea):
    r = pearson(pred, teacher_z)
    return {
        "group": label,
        "n": int(np.sum(np.isfinite(pred) & np.isfinite(teacher_z))),
        "pearson_proxy_teacher_z": r,
        "spearman_proxy_teacher_z": spearman(pred, teacher_z),
        "r2_proxy_teacher_z": r * r if math.isfinite(r) else math.nan,
        "pearson_proxy_ea_years": pearson(pred, ea),
        "spearman_proxy_ea_years": spearman(pred, ea),
    }


def read_final_model_kinship_exclusions(path, threshold, seed_iids, candidate_train_iids):
    """Return candidate training IIDs directly related to final applied seeds."""
    if not path:
        return set(), [], 0
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Missing final-model kinship table: {path}")

    seed_iids = set(seed_iids)
    candidate_train_iids = set(candidate_train_iids)
    excluded = set()
    edge_rows = []
    total_edges_ge_threshold = 0

    with open(path) as handle:
        header = handle.readline().strip().split()
        cols = {name.lstrip("#"): idx for idx, name in enumerate(header)}
        required = {"IID1", "IID2", "KINSHIP"}
        missing = required - set(cols)
        if missing:
            raise ValueError(f"{path} missing required columns: {sorted(missing)}")
        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            iid1 = parts[cols["IID1"]]
            iid2 = parts[cols["IID2"]]
            try:
                kinship = float(parts[cols["KINSHIP"]])
            except ValueError:
                continue
            if kinship < threshold:
                continue
            total_edges_ge_threshold += 1
            if iid1 in seed_iids and iid2 in candidate_train_iids:
                excluded.add(iid2)
                edge_rows.append({"seed_iid": iid1, "excluded_iid": iid2, "kinship": kinship})
            elif iid2 in seed_iids and iid1 in candidate_train_iids:
                excluded.add(iid1)
                edge_rows.append({"seed_iid": iid2, "excluded_iid": iid1, "kinship": kinship})
    edge_rows.sort(key=lambda r: (r["excluded_iid"], r["seed_iid"], -float(r["kinship"])))
    return excluded, edge_rows, total_edges_ge_threshold


def leakage_report_rows(feature_state, feature_columns):
    rows = []
    for qid, meta in sorted(feature_state["question_meta"].items()):
        text = f"{meta.get('question_code', '')} {meta.get('question', '')}"
        hit = bool(LEAKAGE_DENYLIST_RE.search(text))
        whitelisted = qid in LEAKAGE_WHITELIST_IDS
        included = qid in MAIN_PRIMARY_IDS or meta.get("question_code") in BHP_CODES
        if hit:
            rows.append({
                "question_concept_id": qid,
                "question_code": meta.get("question_code", ""),
                "question": meta.get("question", ""),
                "included": int(included),
                "whitelisted": int(whitelisted),
                "status": "allowed" if whitelisted else "denylist_hit",
            })
    for col in feature_columns:
        if "fraction_high_school_edu" in col:
            rows.append({
                "question_concept_id": "zip3",
                "question_code": "fraction_high_school_edu",
                "question": "area-level education",
                "included": 1,
                "whitelisted": 0,
                "status": "denylist_hit",
            })
    return rows


def build_feature_manifest(feature_state, feature_columns):
    rows = []
    for qid, meta in sorted(feature_state["question_meta"].items()):
        qcode = meta.get("question_code", "")
        include = qid in MAIN_PRIMARY_IDS or qcode in BHP_CODES
        encoding = feature_state["question_encoding"].get(qid)
        if encoding is None:
            encoding = "numeric" if qid in NUMERIC_IDS else "one_hot"
        if qid in PFHH_ALLOWLIST_ANSWER_IDS:
            encoding = "allowlisted_one_hot"
        rows.append({
            "question_concept_id": qid,
            "survey": meta.get("survey", ""),
            "item_name": meta.get("question", ""),
            "encoding": encoding,
            "include_exclude": "include" if include else "exclude",
            "external_analog": "",
            "notes": "primary_full feature" if include else "not used",
        })
    for col in AREA_SES_COLS:
        rows.append({
            "question_concept_id": "zip3_ses",
            "survey": "zip3_ses_map",
            "item_name": col,
            "encoding": "numeric",
            "include_exclude": "include",
            "external_analog": "Townsend deprivation index" if col == "deprivation_index" else "",
            "notes": "derived area-SES; raw ZIP not used",
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ea-query", required=True)
    parser.add_argument("--main-survey", required=True)
    parser.add_argument("--bhp-survey", required=True)
    parser.add_argument("--area-ses", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--europeans", required=True)
    parser.add_argument("--fit-pca-iids", required=True)
    parser.add_argument("--sex-covar", required=True)
    parser.add_argument("--exclude-iids", required=True)
    parser.add_argument("--fam", required=True)
    parser.add_argument("--sscore", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-pcs", type=int, default=10)
    parser.add_argument("--min-age-at-basics", type=float, default=26.0)
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--num-boost-round", type=int, default=2000)
    parser.add_argument("--early-stopping-rounds", type=int, default=50)
    parser.add_argument("--cv-folds", type=int, default=4)
    parser.add_argument("--final-kinship-holdout-kin0", default="")
    parser.add_argument("--final-kinship-holdout-threshold", type=float, default=0.0441941)
    args = parser.parse_args()

    start = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    model_dir = os.path.join(args.out_dir, "xgboost_models")
    os.makedirs(model_dir, exist_ok=True)
    log_lines = []
    log(log_lines, "=== AoU ses_ea_proxy primary setup ===")
    log(log_lines, f"out_dir: {args.out_dir}")
    log(log_lines, f"outer_folds: {args.outer_folds}")
    log(log_lines, f"seed: {args.seed}")
    log(log_lines, f"threads: {args.threads}")

    try:
        import xgboost as xgb
    except ImportError as exc:
        raise RuntimeError("xgboost is required; install it before running this helper") from exc

    europeans = read_keep_iids(args.europeans)
    fit_pca = read_keep_iids(args.fit_pca_iids)
    sex_map = load_sex_covar(args.sex_covar)
    excluded_iids = read_keep_iids(args.exclude_iids)
    fid_by_iid = load_fam_fids(args.fam)
    pc_data, pc_headers = load_projected_pcs(args.sscore, args.n_pcs)
    ea_rows = load_ea_rows(args.ea_query)
    area_ses = load_area_ses(args.area_ses)
    metadata = load_local_metadata(args.metadata)

    log(log_lines, "")
    log(log_lines, "=== Input counts ===")
    for name, value in [
        ("classified Europeans", len(europeans)),
        ("fit_pca_iids", len(fit_pca)),
        ("sex covariate rows", len(sex_map)),
        ("sample-QC exclusion IIDs", len(excluded_iids)),
        ("fam rows", len(fid_by_iid)),
        ("projected PC rows", len(pc_data)),
        ("EA query rows", len(ea_rows)),
        ("area SES rows", len(area_ses)),
    ]:
        log(log_lines, f"{name}: {value}")

    candidates = set(europeans) - excluded_iids
    after_sample_qc = set(candidates)
    missing_ea = candidates - set(ea_rows)
    candidates &= set(ea_rows)
    below_min_age = {iid for iid in candidates if ea_rows[iid]["age_at_basics"] < args.min_age_at_basics}
    candidates -= below_min_age
    missing_sex = candidates - set(sex_map)
    candidates &= set(sex_map)
    missing_fam = candidates - set(fid_by_iid)
    candidates &= set(fid_by_iid)
    missing_pcs = candidates - set(pc_data)
    candidates &= set(pc_data)

    eligible_iids = sort_iids(candidates)
    train_iids = sort_iids(set(eligible_iids) & fit_pca)
    applied_iids = sort_iids(set(eligible_iids) - set(train_iids))
    if len(train_iids) < args.outer_folds:
        raise RuntimeError("Too few eligible fit_pca_iids for cross-fitting")
    if not eligible_iids:
        raise RuntimeError("No eligible samples remain")

    log(log_lines, "")
    log(log_lines, "=== Eligibility counts ===")
    log(log_lines, f"Europeans removed by sample-QC exclusion: {len(europeans & excluded_iids)}")
    log(log_lines, f"Europeans after sample-QC exclusion: {len(after_sample_qc)}")
    log(log_lines, f"Europeans missing codeable EA: {len(missing_ea)}")
    log(log_lines, f"EA candidates below age {args.min_age_at_basics:g}: {len(below_min_age)}")
    log(log_lines, f"EA/age candidates missing confirmed genetic sex: {len(missing_sex)}")
    log(log_lines, f"EA/age/sex candidates missing fam row: {len(missing_fam)}")
    log(log_lines, f"EA/age/sex/fam candidates missing PCs: {len(missing_pcs)}")
    log(log_lines, f"Final eligible classified EUR samples: {len(eligible_iids)}")
    log(log_lines, f"Eligible fit_pca training samples: {len(train_iids)}")
    log(log_lines, f"Eligible applied classified-EUR extras: {len(applied_iids)}")

    final_excluded_iids, final_exclusion_edges, final_total_edges_ge_threshold = read_final_model_kinship_exclusions(
        args.final_kinship_holdout_kin0,
        args.final_kinship_holdout_threshold,
        applied_iids,
        train_iids,
    )
    final_train_iids = sort_iids(set(train_iids) - final_excluded_iids)
    if args.final_kinship_holdout_kin0 and not final_train_iids:
        raise RuntimeError("No final-model training samples remain after kinship holdout")
    final_train_iid_set = set(final_train_iids)
    final_train_allowed_by_iid = {iid: iid in final_train_iid_set for iid in eligible_iids}
    log(log_lines, "")
    log(log_lines, "=== Final applied-model kinship holdout ===")
    log(log_lines, f"Kinship table: {args.final_kinship_holdout_kin0 or 'not used'}")
    log(log_lines, f"Kinship threshold: {args.final_kinship_holdout_threshold:g}")
    log(log_lines, f"Applied seed samples: {len(applied_iids)}")
    log(log_lines, f"KING edges at or above threshold: {final_total_edges_ge_threshold}")
    log(log_lines, f"Eligible fit_pca samples excluded as applied relatives: {len(final_excluded_iids)}")
    log(log_lines, f"Final applied-model training samples: {len(final_train_iids)}")

    feature_state = {
        "survey_seen": defaultdict(set),
        "survey_age": defaultdict(dict),
        "question_meta": {},
        "answered_qids": defaultdict(set),
        "selected_answers": defaultdict(lambda: defaultdict(set)),
        "answer_names": {},
        "numeric_values": defaultdict(dict),
        "ordinal_answer_values": {},
        "branch_ordinal_values": defaultdict(dict),
        "branch_numeric_values": defaultdict(dict),
        "pmi_missing_iids": defaultdict(set),
        "nonresponse": Counter(),
        "question_encoding": {},
    }
    main_rows = add_survey_feature_rows(args.main_survey, feature_state, metadata, bhp=False)
    bhp_rows = add_survey_feature_rows(args.bhp_survey, feature_state, metadata, bhp=True)
    branch_recode_rows = apply_lifestyle_branch_recodes(feature_state)
    branch_recode_total = sum(int(row["recoded_samples"]) for row in branch_recode_rows)
    log(log_lines, "")
    log(log_lines, "=== Feature extraction counts ===")
    log(log_lines, f"Main survey extracted rows: {main_rows}")
    log(log_lines, f"BHP extracted rows: {bhp_rows}")
    log(log_lines, f"Lifestyle branch recode rules with samples: {len(branch_recode_rows)}")
    log(log_lines, f"Lifestyle branch recoded sample-feature cells: {branch_recode_total}")
    log(log_lines, f"Survey questions with extracted rows: {len(feature_state['question_meta'])}")
    log(log_lines, f"Numeric feature questions: {len(feature_state['numeric_values'])}")
    log(log_lines, f"Categorical feature questions: {len(feature_state['selected_answers'])}")

    X_all, feature_columns = build_feature_matrix(eligible_iids, sex_map, area_ses, feature_state)
    feature_columns_blob = "\n".join(feature_columns).encode("utf-8")
    feature_columns_sha256 = hashlib.sha256(feature_columns_blob).hexdigest()
    with open(os.path.join(args.out_dir, "xgboost_feature_columns.json"), "w") as f:
        json.dump(feature_columns, f, indent=2)
    write_tsv(
        os.path.join(args.out_dir, "xgboost_feature_columns.tsv"),
        [{"feature_index": i, "feature": feature} for i, feature in enumerate(feature_columns)],
        ["feature_index", "feature"],
    )
    iid_to_row = {iid: idx for idx, iid in enumerate(eligible_iids)}
    train_idx = np.asarray([iid_to_row[iid] for iid in train_iids], dtype=np.int64)
    applied_idx = np.asarray([iid_to_row[iid] for iid in applied_iids], dtype=np.int64)
    final_train_idx = np.asarray([iid_to_row[iid] for iid in final_train_iids], dtype=np.int64)

    ea = np.asarray([ea_rows[iid]["ea_years"] for iid in eligible_iids], dtype=np.float64)
    yob = np.asarray([ea_rows[iid]["yob"] for iid in eligible_iids], dtype=np.float64)
    sex = np.asarray([sex_map[iid] for iid in eligible_iids], dtype=np.float64)
    mean_yob = float(np.mean(yob[train_idx]))
    yob_c = yob - mean_yob
    sex_c = sex - 0.5
    pcs = np.asarray([pc_data[iid] for iid in eligible_iids], dtype=np.float64)

    rng = np.random.default_rng(args.seed)
    shuffled_train = np.asarray(train_iids, dtype=object)
    fold_ids = rng.integers(0, args.outer_folds, size=len(shuffled_train))
    # Ensure no empty folds for unusual small test runs.
    if len(set(fold_ids.tolist())) < args.outer_folds:
        fold_ids = np.arange(len(shuffled_train)) % args.outer_folds
        rng.shuffle(fold_ids)
    fold_by_iid = {iid: int(fold) for iid, fold in zip(shuffled_train, fold_ids)}

    params = {
        "objective": "reg:squarederror",
        "eta": 0.05,
        "max_depth": 6,
        "min_child_weight": 20,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "lambda": 1.0,
        "alpha": 0.0,
        "tree_method": "hist",
        "seed": args.seed,
        "nthread": args.threads,
        "eval_metric": "rmse",
    }

    pred_raw = np.full(len(eligible_iids), np.nan, dtype=np.float64)
    teacher_z = np.full(len(eligible_iids), np.nan, dtype=np.float64)
    best_rounds = []
    fold_metric_rows = []
    fit_rows = []
    model_manifest_rows = []

    log(log_lines, "")
    log(log_lines, "=== XGBoost OOF training ===")
    for fold in range(args.outer_folds):
        test_iids = [iid for iid in train_iids if fold_by_iid[iid] == fold]
        train_pool_iids = [iid for iid in train_iids if fold_by_iid[iid] != fold]
        test_idx = np.asarray([iid_to_row[iid] for iid in test_iids], dtype=np.int64)
        train_pool_idx = np.asarray([iid_to_row[iid] for iid in train_pool_iids], dtype=np.int64)
        train_mask = np.zeros(len(eligible_iids), dtype=bool)
        apply_mask = np.zeros(len(eligible_iids), dtype=bool)
        train_mask[train_pool_idx] = True
        apply_mask[train_pool_idx] = True
        apply_mask[test_idx] = True
        z_all, coef, resid_mean, resid_sd = residualize_and_z(ea, yob_c, sex_c, train_mask, apply_mask)
        y_train = z_all[train_pool_idx]
        teacher_z[test_idx] = z_all[test_idx]
        dtrain = xgb.DMatrix(X_all[train_pool_idx, :], label=y_train, feature_names=feature_columns, missing=np.nan)
        cv = xgb.cv(
            params,
            dtrain,
            num_boost_round=args.num_boost_round,
            nfold=min(args.cv_folds, max(2, len(train_pool_idx) // 2)),
            early_stopping_rounds=args.early_stopping_rounds,
            seed=args.seed + fold,
            verbose_eval=False,
        )
        best_round = int(len(cv))
        best_rounds.append(best_round)
        model = xgb.train(params, dtrain, num_boost_round=best_round, verbose_eval=False)
        model.set_attr(
            role="oof_fold",
            fold_id=str(fold),
            feature_columns_sha256=feature_columns_sha256,
            train_pool_samples=str(len(train_pool_idx)),
            prediction_samples=str(len(test_idx)),
            num_boost_round=str(best_round),
        )
        model_file = os.path.join("xgboost_models", f"fold_{fold}.json")
        model.save_model(os.path.join(args.out_dir, model_file))
        dtest = xgb.DMatrix(X_all[test_idx, :], feature_names=feature_columns, missing=np.nan)
        pred_raw[test_idx] = model.predict(dtest)
        row = compute_metric_row(f"fold_{fold}", pred_raw[test_idx], teacher_z[test_idx], ea[test_idx])
        row.update({
            "best_rounds": best_round,
            "train_pool_samples": len(train_pool_idx),
            "test_samples": len(test_idx),
            "ols_intercept": coef[0],
            "ols_yob_c": coef[1],
            "ols_sex_c": coef[2],
            "ols_yob_c_sex_c_inter": coef[3],
            "resid_train_mean": resid_mean,
            "resid_train_sd": resid_sd,
        })
        fold_metric_rows.append(row)
        model_manifest_rows.append({
            "model_name": f"fold_{fold}",
            "role": "oof_fold",
            "fold_id": fold,
            "model_file": model_file,
            "feature_columns_json": "xgboost_feature_columns.json",
            "feature_columns_sha256": feature_columns_sha256,
            "num_boost_round": best_round,
            "train_pool_samples": len(train_pool_idx),
            "prediction_samples": len(test_idx),
            "ols_intercept": coef[0],
            "ols_yob_c": coef[1],
            "ols_sex_c": coef[2],
            "ols_yob_c_sex_c_inter": coef[3],
            "resid_train_mean": resid_mean,
            "resid_train_sd": resid_sd,
        })
        log(log_lines, f"fold {fold}: train_pool={len(train_pool_idx)} test={len(test_idx)} best_rounds={best_round} r={row['pearson_proxy_teacher_z']:.6f}")

    final_rounds = int(statistics.median(best_rounds))
    train_mask = np.zeros(len(eligible_iids), dtype=bool)
    apply_mask = np.zeros(len(eligible_iids), dtype=bool)
    train_mask[final_train_idx] = True
    apply_mask[final_train_idx] = True
    apply_mask[applied_idx] = True
    final_teacher_z, final_coef, final_resid_mean, final_resid_sd = residualize_and_z(
        ea, yob_c, sex_c, train_mask, apply_mask
    )
    teacher_z[applied_idx] = final_teacher_z[applied_idx]
    dtrain_final = xgb.DMatrix(
        X_all[final_train_idx, :],
        label=final_teacher_z[final_train_idx],
        feature_names=feature_columns,
        missing=np.nan,
    )
    final_model = xgb.train(params, dtrain_final, num_boost_round=final_rounds, verbose_eval=False)
    final_model.set_attr(
        role="applied_final_model",
        fold_id="final_model",
        feature_columns_sha256=feature_columns_sha256,
        train_pool_samples=str(len(final_train_idx)),
        prediction_samples=str(len(applied_idx)),
        num_boost_round=str(final_rounds),
        final_kinship_holdout_kin0=str(args.final_kinship_holdout_kin0),
        final_kinship_holdout_threshold=str(args.final_kinship_holdout_threshold),
        final_kinship_excluded_fit_pca_samples=str(len(final_excluded_iids)),
    )
    final_model_file = os.path.join("xgboost_models", "final_model.json")
    final_model.save_model(os.path.join(args.out_dir, final_model_file))
    model_manifest_rows.append({
        "model_name": "final_model",
        "role": "applied_final_model",
        "fold_id": "final_model",
        "model_file": final_model_file,
        "feature_columns_json": "xgboost_feature_columns.json",
        "feature_columns_sha256": feature_columns_sha256,
        "num_boost_round": final_rounds,
        "train_pool_samples": len(final_train_idx),
        "prediction_samples": len(applied_idx),
        "ols_intercept": final_coef[0],
        "ols_yob_c": final_coef[1],
        "ols_sex_c": final_coef[2],
        "ols_yob_c_sex_c_inter": final_coef[3],
        "resid_train_mean": final_resid_mean,
        "resid_train_sd": final_resid_sd,
        "final_kinship_holdout_kin0": args.final_kinship_holdout_kin0,
        "final_kinship_holdout_threshold": args.final_kinship_holdout_threshold,
        "final_applied_seed_samples": len(applied_iids),
        "final_kinship_excluded_fit_pca_samples": len(final_excluded_iids),
        "final_model_train_allowed_samples": len(final_train_iids),
    })
    if len(applied_idx):
        dapplied = xgb.DMatrix(X_all[applied_idx, :], feature_names=feature_columns, missing=np.nan)
        pred_raw[applied_idx] = final_model.predict(dapplied)
    log(log_lines, f"final model: train={len(train_idx)} applied={len(applied_idx)} rounds={final_rounds}")

    pred_mean = float(np.nanmean(pred_raw))
    pred_sd = float(np.nanstd(pred_raw, ddof=1))
    if not math.isfinite(pred_sd) or pred_sd <= 0:
        raise RuntimeError("Invalid prediction SD")
    proxy_z = (pred_raw - pred_mean) / pred_sd

    for iid in train_iids:
        idx = iid_to_row[iid]
        fit_rows.append({
            "FID": fid_by_iid[iid],
            "IID": iid,
            "role": "oof",
            "fold_id": fold_by_iid[iid],
            "final_model_train_allowed": int(final_train_allowed_by_iid[iid]),
            "ea_years": ea[idx],
            "teacher_z": teacher_z[idx],
            "score_raw": pred_raw[idx],
            "ses_ea_proxy_z": proxy_z[idx],
        })
    applied_rows = []
    for iid in applied_iids:
        idx = iid_to_row[iid]
        applied_rows.append({
            "FID": fid_by_iid[iid],
            "IID": iid,
            "role": "applied",
            "fold_id": "final_model",
            "final_model_train_allowed": 0,
            "ea_years": ea[idx],
            "teacher_z": teacher_z[idx],
            "score_raw": pred_raw[idx],
            "ses_ea_proxy_z": proxy_z[idx],
        })

    applied_metrics = []
    if len(applied_idx):
        row = compute_metric_row("applied_6th_model", proxy_z[applied_idx], teacher_z[applied_idx], ea[applied_idx])
        row.update({
            "best_rounds": final_rounds,
            "train_pool_samples": len(final_train_idx),
            "test_samples": len(applied_idx),
            "ols_intercept": final_coef[0],
            "ols_yob_c": final_coef[1],
            "ols_sex_c": final_coef[2],
            "ols_yob_c_sex_c_inter": final_coef[3],
            "resid_train_mean": final_resid_mean,
            "resid_train_sd": final_resid_sd,
        })
        applied_metrics.append(row)

    covariate_rows = []
    group_masks = {}
    for fold in range(args.outer_folds):
        mask = np.zeros(len(eligible_iids), dtype=bool)
        mask[[iid_to_row[iid] for iid in train_iids if fold_by_iid[iid] == fold]] = True
        group_masks[f"fold_{fold}"] = mask
    mask = np.zeros(len(eligible_iids), dtype=bool)
    mask[train_idx] = True
    group_masks["oof_overall"] = mask
    mask = np.zeros(len(eligible_iids), dtype=bool)
    mask[applied_idx] = True
    group_masks["applied_6th_model"] = mask
    group_masks["combined_overall"] = np.isfinite(proxy_z)
    covars = {
        "teacher_z": teacher_z,
        "ea_years": ea,
        "yob_c": yob_c,
        "sex_c": sex_c,
    }
    for i, pc in enumerate(pc_headers):
        covars[pc] = pcs[:, i]
    for group, mask in group_masks.items():
        for covar, values in covars.items():
            covariate_rows.append({
                "group": group,
                "covariate": covar,
                "n": int(np.sum(mask & np.isfinite(proxy_z) & np.isfinite(values))),
                "pearson": pearson(proxy_z[mask], values[mask]),
                "spearman": spearman(proxy_z[mask], values[mask]),
            })

    all_score_rows = fit_rows + applied_rows
    all_score_rows.sort(key=lambda r: (0, int(r["IID"])) if str(r["IID"]).isdigit() else (1, r["IID"]))
    score_header = [
        "FID", "IID", "role", "fold_id", "final_model_train_allowed",
        "ea_years", "teacher_z", "score_raw", "ses_ea_proxy_z",
    ]
    write_tsv(os.path.join(args.out_dir, "oof_scores.tsv"), fit_rows, score_header)
    write_tsv(os.path.join(args.out_dir, "applied_scores.tsv"), applied_rows, score_header)
    write_tsv(os.path.join(args.out_dir, "all_scores.tsv"), all_score_rows, score_header)

    with open(os.path.join(args.out_dir, "phen.txt"), "w") as f:
        f.write("FID\tIID\tses_ea_proxy_z\n")
        for row in all_score_rows:
            f.write(f"{row['FID']}\t{row['IID']}\t{float(row['ses_ea_proxy_z']):.10g}\n")

    with open(os.path.join(args.out_dir, "training_iids.txt"), "w") as f:
        for row in all_score_rows:
            f.write(f"{row['FID']} {row['IID']}\n")

    covar_by_iid = {}
    for iid in eligible_iids:
        idx = iid_to_row[iid]
        inter = yob_c[idx] * sex_c[idx]
        covar_by_iid[iid] = [yob_c[idx], sex_c[idx], inter] + pc_data[iid]
    covar_cols = ["yob_c", "sex_c", "yob_c_sex_c_inter"] + pc_headers
    with open(os.path.join(args.out_dir, "base_covar.txt"), "w") as f:
        f.write("FID\tIID\tyob_c\tsex_c\tyob_c_sex_c_inter\n")
        for row in all_score_rows:
            vals = covar_by_iid[row["IID"]]
            f.write(f"{row['FID']}\t{row['IID']}\t{vals[0]:.12g}\t{vals[1]:.1f}\t{vals[2]:.12g}\n")
    with open(os.path.join(args.out_dir, "covar.txt"), "w") as f:
        f.write("FID\tIID\t" + "\t".join(covar_cols) + "\n")
        for row in all_score_rows:
            vals = covar_by_iid[row["IID"]]
            f.write(f"{row['FID']}\t{row['IID']}\t" + "\t".join(f"{v:.12g}" for v in vals) + "\n")

    fold_assignment_rows = []
    for iid in train_iids:
        fold_assignment_rows.append({"FID": fid_by_iid[iid], "IID": iid, "role": "oof", "fold_id": fold_by_iid[iid]})
    for iid in applied_iids:
        fold_assignment_rows.append({"FID": fid_by_iid[iid], "IID": iid, "role": "applied", "fold_id": "final_model"})
    write_tsv(os.path.join(args.out_dir, "fold_assignment.tsv"), fold_assignment_rows, ["FID", "IID", "role", "fold_id"])

    metric_header = [
        "group", "n", "pearson_proxy_teacher_z", "spearman_proxy_teacher_z", "r2_proxy_teacher_z",
        "pearson_proxy_ea_years", "spearman_proxy_ea_years", "best_rounds", "train_pool_samples",
        "test_samples", "ols_intercept", "ols_yob_c", "ols_sex_c", "ols_yob_c_sex_c_inter",
        "resid_train_mean", "resid_train_sd",
    ]
    write_tsv(os.path.join(args.out_dir, "fold_metrics.tsv"), fold_metric_rows, metric_header)
    write_tsv(os.path.join(args.out_dir, "applied_metrics.tsv"), applied_metrics, metric_header)
    write_tsv(os.path.join(args.out_dir, "proxy_covariate_correlations.tsv"), covariate_rows,
              ["group", "covariate", "n", "pearson", "spearman"])
    write_tsv(os.path.join(args.out_dir, "xgboost_model_manifest.tsv"), model_manifest_rows,
              ["model_name", "role", "fold_id", "model_file", "feature_columns_json",
               "feature_columns_sha256", "num_boost_round", "train_pool_samples",
               "prediction_samples", "ols_intercept", "ols_yob_c", "ols_sex_c",
               "ols_yob_c_sex_c_inter", "resid_train_mean", "resid_train_sd",
               "final_kinship_holdout_kin0", "final_kinship_holdout_threshold",
               "final_applied_seed_samples", "final_kinship_excluded_fit_pca_samples",
               "final_model_train_allowed_samples"])

    with open(os.path.join(args.out_dir, "final_model_train_iids.txt"), "w") as f:
        for iid in final_train_iids:
            f.write(f"{fid_by_iid[iid]} {iid}\n")
    with open(os.path.join(args.out_dir, "final_model_excluded_related_to_applied_iids.txt"), "w") as f:
        for iid in sort_iids(final_excluded_iids):
            f.write(f"{fid_by_iid.get(iid, iid)} {iid}\n")
    write_tsv(
        os.path.join(args.out_dir, "final_model_excluded_related_to_applied_edges.tsv"),
        final_exclusion_edges,
        ["seed_iid", "excluded_iid", "kinship"],
    )
    write_tsv(
        os.path.join(args.out_dir, "final_model_kinholdout_summary.tsv"),
        [
            {"metric": "final_kinship_holdout_kin0", "value": args.final_kinship_holdout_kin0},
            {"metric": "final_kinship_holdout_threshold", "value": args.final_kinship_holdout_threshold},
            {"metric": "final_applied_seed_samples", "value": len(applied_iids)},
            {"metric": "final_candidate_fit_pca_samples", "value": len(train_iids)},
            {"metric": "final_kinship_edges_ge_threshold", "value": final_total_edges_ge_threshold},
            {"metric": "final_kinship_excluded_fit_pca_samples", "value": len(final_excluded_iids)},
            {"metric": "final_model_train_allowed_samples", "value": len(final_train_iids)},
        ],
        ["metric", "value"],
    )

    importance_rows = []
    importance_gain = final_model.get_score(importance_type="gain")
    importance_cover = final_model.get_score(importance_type="cover")
    for feature in sorted(set(importance_gain) | set(importance_cover)):
        importance_rows.append({
            "feature": feature,
            "gain": importance_gain.get(feature, 0.0),
            "cover": importance_cover.get(feature, 0.0),
        })
    importance_rows.sort(key=lambda r: (-float(r["gain"]), r["feature"]))
    write_tsv(os.path.join(args.out_dir, "feature_importance.tsv"), importance_rows, ["feature", "gain", "cover"])

    feature_manifest = build_feature_manifest(feature_state, feature_columns)
    write_tsv(os.path.join(args.out_dir, "feature_manifest.resolved.tsv"), feature_manifest,
              ["question_concept_id", "survey", "item_name", "encoding", "include_exclude", "external_analog", "notes"])
    leakage_rows = leakage_report_rows(feature_state, feature_columns)
    write_tsv(os.path.join(args.out_dir, "leakage_denylist_hits.tsv"), leakage_rows,
              ["question_concept_id", "question_code", "question", "included", "whitelisted", "status"])
    bad_leakage = [r for r in leakage_rows if r["included"] and not r["whitelisted"] and r["status"] == "denylist_hit"]
    if bad_leakage:
        raise RuntimeError(f"Leakage denylist found unapproved included features: {bad_leakage[:5]}")

    feature_count_rows = [
        {"metric": "feature_columns", "value": len(feature_columns)},
        {"metric": "numeric_feature_columns", "value": sum(1 for c in feature_columns if c.endswith("_num") or c.startswith("zip3_") or c.startswith("age_at_") or c == "genetic_sex_01")},
        {"metric": "ordinal_feature_columns", "value": sum(1 for c in feature_columns if c.endswith("_ord"))},
        {"metric": "one_hot_feature_columns", "value": sum(1 for c in feature_columns if re.search(r"_a[0-9]+_", c))},
        {"metric": "nonresponse_indicator_columns", "value": sum(1 for c in feature_columns if c.endswith("_nonresponse"))},
        {"metric": "survey_took_indicators", "value": sum(1 for c in feature_columns if c.startswith("took_"))},
        {"metric": "lifestyle_branch_recode_rules_with_samples", "value": len(branch_recode_rows)},
        {"metric": "lifestyle_branch_recoded_sample_feature_cells", "value": branch_recode_total},
    ]
    write_tsv(os.path.join(args.out_dir, "feature_counts.tsv"), feature_count_rows, ["metric", "value"])

    pmi_rows = []
    answer_label = {
        PMI_SKIP_ANSWER_ID: "pmi_skip",
        PMI_PREFER_NOT_ANSWER_ID: "prefer_not_to_answer",
        PMI_DONT_KNOW_ANSWER_ID: "dont_know",
    }
    for (qid, answer_id), n_rows in sorted(feature_state["nonresponse"].items()):
        meta = feature_state["question_meta"].get(qid, {})
        pmi_rows.append({
            "question_concept_id": qid,
            "survey": meta.get("survey", ""),
            "item_name": meta.get("question", ""),
            "answer_concept_id": answer_id,
            "missing_kind": answer_label.get(answer_id, "other_missing"),
            "answer_rows": n_rows,
            "distinct_samples": len(feature_state["pmi_missing_iids"].get((qid, answer_id), set())),
            "value_encoding": "NaN",
            "nonresponse_indicator": int(qid in NONRESPONSE_INDICATOR_QIDS and answer_id in PMI_NONRESPONSE_ANSWER_IDS),
        })
    write_tsv(os.path.join(args.out_dir, "pmi_missingness_counts.tsv"), pmi_rows,
              ["question_concept_id", "survey", "item_name", "answer_concept_id", "missing_kind",
               "answer_rows", "distinct_samples", "value_encoding", "nonresponse_indicator"])

    write_tsv(os.path.join(args.out_dir, "branch_recoding_summary.tsv"), branch_recode_rows,
              ["parent_question_concept_id", "parent_answer_concept_id", "child_question_concept_id",
               "encoding", "value", "recoded_samples"])

    missing_policy_rows = [
        {
            "case": "did_not_take_survey",
            "handling": "survey-block item features remain NaN; took_<survey>=0; age_at_<survey>=NaN",
        },
        {
            "case": "took_survey",
            "handling": "took_<survey>=1; age_at_<survey> is the response age from the latest retained survey instance",
        },
        {
            "case": "latest_instance_multi_select",
            "handling": "select latest person-question timestamp, keep all answer rows at that timestamp, one-hot each answer_concept_id",
        },
        {
            "case": "pmi_skip_903096",
            "handling": "value remains NaN for native XGBoost missing routing; counted separately; no curated nonresponse flag",
        },
        {
            "case": "pmi_prefer_not_903079_or_dont_know_903087",
            "handling": "value remains NaN; curated high-value questions also get *_nonresponse=1; valid answers get *_nonresponse=0; not asked remains NaN",
        },
        {
            "case": "lifestyle_branch_not_applicable",
            "handling": "never/no parent answers set downstream smoking/alcohol/substance frequency or count features to 0 instead of NaN",
        },
        {
            "case": "xgboost_missing",
            "handling": "all matrices are built with DMatrix(missing=np.nan), so unobserved survey blocks route natively through trees",
        },
    ]
    write_tsv(os.path.join(args.out_dir, "missing_data_handling.tsv"), missing_policy_rows,
              ["case", "handling"])

    missing_rows = []
    for j, col in enumerate(feature_columns):
        vals = X_all[:, j]
        missing_rows.append({
            "feature": col,
            "missing": int(np.sum(~np.isfinite(vals))),
            "nonmissing": int(np.sum(np.isfinite(vals))),
            "missing_fraction": float(np.mean(~np.isfinite(vals))),
        })
    missing_rows.sort(key=lambda r: (-r["missing_fraction"], r["feature"]))
    write_tsv(os.path.join(args.out_dir, "feature_missingness.tsv"), missing_rows,
              ["feature", "missing", "nonmissing", "missing_fraction"])

    summary_path = os.path.join(args.out_dir, "ses_ea_proxy_gwas.summary.tsv")
    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"classified_europeans\t{len(europeans)}\n")
        f.write(f"fit_pca_iids\t{len(fit_pca)}\n")
        f.write(f"sample_qc_exclusion_iids\t{len(excluded_iids)}\n")
        f.write(f"classified_europeans_after_sample_qc\t{len(after_sample_qc)}\n")
        f.write(f"europeans_missing_codeable_ea\t{len(missing_ea)}\n")
        f.write(f"ea_candidates_below_min_age_at_basics\t{len(below_min_age)}\n")
        f.write(f"ea_age_candidates_missing_sex_covar\t{len(missing_sex)}\n")
        f.write(f"ea_age_sex_candidates_missing_fam\t{len(missing_fam)}\n")
        f.write(f"ea_age_sex_fam_candidates_missing_pcs\t{len(missing_pcs)}\n")
        f.write(f"eligible_classified_eur_samples\t{len(eligible_iids)}\n")
        f.write(f"eligible_fit_pca_oof_samples\t{len(train_iids)}\n")
        f.write(f"eligible_applied_6th_model_samples\t{len(applied_iids)}\n")
        f.write(f"final_kinship_holdout_kin0\t{args.final_kinship_holdout_kin0}\n")
        f.write(f"final_kinship_holdout_threshold\t{args.final_kinship_holdout_threshold:.10g}\n")
        f.write(f"final_kinship_edges_ge_threshold\t{final_total_edges_ge_threshold}\n")
        f.write(f"final_kinship_excluded_fit_pca_samples\t{len(final_excluded_iids)}\n")
        f.write(f"final_model_train_allowed_samples\t{len(final_train_iids)}\n")
        f.write(f"mean_yob_fit_pca\t{mean_yob:.10g}\n")
        f.write(f"outer_folds\t{args.outer_folds}\n")
        f.write(f"seed\t{args.seed}\n")
        f.write(f"feature_columns\t{len(feature_columns)}\n")
        f.write(f"feature_columns_sha256\t{feature_columns_sha256}\n")
        f.write(f"xgboost_saved_models\t{len(model_manifest_rows)}\n")
        f.write(f"lifestyle_branch_recode_rules_with_samples\t{len(branch_recode_rows)}\n")
        f.write(f"lifestyle_branch_recoded_sample_feature_cells\t{branch_recode_total}\n")
        f.write(f"pmi_missing_answer_rows\t{sum(feature_state['nonresponse'].values())}\n")
        f.write(f"prediction_raw_mean\t{pred_mean:.10g}\n")
        f.write(f"prediction_raw_sd\t{pred_sd:.10g}\n")
        f.write(f"oof_overall_pearson_proxy_teacher_z\t{pearson(proxy_z[train_idx], teacher_z[train_idx]):.10g}\n")
        f.write(f"oof_overall_spearman_proxy_teacher_z\t{spearman(proxy_z[train_idx], teacher_z[train_idx]):.10g}\n")
        if len(applied_idx):
            f.write(f"applied_pearson_proxy_teacher_z\t{pearson(proxy_z[applied_idx], teacher_z[applied_idx]):.10g}\n")
            f.write(f"applied_spearman_proxy_teacher_z\t{spearman(proxy_z[applied_idx], teacher_z[applied_idx]):.10g}\n")
        f.write("covar_cols\t" + ",".join(covar_cols) + "\n")

    runtime_manifest = {
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "xgboost_version": xgb.__version__,
        "params": params,
        "outer_folds": args.outer_folds,
        "seed": args.seed,
        "best_rounds": best_rounds,
        "final_rounds": final_rounds,
        "feature_columns": len(feature_columns),
        "feature_columns_sha256": feature_columns_sha256,
        "saved_xgboost_models": len(model_manifest_rows),
        "eligible_fit_pca_oof_samples": len(train_iids),
        "eligible_applied_samples": len(applied_iids),
        "final_kinship_holdout_kin0": args.final_kinship_holdout_kin0,
        "final_kinship_holdout_threshold": args.final_kinship_holdout_threshold,
        "final_kinship_edges_ge_threshold": final_total_edges_ge_threshold,
        "final_kinship_excluded_fit_pca_samples": len(final_excluded_iids),
        "final_model_train_allowed_samples": len(final_train_iids),
        "lifestyle_branch_recode_rules_with_samples": len(branch_recode_rows),
        "lifestyle_branch_recoded_sample_feature_cells": branch_recode_total,
        "pmi_missing_answer_rows": sum(feature_state["nonresponse"].values()),
    }
    with open(os.path.join(args.out_dir, "runtime_manifest.json"), "w") as f:
        json.dump(runtime_manifest, f, indent=2, sort_keys=True)

    with open(os.path.join(args.out_dir, "ses_ea_proxy_gwas_log.txt"), "w") as f:
        f.write("\n".join(log_lines) + "\n")

    print("\n=== Verification checks ===")
    passed = True

    def check(name, condition):
        nonlocal passed
        passed = passed and bool(condition)
        print(("PASS" if condition else "FAIL") + f": {name}")

    all_iids = [row["IID"] for row in all_score_rows]
    check("all score IIDs are classified European", set(all_iids) <= europeans)
    check("no score IID is in sample-QC exclusion list", not (set(all_iids) & excluded_iids))
    check("all score IIDs have confirmed genetic sex", set(all_iids) <= set(sex_map))
    check("all score IIDs have codeable EA and age >= threshold",
          all(iid in ea_rows and ea_rows[iid]["age_at_basics"] >= args.min_age_at_basics for iid in all_iids))
    check("all score IIDs have genotype FID", set(all_iids) <= set(fid_by_iid))
    check("all score IIDs have requested PCs", set(all_iids) <= set(pc_data))
    check("OOF and applied sets are disjoint", not (set(train_iids) & set(applied_iids)))
    check("final-model training set is subset of OOF fit_pca set", set(final_train_iids) <= set(train_iids))
    check("final-model training set excludes applied relatives", not (set(final_train_iids) & final_excluded_iids))
    check("all eligible fit_pca samples have OOF scores", np.all(np.isfinite(pred_raw[train_idx])))
    check("all eligible applied samples have sixth-model scores", np.all(np.isfinite(pred_raw[applied_idx])) if len(applied_idx) else True)
    check("phenotype rows match all score rows", len(all_score_rows) == len(eligible_iids))
    check("proxy_z is finite for every row", np.all(np.isfinite(proxy_z)))
    check("six XGBoost model files were saved",
          len(model_manifest_rows) == args.outer_folds + 1 and
          all(os.path.exists(os.path.join(args.out_dir, row["model_file"])) for row in model_manifest_rows))
    if not passed:
        sys.exit(1)
    print("All verification checks passed.")


if __name__ == "__main__":
    main()
