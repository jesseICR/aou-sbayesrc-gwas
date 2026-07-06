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
import json
import math
import re
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


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def norm_q(text: str) -> str:
    text = re.sub(r"<[^>]+>", " ", str(text or ""))
    text = re.sub(r"\s+", " ", text).strip().lower().rstrip("?.").strip()
    return text


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


def load_question_manifest(path: Path) -> dict[str, dict]:
    """Map normalized question text -> its manifest disposition row."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out.setdefault(norm_q(row["field_label"]), row)
    return out


def load_ordinal_lookup(path: Path) -> dict[tuple[str, str], float]:
    """(normalized question text, normalized answer label) -> ordinal value."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (norm_q(row["field_label"]), R.norm(row["answer_label"]))
            out[key] = float(row["ordinal_value"])
    return out


# --------------------------------------------------------------------------- #
# survey ingest: latest response per (person, question)
# --------------------------------------------------------------------------- #
def read_survey_rows(paths: list[Path], keep: set[str]):
    """Yield dict rows for retained samples from one or more survey CSVs."""
    for path in paths:
        if not path or not path.exists() or path.stat().st_size == 0:
            continue
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                pid = (row.get("person_id") or "").strip()
                if pid in keep:
                    yield row


def build_latest_responses(survey_paths, keep):
    """Return dict[qid] -> {question, pid -> (age, [(ans_text)...])} using latest datetime."""
    # per (pid, qid): keep max datetime, collect answers at that datetime
    latest = {}  # (pid, qid) -> [datetime, question_text, age, set(answer_text)]
    for row in read_survey_rows(survey_paths, keep):
        pid = row["person_id"].strip()
        qid = (row.get("question_concept_id") or "").strip()
        if not qid:
            continue
        dt = (row.get("survey_datetime") or "").strip()
        ans = (row.get("answer") or "").strip()
        try:
            age = float(row.get("age_at_survey") or "nan")
        except ValueError:
            age = float("nan")
        qtext = row.get("question") or ""
        k = (pid, qid)
        cur = latest.get(k)
        if cur is None or dt > cur[0]:
            latest[k] = [dt, qtext, age, {ans}]
        elif dt == cur[0]:
            cur[3].add(ans)

    questions = defaultdict(lambda: {"question": "", "responses": {}})
    for (pid, qid), (dt, qtext, age, answers) in latest.items():
        q = questions[qid]
        if not q["question"]:
            q["question"] = qtext
        q["responses"][pid] = (age, answers)
    return questions


# --------------------------------------------------------------------------- #
# phenotype builders -> {"kind": binary|quant, "values": {iid: (y, age)}}
# --------------------------------------------------------------------------- #
def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")[:40] or "x"


def build_survey_phenotypes(questions, qman, ord_lookup):
    """Yield (pheno_id, trait_type, kind, {iid: (y, age)}, meta)."""
    for qid, q in questions.items():
        qtext = q["question"]
        man = qman.get(norm_q(qtext))
        if man is None:
            continue  # question not in our included/classified manifest
        disp = man["disposition"]
        if disp.startswith("excluded") or disp == "numeric":
            # numeric handled separately below via value_as_number
            if disp != "numeric":
                continue
        responses = q["responses"]

        is_multi = man["phenotype_class"] == "multi_select"
        # ---- binary one-vs-rest per observed valid answer ------------------
        if disp in ("ordinal_and_binary", "binary_only", "nominal_binary", "flagged_review"):
            # collect the valid (non-missing) answer universe for this question
            valid_answers = set()
            for _, (_, answers) in responses.items():
                for a in answers:
                    if not R.is_missing(a):
                        valid_answers.add(a)
            for ans in sorted(valid_answers):
                values = {}
                for pid, (age, answers) in responses.items():
                    non_missing = {a for a in answers if not R.is_missing(a)}
                    if not non_missing:
                        continue
                    if ans in non_missing:
                        values[pid] = (1.0, age)
                    else:
                        # single-select control = answered another valid option;
                        # checkbox control = question shown and option not selected.
                        values[pid] = (0.0, age)
                pid_ = f"bin_{qid}__{slug(ans)}"
                yield pid_, "binary", "binary", values, {
                    "question_concept_id": qid,
                    "question": qtext,
                    "answer": ans,
                    "ordinal_rule": "",
                }
        # ---- ordinal -------------------------------------------------------
        if disp == "ordinal_and_binary":
            values = {}
            for pid, (age, answers) in responses.items():
                non_missing = [a for a in answers if not R.is_missing(a)]
                if len(non_missing) != 1:
                    continue
                v = ord_lookup.get((norm_q(qtext), R.norm(non_missing[0])))
                if v is not None:
                    values[pid] = (float(v), age)
            yield f"ord_{qid}", "ordinal", "quant", values, {
                "question_concept_id": qid,
                "question": qtext,
                "answer": "",
                "ordinal_rule": man["ordinal_rule"],
            }


def build_numeric_phenotypes(questions, qman):
    for qid, q in questions.items():
        man = qman.get(norm_q(q["question"]))
        if man is None or man["disposition"] != "numeric":
            continue
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
        yield f"num_{qid}", "numeric", "quant", values, {
            "question_concept_id": qid,
            "question": q["question"],
            "answer": "",
            "ordinal_rule": "",
        }


# Group -> PFHH category-screen question_concept_id used to recover controls.
PFHH_SCREEN_QID = {
    "Brain and nervous system": "43529272",
    "Mental health or substance use": "43529217",
    "Added skeletal/pain/injury": "702786",
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


def _pfhh_screen_completers(questions):
    out = {}
    for grp, sqid in PFHH_SCREEN_QID.items():
        pids = {}
        q = questions.get(sqid)
        if q:
            for pid, (age, answers) in q["responses"].items():
                if any(not R.is_missing(a) for a in answers):
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
            if any(R.norm(a) == "self" for a in answers):
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
            non_missing = [a for a in answers if not R.is_missing(a)]
            if not non_missing:
                burden.pop(pid, None)  # answered only PMI -> drop from denominator
                continue
            score = sum(PFHH_RELATION_WEIGHT.get(R.norm(a), 0.0) for a in non_missing)
            burden[pid] = (score, age)
        yield f"pfhh_burden_{cond}", "pfhh_sumscore", "quant", burden, {
            "question_concept_id": qid,
            "question": row.get("question", ""),
            "answer": "relatedness-weighted burden (self=1, 1st-deg=0.5, grandparent=0.25)",
            "ordinal_rule": "pfhh_relatedness_burden",
        }


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


# --------------------------------------------------------------------------- #
# residualize + write + PLINK2
# --------------------------------------------------------------------------- #
def prepare_and_write(pheno_id, kind, values, sex, pcs, outdir):
    """Return (pheno_path, keep_path, n, n_cases, n_controls) or None if it fails QC."""
    rows = []
    for iid, (y, age) in values.items():
        if iid not in sex or iid not in pcs or math.isnan(y) or (age is None) or math.isnan(age):
            continue
        rows.append((iid, y, age, sex[iid], pcs[iid]))

    if kind == "binary":
        ncase = sum(1 for _, y, *_ in rows if y == 1.0)
        nctrl = sum(1 for _, y, *_ in rows if y == 0.0)
        if ncase < MIN_CASES or nctrl < MIN_CONTROLS:
            return None
    else:
        if len(rows) < MIN_QUANT_N:
            return None
        levels = len({round(y, 6) for _, y, *_ in rows})
        if levels < MIN_ORDINAL_LEVELS:
            return None
        ncase = nctrl = 0

    iids = [r[0] for r in rows]
    y = np.array([r[1] for r in rows], dtype=float)
    age = np.array([r[2] for r in rows], dtype=float)
    sex_c = np.array([r[3] for r in rows], dtype=float) - 0.5
    pc = np.array([r[4] for r in rows], dtype=float)
    age_c = age - age.mean()
    covars = np.column_stack([age_c, sex_c, age_c * sex_c, pc])

    if kind == "binary":
        pheno_vec = residualize(y, covars)          # LPM residual, no INT
    else:
        pheno_vec = residualize(rint(y), covars)    # INT then residualize

    pdir = outdir / "phenotypes"
    pdir.mkdir(parents=True, exist_ok=True)
    name = f"{pheno_id}_resid"
    pheno_path = pdir / f"{pheno_id}.resid.pheno.tsv"
    keep_path = pdir / f"{pheno_id}.keep.tsv"
    with open(pheno_path, "w") as f:
        f.write(f"FID\tIID\t{name}\n")
        for iid, v in zip(iids, pheno_vec):
            f.write(f"{iid}\t{iid}\t{v:.17g}\n")
    with open(keep_path, "w") as f:
        for iid in iids:
            f.write(f"{iid}\t{iid}\n")
    return pheno_path, keep_path, name, len(rows), ncase, nctrl


def run_plink2(plink2, bfile, keep_path, pheno_path, pheno_name, out_prefix, force):
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    expected = out_prefix.parent / f"{out_prefix.name}.{pheno_name}.glm.linear"
    if expected.exists() and expected.stat().st_size > 0 and not force:
        return expected, 0.0
    cmd = [
        plink2,
        "--bfile", str(bfile),
        "--keep", str(keep_path),
        "--pheno", str(pheno_path),
        "--pheno-name", pheno_name,
        "--glm", "allow-no-covars",
        "--no-input-missing-phenotype",
        "--out", str(out_prefix),
    ]
    start = time.time()
    res = subprocess.run(cmd, text=True, capture_output=True)
    elapsed = time.time() - start
    (out_prefix.parent / f"{out_prefix.name}.plink2.log").write_text(res.stderr + "\n" + res.stdout)
    if res.returncode != 0 or not expected.exists():
        raise RuntimeError(f"PLINK2 failed for {out_prefix}; see log.")
    return expected, elapsed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bfile", type=Path, required=True)
    ap.add_argument("--keep", type=Path, required=True)
    ap.add_argument("--sex", type=Path, required=True)
    ap.add_argument("--pcs", type=Path, required=True)
    ap.add_argument("--survey-csv", type=Path, required=True)
    ap.add_argument("--bhp-csv", type=Path, default=None)
    ap.add_argument("--measurements-csv", type=Path, default=None)
    ap.add_argument("--question-manifest", type=Path, required=True)
    ap.add_argument("--ordinal-manifest", type=Path, required=True)
    ap.add_argument("--pfhh-allowlist", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--phenotypes", default="", help="comma-separated pheno_id filter (smoke test).")
    ap.add_argument("--plink2-bin", default=shutil.which("plink2") or "plink2")
    ap.add_argument("--skip-gwas", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    keep = set(load_keep(args.keep))
    sex = load_sex(args.sex)
    pcs = load_pcs(args.pcs)
    log(f"keep={len(keep)}  sex={len(sex)}  pcs={len(pcs)}")

    qman = load_question_manifest(args.question_manifest)
    ord_lookup = load_ordinal_lookup(args.ordinal_manifest)
    log(f"manifest questions={len(qman)}  ordinal answer maps={len(ord_lookup)}")

    survey_paths = [args.survey_csv]
    if args.bhp_csv:
        survey_paths.append(args.bhp_csv)
    log("Building latest-response table ...")
    questions = build_latest_responses(survey_paths, keep)
    log(f"questions with responses: {len(questions)}")

    only = {p.strip() for p in args.phenotypes.split(",") if p.strip()}

    builders = [
        build_survey_phenotypes(questions, qman, ord_lookup),
        build_numeric_phenotypes(questions, qman),
        build_pfhh_phenotypes(questions, args.pfhh_allowlist),
        build_measurement_phenotypes(args.measurements_csv, keep),
    ]

    manifest_rows = []
    metadir = args.outdir / "metadata"
    metadir.mkdir(parents=True, exist_ok=True)

    for gen in builders:
        for pheno_id, trait_type, kind, values, meta in gen:
            if only and pheno_id not in only:
                continue
            prep = prepare_and_write(pheno_id, kind, values, sex, pcs, args.outdir)
            if prep is None:
                continue
            pheno_path, keep_path, pheno_name, n, ncase, nctrl = prep
            row = {
                "pheno_id": pheno_id,
                "trait_type": trait_type,
                "kind": kind,
                "n": n,
                "n_cases": ncase,
                "n_controls": nctrl,
                "ordinal_rule": meta.get("ordinal_rule", ""),
                "question_concept_id": meta.get("question_concept_id", ""),
                "question": meta.get("question", ""),
                "answer": meta.get("answer", ""),
                "pheno_path": str(pheno_path),
            }
            if not args.skip_gwas:
                out_prefix = args.outdir / "gwas" / pheno_id / pheno_id
                glm, elapsed = run_plink2(
                    args.plink2_bin, args.bfile, keep_path, pheno_path, pheno_name, out_prefix, args.force
                )
                row["glm"] = str(glm)
                row["gwas_seconds"] = round(elapsed, 1)
            manifest_rows.append(row)
            if len(manifest_rows) % 100 == 0:
                log(f"  {len(manifest_rows)} phenotypes done")

    man_path = metadir / "phenotype_manifest.tsv"
    if manifest_rows:
        cols = list(manifest_rows[0].keys())
        with open(man_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(manifest_rows)
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
                "skip_gwas": bool(args.skip_gwas),
            },
            indent=2,
        )
        + "\n"
    )
    log(f"Wrote {man_path} ({len(manifest_rows)} phenotypes passed QC)")


if __name__ == "__main__":
    main()
