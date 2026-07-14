#!/usr/bin/env python3
"""Export pan-AoU GWAS results as compact row-aligned Parquet files."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import time
from pathlib import Path
from numbers import Integral, Real

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


DEFAULT_PLINK2 = "/opt/workbench-tools/binaries/bin/plink2"
DEFAULT_SNP_COUNT = 1_140_557
SURVEY_ABBREVIATIONS = {
    "Behavioral Health and Personality": "BHP",
    "COVID-19 Participant Experience (COPE)": "COPE",
    "Emotional Health History and Well-Being": "EHW",
    "Family Health History": "FHH",
    "Healthcare Access and Utilization": "HCAU",
    "Life Functioning": "LifeFunctioning",
    "Lifestyle": "Lifestyle",
    "Minute Survey on COVID-19 Vaccines": "COVIDVax",
    "Overall Health": "OverallHealth",
    "Personal Medical History": "PMH",
    "Personal and Family Health History": "PFHH",
    "Social Determinants of Health": "SDOH",
    "The Basics": "Basics",
}

POOLED_SOURCE_DESCRIPTORS = {
    "phq_gad_ehhwb_cope_pooled_v1": (
        "EHW+COPE",
        "Emotional Health History and Well-Being + COVID-19 Participant Experience (COPE)",
    ),
    "pss_sdoh_cope_pooled_v1": (
        "SDOH+COPE",
        "Social Determinants of Health + COVID-19 Participant Experience (COPE)",
    ),
    "mos_ss_sdoh_cope_pooled_v1": (
        "SDOH+COPE",
        "Social Determinants of Health + COVID-19 Participant Experience (COPE)",
    ),
    "ucla_sdoh_cope_pooled_v1": (
        "SDOH+COPE",
        "Social Determinants of Health + COVID-19 Participant Experience (COPE)",
    ),
    "eds_sdoh_cope_harmonized_4level_v1": (
        "SDOH+COPE",
        "Social Determinants of Health + COVID-19 Participant Experience (COPE)",
    ),
    "auditc_lifestyle_cope_pooled_v1": (
        "Lifestyle+COPE",
        "Lifestyle + COVID-19 Participant Experience (COPE)",
    ),
    "auditc_population_zero_pooled_v1": (
        "Lifestyle+COPE",
        "Lifestyle + COVID-19 Participant Experience (COPE)",
    ),
    "subjective_wellbeing_ehw_cope_pooled_v1": (
        "EHW+COPE",
        "Emotional Health History and Well-Being + COVID-19 Participant Experience (COPE)",
    ),
}


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean_cell(row.get(field, "")) for field in fields})


def clean_cell(value: object) -> object:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    if isinstance(value, Integral) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, Real) and not isinstance(value, bool):
        if math.isfinite(float(value)) and float(value).is_integer():
            return int(value)
    return value


def norm_label(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().lower().split())


def parse_options(options: object) -> list[tuple[str, str]]:
    if options is None:
        return []
    try:
        if pd.isna(options):
            return []
    except (TypeError, ValueError):
        pass
    parsed: list[tuple[str, str]] = []
    for raw in str(options).split(" | "):
        raw = raw.strip()
        if not raw:
            continue
        if "=" in raw:
            code, label = raw.split("=", 1)
            parsed.append((code.strip(), label.strip()))
        else:
            parsed.append(("", raw))
    return parsed


def source_abbrev(source_name: object, pheno_id: str = "") -> str:
    source = "" if source_name is None or pd.isna(source_name) else str(source_name)
    if source in SURVEY_ABBREVIATIONS:
        return SURVEY_ABBREVIATIONS[source]
    if pheno_id.startswith("cog_"):
        return "EAProxy_ETM"
    if pheno_id.startswith("zip3_"):
        return "ZIP3_SES"
    if pheno_id in {"height_cm", "bmi_kg_m2", "systolic_bp_mmhg", "diastolic_bp_mmhg",
                    "mean_arterial_pressure_mmhg", "pulse_pressure_mmhg", "heart_rate_bpm"}:
        return "PhysicalMeasurements"
    if pheno_id.startswith("pfhh_burden_"):
        return "PFHH"
    if pheno_id.startswith("comp_"):
        return "Composite"
    return source or "Derived"


def construction_method(row: pd.Series, qrow: dict[str, object] | None) -> str:
    pheno_id = str(row.get("pheno_id", ""))
    trait_type = str(row.get("trait_type", ""))
    if pheno_id.startswith("bin_") or trait_type == "binary":
        return "binary_one_vs_rest"
    if pheno_id.startswith("ord_") or trait_type == "ordinal":
        return "ordinal_numeric"
    if pheno_id.startswith("num_"):
        return "numeric_text"
    if pheno_id.startswith("comp_"):
        return "composite_score"
    if pheno_id.startswith("pfhh_burden_"):
        return "pfhh_family_burden_score"
    if pheno_id.startswith("cog_"):
        return "external_score"
    if pheno_id.startswith("zip3_"):
        return "zip3_ses"
    if qrow and qrow.get("field_type") == "slider":
        return "slider_numeric"
    if trait_type in {"numeric", "quant"} or str(row.get("kind", "")) == "quant":
        return "continuous"
    return trait_type or "unknown"


def selection_semantics(row: pd.Series, qrow: dict[str, object] | None) -> str:
    method = construction_method(row, qrow)
    field_type = str(qrow.get("field_type", "")) if qrow else ""
    pheno_class = str(qrow.get("phenotype_class", "")) if qrow else ""
    if method == "binary_one_vs_rest":
        if field_type == "checkbox" or pheno_class == "multi_select":
            return "multi_select_select_all_that_apply"
        if field_type in {"radio", "dropdown"} or pheno_class == "single_select":
            return "single_select_one_response"
        return "binary_target_vs_eligible_non_target"
    if method == "ordinal_numeric":
        return "ordinal_response_coded_numeric"
    if method in {"numeric_text", "slider_numeric", "continuous"}:
        return "continuous_or_numeric_response"
    if method == "composite_score":
        return "derived_multi_item_composite"
    if method == "pfhh_family_burden_score":
        return "derived_family_history_burden"
    if method == "external_score":
        return "precomputed_external_score"
    if method == "zip3_ses":
        return "zip3_contextual_ses"
    return method


def case_definition(row: pd.Series, semantics: str) -> str:
    answer = clean_cell(row.get("answer", ""))
    if semantics == "single_select_one_response":
        return f"selected target answer: {answer}"
    if semantics == "multi_select_select_all_that_apply":
        return f"selected target answer: {answer}"
    if semantics == "binary_target_vs_eligible_non_target":
        return f"target answer present: {answer}"
    return ""


def control_definition(row: pd.Series, semantics: str) -> str:
    if semantics == "single_select_one_response":
        return "selected a different valid answer for the same question"
    if semantics == "multi_select_select_all_that_apply":
        return "eligible respondent did not select the target answer"
    if semantics == "binary_target_vs_eligible_non_target":
        return "eligible respondent did not have the target answer"
    return ""


def load_question_maps(question_manifest: Path) -> tuple[dict[str, dict[str, object]], dict[tuple[str, str], dict[str, object]]]:
    if not question_manifest.exists():
        return {}, {}
    q = read_tsv(question_manifest)
    by_item: dict[str, dict[str, object]] = {}
    by_item_qid: dict[tuple[str, str], dict[str, object]] = {}
    for _, row in q.iterrows():
        item = str(row.get("item_concept", "") or "")
        qid = str(row.get("question_concept_id", "") or "")
        record = {k: clean_cell(v) for k, v in row.items()}
        if item and item not in by_item:
            by_item[item] = record
        if item and qid:
            by_item_qid[(item, qid)] = record
    return by_item, by_item_qid


def find_question_row(
    row: pd.Series,
    by_item: dict[str, dict[str, object]],
    by_item_qid: dict[tuple[str, str], dict[str, object]],
) -> dict[str, object] | None:
    item = str(row.get("item_concept", "") or "")
    qid = str(row.get("question_concept_id", "") or "")
    if item and qid and (item, qid) in by_item_qid:
        return by_item_qid[(item, qid)]
    if item and item in by_item:
        return by_item[item]
    return None


def pipe_join(values: list[object]) -> str:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        value = clean_cell(value)
        if value == "":
            continue
        text = str(value)
        if text not in seen:
            seen.add(text)
            out.append(text)
    return " | ".join(out)


def covariate_mode_description(row: pd.Series) -> str:
    mode = clean_cell(row.get("covar_mode", "")) or "full"
    if mode == "sexpc":
        covars = "sex_c + PC1..PC10"
    elif mode == "agepc":
        covars = "age_c + PC1..PC10"
    else:
        covars = "age_c + sex_c + age_c:sex_c + PC1..PC10"
    extra = clean_cell(row.get("extra_covariates", ""))
    if extra:
        covars = f"{covars} + {extra}"
    return f"pre-residualized on {covars}; PLINK2 --glm allow-no-covars"


def sample_set_description(row: pd.Series) -> str:
    sex_filter = clean_cell(row.get("sex_filter", "")) or "all"
    if sex_filter == "female":
        return "unrelated European female-only pan-AoU GWAS keep-list"
    if sex_filter == "male":
        return "unrelated European male-only pan-AoU GWAS keep-list"
    return "unrelated European pan-AoU GWAS keep-list"


def build_metadata(
    manifest: pd.DataFrame,
    skipped: pd.DataFrame,
    question_manifest: Path,
    item_inventory: Path,
    ordinal_mapping: Path,
    composite_items: Path,
    external_scores: Path,
    out_dir: Path,
) -> None:
    by_item, by_item_qid = load_question_maps(question_manifest)
    inventory = read_tsv(item_inventory) if item_inventory.exists() else pd.DataFrame()
    ordinal = read_tsv(ordinal_mapping) if ordinal_mapping.exists() else pd.DataFrame()
    composite = read_tsv(composite_items) if composite_items.exists() else pd.DataFrame()
    external = read_tsv(external_scores, comment="#") if external_scores.exists() else pd.DataFrame()

    options_by_item: dict[str, list[tuple[str, str]]] = {}
    if not inventory.empty:
        inventory_by_item: dict[str, dict[str, object]] = {}
        for _, row in inventory.iterrows():
            item = str(row.get("item_concept", "") or "")
            if item and item not in inventory_by_item:
                inventory_by_item[item] = {k: clean_cell(v) for k, v in row.items()}
            opts = parse_options(row.get("options", ""))
            if item and opts and item not in options_by_item:
                options_by_item[item] = opts
    else:
        inventory_by_item = {}

    ordinal_by_item_label: dict[tuple[str, str], object] = {}
    ordinal_summary_by_item: dict[str, list[str]] = {}
    if not ordinal.empty:
        for _, row in ordinal.iterrows():
            item = str(row.get("item_concept", "") or "")
            label = str(row.get("answer_label", "") or "")
            value = clean_cell(row.get("ordinal_value", ""))
            if not item or label == "":
                continue
            ordinal_by_item_label[(item, norm_label(label))] = value
            ordinal_summary_by_item.setdefault(item, []).append(f"{label}={value}")

    binary_run: dict[tuple[str, str], list[str]] = {}
    binary_skip: dict[tuple[str, str], list[str]] = {}
    ordinal_run_by_item: dict[str, list[str]] = {}
    for _, row in manifest.iterrows():
        item = str(row.get("item_concept", "") or "")
        answer = norm_label(row.get("answer", ""))
        pheno_id = str(row.get("pheno_id", ""))
        if pheno_id.startswith("bin_") and item and answer:
            binary_run.setdefault((item, answer), []).append(pheno_id)
        if pheno_id.startswith("ord_") and item:
            ordinal_run_by_item.setdefault(item, []).append(pheno_id)
    for _, row in skipped.iterrows():
        item = str(row.get("item_concept", "") or "")
        answer = norm_label(row.get("answer", ""))
        pheno_id = str(row.get("pheno_id", ""))
        reason = str(row.get("skip_reason", ""))
        if pheno_id.startswith("bin_") and item and answer:
            binary_skip.setdefault((item, answer), []).append(f"{pheno_id}:{reason}")

    catalog_rows: list[dict[str, object]] = []
    for _, row in manifest.iterrows():
        qrow = find_question_row(row, by_item, by_item_qid)
        item = str(row.get("item_concept", "") or "")
        invrow = inventory_by_item.get(item, {})
        opts = options_by_item.get(item, [])
        valid_options = [label for _, label in opts]
        method = construction_method(row, qrow)
        semantics = selection_semantics(row, qrow)
        construction_id = clean_cell(row.get("construction_id", ""))
        pooled_source = POOLED_SOURCE_DESCRIPTORS.get(str(construction_id))
        survey_name = qrow.get("survey", "") if qrow else invrow.get("survey", "")
        pheno_id = str(row["pheno_id"])
        if pooled_source:
            source, survey_name = pooled_source
        else:
            source = source_abbrev(survey_name, pheno_id)
        catalog_rows.append({
            "pheno_id": pheno_id,
            "trait_label": row.get("question", "") or pheno_id,
            "trait_type": row.get("trait_type", ""),
            "kind": row.get("kind", ""),
            "n": row.get("n", ""),
            "n_cases": row.get("n_cases", ""),
            "n_controls": row.get("n_controls", ""),
            "sex_filter": row.get("sex_filter", "all") or "all",
            "extra_covariates": row.get("extra_covariates", ""),
            "construction_id": construction_id,
            "source_abbrev": source,
            "source_name": survey_name,
            "item_concept": item,
            "question_concept_id": row.get("question_concept_id", ""),
            "question": row.get("question", ""),
            "target_answer": row.get("answer", ""),
            "field_type": qrow.get("field_type", "") if qrow else invrow.get("field_type", ""),
            "phenotype_class": qrow.get("phenotype_class", "") if qrow else invrow.get("phenotype_class", ""),
            "construction_method": method,
            "selection_semantics": semantics,
            "case_definition": case_definition(row, semantics),
            "control_definition": control_definition(row, semantics),
            "valid_response_options": pipe_join(valid_options),
            "excluded_response_options": "",
            "ordinal_rule": row.get("ordinal_rule", ""),
            "ordinal_coding_summary": pipe_join(ordinal_summary_by_item.get(item, [])),
            "branching_logic_present": invrow.get("has_branching", ""),
            "transform": "IRNT/residualized phenotype as written by pan_aou_gwas",
            "covariate_mode": covariate_mode_description(row),
            "genotype_panel": "HapMap3 HQ bfile",
            "sample_set": sample_set_description(row),
            "source_glm": row.get("glm", ""),
            "source_sumstats": row.get("sumstats", ""),
            "export_parquet": str(out_dir / "gwas" / f"{pheno_id}.parquet"),
        })

    catalog_fields = [
        "pheno_id", "trait_label", "trait_type", "kind", "n", "n_cases", "n_controls",
        "sex_filter", "extra_covariates", "construction_id",
        "source_abbrev", "source_name", "item_concept", "question_concept_id", "question",
        "target_answer", "field_type", "phenotype_class", "construction_method",
        "selection_semantics", "case_definition", "control_definition",
        "valid_response_options", "excluded_response_options", "ordinal_rule",
        "ordinal_coding_summary", "branching_logic_present", "transform", "covariate_mode",
        "genotype_panel", "sample_set", "source_glm", "source_sumstats", "export_parquet",
    ]
    write_tsv(out_dir / "phenotype_catalog.tsv", catalog_rows, catalog_fields)

    response_rows: list[dict[str, object]] = []
    if not inventory.empty:
        for _, row in inventory.iterrows():
            item = str(row.get("item_concept", "") or "")
            survey_name = row.get("survey", "")
            qid = clean_cell(row.get("question_concept_id", ""))
            question = clean_cell(row.get("field_label", ""))
            for answer_code, answer_label in parse_options(row.get("options", "")):
                key = (item, norm_label(answer_label))
                response_rows.append({
                    "source_abbrev": source_abbrev(survey_name),
                    "source_name": survey_name,
                    "item_concept": item,
                    "question_concept_id": qid,
                    "field_type": row.get("field_type", ""),
                    "phenotype_class": row.get("phenotype_class", ""),
                    "question": question,
                    "answer_code": answer_code,
                    "answer_label": answer_label,
                    "binary_pheno_ids_run": pipe_join(binary_run.get(key, [])),
                    "binary_skip_reasons": pipe_join(binary_skip.get(key, [])),
                    "ordinal_pheno_ids_run": pipe_join(ordinal_run_by_item.get(item, [])),
                    "ordinal_value_if_applicable": ordinal_by_item_label.get(key, ""),
                    "selection_semantics": (
                        "multi_select_select_all_that_apply"
                        if row.get("field_type") == "checkbox" or row.get("phenotype_class") == "multi_select"
                        else "single_select_one_response"
                        if row.get("field_type") in {"radio", "dropdown"} or row.get("phenotype_class") == "single_select"
                        else ""
                    ),
                })

    response_fields = [
        "source_abbrev", "source_name", "item_concept", "question_concept_id", "field_type",
        "phenotype_class", "question", "answer_code", "answer_label", "binary_pheno_ids_run",
        "binary_skip_reasons", "ordinal_pheno_ids_run", "ordinal_value_if_applicable",
        "selection_semantics",
    ]
    write_tsv(out_dir / "response_options.tsv", response_rows, response_fields)

    ordinal_rows: list[dict[str, object]] = []
    if not ordinal.empty:
        for _, row in ordinal.iterrows():
            ordinal_rows.append({k: row.get(k, "") for k in ordinal.columns})
    write_tsv(out_dir / "ordinal_coding.tsv", ordinal_rows, list(ordinal.columns) if not ordinal.empty else [
        "survey", "item_concept", "question_concept_id", "ordinal_rule", "ordinal_source",
        "confidence", "answer_label", "answer_label_normalized", "ordinal_value", "field_label",
    ])

    survey_names = set(SURVEY_ABBREVIATIONS)
    if not inventory.empty:
        survey_names.update(str(x) for x in inventory["survey"].dropna().unique())
    source_rows = [{"source_abbrev": abbr, "source_name": name} for name, abbr in sorted(SURVEY_ABBREVIATIONS.items())]
    source_rows.extend([
        {"source_abbrev": "EAProxy_ETM", "source_name": "Precomputed EA-proxy and Exploring-the-Mind cognitive scores"},
        {"source_abbrev": "ZIP3_SES", "source_name": "AoU ZIP3 socioeconomic context table"},
        {"source_abbrev": "PhysicalMeasurements", "source_name": "AoU physical measurement table"},
        {"source_abbrev": "Composite", "source_name": "Derived multi-item survey composite"},
    ])
    source_rows.extend(
        {"source_abbrev": abbrev, "source_name": name}
        for abbrev, name in sorted(set(POOLED_SOURCE_DESCRIPTORS.values()))
    )
    seen_sources = {(r["source_abbrev"], r["source_name"]) for r in source_rows}
    for name in sorted(survey_names):
        row = {"source_abbrev": source_abbrev(name), "source_name": name}
        key = (row["source_abbrev"], row["source_name"])
        if key not in seen_sources:
            seen_sources.add(key)
            source_rows.append(row)
    write_tsv(out_dir / "survey_sources.tsv", source_rows, ["source_abbrev", "source_name"])

    composite_path = out_dir / "composite_items.tsv"
    if not composite.empty:
        composite.to_csv(composite_path, sep="\t", index=False)
    external_path = out_dir / "external_scores.tsv"
    if not external.empty:
        external.to_csv(external_path, sep="\t", index=False)


def read_bim(bfile: Path) -> pd.DataFrame:
    bim = Path(f"{bfile}.bim")
    if not bim.exists():
        raise FileNotFoundError(f"missing BIM: {bim}")
    return pd.read_csv(
        bim,
        sep="\t",
        header=None,
        names=["chrom", "rsid", "cm", "pos", "A1", "A2"],
        dtype={"chrom": str, "rsid": str, "pos": np.int64, "A1": str, "A2": str},
    )


def load_sbayes_subset(sbayesrc: Path, rsids: set[str]) -> pd.DataFrame:
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        sbayesrc,
        usecols=["chrom", "pos", "ref", "alt", "rsid"],
        dtype={"chrom": str, "pos": np.int64, "ref": str, "alt": str, "rsid": str},
        chunksize=500_000,
    ):
        sub = chunk[chunk["rsid"].isin(rsids)]
        if len(sub):
            chunks.append(sub.copy())
    if not chunks:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "rsid"])
    return pd.concat(chunks, ignore_index=True)


def normalize_keep_for_bfile(bfile: Path, keep: Path, out_dir: Path, force: bool) -> Path:
    """Return a keep file whose FID/IID pairs match the target bfile.

    Some older pan-AoU workdirs contain keep files written as IID/IID, while
    the HapMap3 bfile uses FID=0.  PLINK2 requires exact FID/IID matches, so
    normalize to the bfile's own FID for every requested IID.
    """
    out = out_dir / "reference_freq" / "hapmap3_unrelated_eur.keep"
    if out.exists() and out.stat().st_size > 0 and not force:
        return out
    fam = Path(f"{bfile}.fam")
    if not fam.exists():
        raise FileNotFoundError(f"missing FAM: {fam}")
    fid_by_iid: dict[str, str] = {}
    with open(fam) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                fid_by_iid[parts[1]] = parts[0]
    requested: list[str] = []
    with open(keep) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            requested.append(parts[1] if len(parts) >= 2 else parts[0])
    missing = [iid for iid in requested if iid not in fid_by_iid]
    if missing:
        examples = ",".join(missing[:10])
        raise ValueError(f"{len(missing)} keep-list IIDs are absent from {fam}; examples: {examples}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for iid in requested:
            writer.writerow([fid_by_iid[iid], iid])
    return out


def run_plink_freq(plink2: str, bfile: Path, keep: Path, out_dir: Path, force: bool) -> Path:
    prefix = out_dir / "reference_freq" / "hapmap3_unrelated_eur"
    afreq = prefix.with_suffix(".afreq")
    if afreq.exists() and afreq.stat().st_size > 0 and not force:
        return afreq
    prefix.parent.mkdir(parents=True, exist_ok=True)
    normalized_keep = normalize_keep_for_bfile(bfile, keep, out_dir, force)
    cmd = [
        plink2,
        "--bfile", str(bfile),
        "--keep", str(normalized_keep),
        "--nonfounders",
        "--freq", "cols=chrom,pos,ref,alt1,alt1freq,nobs",
        "--out", str(prefix),
    ]
    start = time.time()
    res = subprocess.run(cmd, text=True, capture_output=True)
    (prefix.with_suffix(".freq.log")).write_text(res.stderr + "\n" + res.stdout)
    if res.returncode != 0:
        raise RuntimeError(f"PLINK2 --freq failed; see {prefix.with_suffix('.freq.log')}")
    if not afreq.exists() or afreq.stat().st_size == 0:
        raise RuntimeError(f"PLINK2 --freq did not create {afreq}")
    print(f"computed reference frequencies in {time.time() - start:.1f}s: {afreq}", flush=True)
    return afreq


def load_a1_freq(afreq: Path, reference: pd.DataFrame) -> pd.Series:
    freq = pd.read_csv(afreq, sep="\t", dtype={"#CHROM": str, "CHROM": str, "ID": str})
    id_col = "ID"
    ref_col = "REF" if "REF" in freq.columns else None
    alt_col = "ALT1" if "ALT1" in freq.columns else "ALT" if "ALT" in freq.columns else None
    freq_col = (
        "ALT1_FREQ" if "ALT1_FREQ" in freq.columns else
        "ALT_FREQS" if "ALT_FREQS" in freq.columns else
        "ALT_FREQ" if "ALT_FREQ" in freq.columns else
        None
    )
    if id_col not in freq.columns or alt_col is None or freq_col is None:
        raise ValueError(f"cannot identify ID/ALT/frequency columns in {afreq}: {freq.columns.tolist()}")
    freq = freq[[id_col, ref_col, alt_col, freq_col] if ref_col else [id_col, alt_col, freq_col]].copy()
    freq.rename(columns={id_col: "rsid", ref_col or alt_col: "freq_ref", alt_col: "freq_alt", freq_col: "freq_value"}, inplace=True)
    merged = reference[["rsid", "A1", "A2"]].merge(freq, on="rsid", how="left", sort=False)
    if merged["freq_value"].isna().any():
        raise ValueError(f"{int(merged['freq_value'].isna().sum())} reference SNPs missing from {afreq}")
    freq_value = pd.to_numeric(merged["freq_value"], errors="coerce")
    if freq_value.isna().any():
        raise ValueError(f"{int(freq_value.isna().sum())} nonnumeric frequencies in {afreq}")
    if "freq_alt" not in merged.columns:
        raise ValueError(f"frequency file lacks an alternate allele column: {afreq}")
    alt = merged["freq_alt"].astype(str)
    a1 = merged["A1"].astype(str)
    a2 = merged["A2"].astype(str)
    if "freq_ref" in merged.columns:
        ref = merged["freq_ref"].astype(str)
        a1_freq = np.where(alt.eq(a1), freq_value, np.where(ref.eq(a1), 1.0 - freq_value, np.nan))
        bad = np.isnan(a1_freq)
        if bad.any():
            raise ValueError(f"{int(bad.sum())} frequency rows do not match canonical alleles")
        return pd.Series(a1_freq.astype(np.float32), name="A1_freq")
    a1_freq = np.where(alt.eq(a1), freq_value, np.where(alt.eq(a2), 1.0 - freq_value, np.nan))
    bad = np.isnan(a1_freq)
    if bad.any():
        raise ValueError(f"{int(bad.sum())} frequency rows do not match canonical alleles")
    return pd.Series(a1_freq.astype(np.float32), name="A1_freq")


def build_reference(args: argparse.Namespace) -> pd.DataFrame:
    bim = read_bim(args.bfile)
    if bim["rsid"].duplicated().any():
        raise ValueError("HapMap3 BIM contains duplicate rsids")
    sbayes = load_sbayes_subset(args.sbayesrc, set(bim["rsid"]))
    if sbayes["rsid"].duplicated().any():
        dupes = sbayes.loc[sbayes["rsid"].duplicated(), "rsid"].head().tolist()
        raise ValueError(f"SBayesRC support file contains duplicate HapMap3 rsids: {dupes}")
    merged = bim.merge(sbayes, on="rsid", how="left", sort=False, indicator=True)
    missing = merged["_merge"].eq("left_only")
    if missing.any():
        examples = ",".join(merged.loc[missing, "rsid"].head(10).tolist())
        raise ValueError(f"{int(missing.sum())} HapMap3 final SNPs missing from SBayesRC support file; examples: {examples}")
    chrom_ok = merged["chrom_x"].astype(str).eq(merged["chrom_y"].astype(str))
    pos_ok = merged["pos_x"].eq(merged["pos_y"])
    alt_ok = merged["A1"].eq(merged["alt"])
    ref_ok = merged["A2"].eq(merged["ref"])
    if not (chrom_ok.all() and pos_ok.all() and alt_ok.all() and ref_ok.all()):
        raise ValueError(
            "HapMap3/SBayesRC validation failed: "
            f"chrom_mismatch={int((~chrom_ok).sum())} pos_mismatch={int((~pos_ok).sum())} "
            f"A1_not_alt={int((~alt_ok).sum())} A2_not_ref={int((~ref_ok).sum())}"
        )

    afreq = args.afreq if args.afreq else run_plink_freq(args.plink2_bin, args.bfile, args.keep, args.out_dir, args.force_freq)
    reference = pd.DataFrame({
        "rsid": merged["rsid"].astype(str),
        "chrom": merged["chrom_x"].astype(str),
        "pos": merged["pos_x"].astype(np.int32),
        "A1": merged["A1"].astype(str),
        "A2": merged["A2"].astype(str),
    })
    reference["A1_freq"] = load_a1_freq(afreq, reference)
    return reference


def write_reference_parquet(reference: pd.DataFrame, out_path: Path) -> None:
    schema = pa.schema([
        ("rsid", pa.string()),
        ("chrom", pa.string()),
        ("pos", pa.int32()),
        ("A1", pa.string()),
        ("A2", pa.string()),
        ("A1_freq", pa.float32()),
    ])
    table = pa.Table.from_pandas(reference, schema=schema, preserve_index=False)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    pq.write_table(table, tmp, compression="zstd", compression_level=19, row_group_size=len(reference))
    os.replace(tmp, out_path)


def source_path(row: pd.Series) -> Path:
    sumstats = Path(str(row.get("sumstats", "")))
    if sumstats.exists() and sumstats.stat().st_size > 0:
        return sumstats
    glm = Path(str(row.get("glm", "")))
    if glm.exists() and glm.stat().st_size > 0:
        return glm
    raise FileNotFoundError(f"missing source sumstats/glm for {row.get('pheno_id', '')}")


def read_source_sumstats(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    if "rsid" in header:
        cols = ["rsid", "chrom", "pos", "allele1", "n", "beta", "se"]
        rename = {}
    else:
        cols = ["ID", "POS", "A1", "OBS_CT", "BETA", "SE"]
        chrom_col = "#CHROM" if "#CHROM" in header else "CHROM"
        cols.insert(0, chrom_col)
        rename = {chrom_col: "chrom", "ID": "rsid", "POS": "pos", "A1": "allele1",
                  "OBS_CT": "n", "BETA": "beta", "SE": "se"}
    df = pd.read_csv(path, sep="\t", usecols=cols, dtype={"rsid": str, "ID": str, "chrom": str, "#CHROM": str, "CHROM": str, "allele1": str, "A1": str})
    if rename:
        df = df.rename(columns=rename)
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(np.int64)
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    df["n"] = pd.to_numeric(df["n"], errors="coerce").astype("Int32")
    return df[["rsid", "chrom", "pos", "allele1", "beta", "se", "n"]]


def parquet_complete(path: Path, expected_rows: int) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        meta = pq.read_metadata(path)
    except Exception:
        return False
    if meta.num_rows != expected_rows:
        return False
    names = meta.schema.names
    return names == ["beta", "se", "N"]


def parquet_stored_metrics(path: Path) -> dict[str, object]:
    meta = pq.read_metadata(path)
    kv = meta.metadata or {}
    out: dict[str, object] = {
        "rows": meta.num_rows,
        "n_flipped_beta_rows": "",
        "n_allele_match_failures": 0,
        "n_beta_missing": "",
        "n_se_missing": "",
        "n_N_missing": "",
        "export_size_bytes": path.stat().st_size,
    }
    for key in ["n_flipped_beta_rows", "n_allele_match_failures", "n_beta_missing", "n_se_missing", "n_N_missing"]:
        raw = kv.get(key.encode())
        if raw is not None:
            out[key] = int(raw.decode())
    return out


def write_gwas_parquet(df: pd.DataFrame, reference: pd.DataFrame, out_path: Path) -> dict[str, object]:
    n_ref = len(reference)
    if len(df) != n_ref:
        raise ValueError(f"row count mismatch: source={len(df)} reference={n_ref}")
    order_ok = (
        df["rsid"].to_numpy(dtype=str) == reference["rsid"].to_numpy(dtype=str)
    ).all() and (
        df["chrom"].to_numpy(dtype=str) == reference["chrom"].to_numpy(dtype=str)
    ).all() and (
        df["pos"].to_numpy(dtype=np.int64) == reference["pos"].to_numpy(dtype=np.int64)
    ).all()
    if not order_ok:
        raise ValueError("source SNP order does not match reference rsid/chrom/pos order")

    src_a1 = df["allele1"].to_numpy(dtype=str)
    ref_a1 = reference["A1"].to_numpy(dtype=str)
    ref_a2 = reference["A2"].to_numpy(dtype=str)
    same = src_a1 == ref_a1
    flip = src_a1 == ref_a2
    fail = ~(same | flip)
    if fail.any():
        bad_idx = np.flatnonzero(fail)[:10].tolist()
        bad = ",".join(f"{reference.iloc[i]['rsid']}:{src_a1[i]}" for i in bad_idx)
        raise ValueError(f"{int(fail.sum())} allele matching failures; examples: {bad}")

    beta = df["beta"].to_numpy(dtype=np.float64)
    beta = np.where(flip, -beta, beta).astype(np.float32)
    se = df["se"].to_numpy(dtype=np.float32)
    n = df["n"].astype("Int32")
    table = pa.Table.from_arrays(
        [
            pa.array(beta, type=pa.float32()),
            pa.array(se, type=pa.float32()),
            pa.Array.from_pandas(n, type=pa.int32()),
        ],
        names=["beta", "se", "N"],
    )
    metrics = {
        "rows": n_ref,
        "n_flipped_beta_rows": int(flip.sum()),
        "n_allele_match_failures": 0,
        "n_beta_missing": int(np.isnan(beta).sum()),
        "n_se_missing": int(np.isnan(se).sum()),
        "n_N_missing": int(n.isna().sum()),
    }
    table = table.replace_schema_metadata({
        key.encode(): str(value).encode()
        for key, value in metrics.items()
        if key != "rows"
    })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    pq.write_table(
        table,
        tmp,
        compression="zstd",
        compression_level=19,
        use_byte_stream_split=["beta", "se"],
        use_dictionary=["N"],
        row_group_size=n_ref,
    )
    os.replace(tmp, out_path)
    metrics["export_size_bytes"] = out_path.stat().st_size
    return metrics


def select_manifest_rows(manifest: pd.DataFrame, pheno_ids: list[str], max_phenotypes: int | None) -> pd.DataFrame:
    selected = manifest
    if pheno_ids:
        wanted = set(pheno_ids)
        selected = selected[selected["pheno_id"].isin(wanted)].copy()
        missing = wanted - set(selected["pheno_id"])
        if missing:
            raise ValueError(f"requested pheno_id(s) absent from manifest: {sorted(missing)}")
    if max_phenotypes is not None:
        selected = selected.head(max_phenotypes).copy()
    return selected.reset_index(drop=True)


def export_gwas_files(args: argparse.Namespace, manifest: pd.DataFrame, reference: pd.DataFrame) -> list[dict[str, object]]:
    qc_rows: list[dict[str, object]] = []
    gwas_dir = args.out_dir / "gwas"
    total = len(manifest)
    for i, row in manifest.iterrows():
        pheno_id = str(row["pheno_id"])
        out_path = gwas_dir / f"{pheno_id}.parquet"
        source = source_path(row)
        if not args.force and parquet_complete(out_path, len(reference)):
            metrics = parquet_stored_metrics(out_path)
            qc_rows.append({
                "pheno_id": pheno_id,
                "status": "skipped_existing",
                "source": str(source),
                "export_parquet": str(out_path),
                **metrics,
                "elapsed_sec": 0,
            })
            continue
        start = time.time()
        df = read_source_sumstats(source)
        metrics = write_gwas_parquet(df, reference, out_path)
        elapsed = time.time() - start
        qc_rows.append({
            "pheno_id": pheno_id,
            "status": "exported",
            "source": str(source),
            "export_parquet": str(out_path),
            **metrics,
            "elapsed_sec": round(elapsed, 3),
        })
        if (i + 1) % args.progress_every == 0 or i + 1 == total:
            print(f"exported {i + 1}/{total}: {pheno_id} ({elapsed:.1f}s, flips={metrics['n_flipped_beta_rows']})", flush=True)
    return qc_rows


def write_run_info(args: argparse.Namespace, manifest: pd.DataFrame, reference: pd.DataFrame) -> None:
    info = {
        "created_at_unix": time.time(),
        "manifest": str(args.manifest),
        "bfile": str(args.bfile),
        "keep": str(args.keep),
        "sbayesrc": str(args.sbayesrc),
        "afreq": str(args.afreq) if args.afreq else "",
        "n_phenotypes_selected": int(len(manifest)),
        "n_reference_snps": int(len(reference)),
        "canonical_a1": "SBayesRC_hg38_alt_and_HapMap3_BIM_allele1",
        "canonical_a2": "SBayesRC_hg38_ref_and_HapMap3_BIM_allele2",
        "per_gwas_columns": {"beta": "float32", "se": "float32", "N": "int32"},
    }
    (args.out_dir / "export_run.json").write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--skipped", type=Path, default=None)
    ap.add_argument("--bfile", type=Path, required=True)
    ap.add_argument("--keep", type=Path, required=True)
    ap.add_argument("--sbayesrc", type=Path, default=repo_root / "data/support/sbayesrc_hg38.csv")
    ap.add_argument("--afreq", type=Path, default=None, help="Precomputed PLINK2 .afreq over the GWAS keep-list.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--plink2-bin", default=DEFAULT_PLINK2)
    ap.add_argument("--question-manifest", type=Path, default=repo_root / "pan_aou_gwas/metadata/survey_question_manifest.tsv")
    ap.add_argument("--item-inventory", type=Path, default=repo_root / "pan_aou_gwas/metadata/survey_item_inventory.tsv")
    ap.add_argument("--ordinal-mapping", type=Path, default=repo_root / "pan_aou_gwas/metadata/ordinal_mapping_manifest.tsv")
    ap.add_argument("--composite-items", type=Path, default=repo_root / "pan_aou_gwas/metadata/composite_items_manifest.tsv")
    ap.add_argument("--external-scores", type=Path, default=repo_root / "pan_aou_gwas/metadata/external_scores.tsv")
    ap.add_argument("--pheno-id", action="append", default=[], help="Restrict export to one phenotype ID; repeatable.")
    ap.add_argument("--max-phenotypes", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--force-freq", action="store_true")
    ap.add_argument("--metadata-only", action="store_true")
    ap.add_argument("--progress-every", type=int, default=25)
    args = ap.parse_args()
    if args.max_phenotypes is not None and args.max_phenotypes < 1:
        raise SystemExit("--max-phenotypes must be >= 1")
    if args.progress_every < 1:
        raise SystemExit("--progress-every must be >= 1")
    if args.skipped is None:
        args.skipped = args.manifest.parent / "skipped_phenotypes.tsv"
    return args


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = read_tsv(args.manifest)
    manifest = select_manifest_rows(manifest, args.pheno_id, args.max_phenotypes)
    skipped = read_tsv(args.skipped) if args.skipped.exists() else pd.DataFrame()

    print(f"building reference for {args.bfile}", flush=True)
    reference = build_reference(args)
    if len(reference) != DEFAULT_SNP_COUNT:
        print(f"WARNING: reference SNP count is {len(reference)}, expected {DEFAULT_SNP_COUNT}", flush=True)
    write_reference_parquet(reference, args.out_dir / "snp_reference.parquet")
    print(f"wrote reference rows={len(reference)}", flush=True)

    build_metadata(
        manifest, skipped, args.question_manifest, args.item_inventory, args.ordinal_mapping,
        args.composite_items, args.external_scores, args.out_dir,
    )
    write_run_info(args, manifest, reference)

    qc_fields = [
        "pheno_id", "status", "source", "export_parquet", "rows", "n_flipped_beta_rows",
        "n_allele_match_failures", "n_beta_missing", "n_se_missing", "n_N_missing",
        "export_size_bytes", "elapsed_sec",
    ]
    if args.metadata_only:
        write_tsv(args.out_dir / "export_qc.tsv", [], qc_fields)
        print("metadata-only export complete", flush=True)
        return
    qc_rows = export_gwas_files(args, manifest, reference)
    write_tsv(args.out_dir / "export_qc.tsv", qc_rows, qc_fields)
    total_bytes = sum(int(r["export_size_bytes"]) for r in qc_rows if clean_cell(r.get("export_size_bytes", "")) != "")
    print(f"done: phenotypes={len(qc_rows)} total_gwas_bytes={total_bytes}", flush=True)


if __name__ == "__main__":
    main()
