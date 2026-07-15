"""Registry and builders for the approved sex-stratified survey GWAS.

The values in this module are constructed before the pipeline applies the
genetic-sex filter.  Each definition therefore carries an explicit
``sex_filter`` and uses age plus PC1--PC10 covariates (``covar_mode=agepc``).
Missing/non-substantive survey responses are never converted to controls.
"""

from __future__ import annotations

import json
import math
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import ordinal_rules as R


GENDER_IDENTITY_QID = "1585838"
CLOSER_GENDER_QID = "1585348"
SEXUAL_ORIENTATION_QID = "1585899"
ACTIVE_DUTY_QID = "1585852"
MARITAL_PRIMARY_QID = "1585892"
MARITAL_COPE_QID = "1332833"
UNDER18_PRIMARY_QID = "1585890"
UNDER18_COPE_QID = "1333023"


@dataclass(frozen=True)
class SexSpecificDefinition:
    phenotype_id: str
    construction_id: str
    trait_type: str
    kind: str
    sex_filter: str
    source_qids: tuple[str, ...]


DEFINITIONS = (
    SexSpecificDefinition(
        "bin_gender_transgender_expanded_male",
        "gender_transgender_or_trans_woman_male_v1",
        "binary",
        "binary",
        "male",
        (GENDER_IDENTITY_QID, CLOSER_GENDER_QID),
    ),
    SexSpecificDefinition(
        "bin_gender_transgender_expanded_female",
        "gender_transgender_or_trans_man_female_v1",
        "binary",
        "binary",
        "female",
        (GENDER_IDENTITY_QID, CLOSER_GENDER_QID),
    ),
    SexSpecificDefinition(
        "bin_thebasics_sexualorientation__gay_male",
        "sexual_orientation_gay_male_v1",
        "binary",
        "binary",
        "male",
        (SEXUAL_ORIENTATION_QID,),
    ),
    SexSpecificDefinition(
        "bin_thebasics_sexualorientation__gay_or_lesbian_female",
        "sexual_orientation_gay_or_lesbian_female_v1",
        "binary",
        "binary",
        "female",
        (SEXUAL_ORIENTATION_QID,),
    ),
    SexSpecificDefinition(
        "bin_activeduty_activedutyservestatus__yes_male",
        "active_duty_yes_male_v1",
        "binary",
        "binary",
        "male",
        (ACTIVE_DUTY_QID,),
    ),
    SexSpecificDefinition(
        "bin_activeduty_activedutyservestatus__yes_female",
        "active_duty_yes_female_v1",
        "binary",
        "binary",
        "female",
        (ACTIVE_DUTY_QID,),
    ),
    SexSpecificDefinition(
        "num_livingsituation_peopleunder18_male",
        "people_under_18_basics_cope_pooled_male_v1",
        "numeric",
        "quant",
        "male",
        (UNDER18_PRIMARY_QID, UNDER18_COPE_QID),
    ),
    SexSpecificDefinition(
        "num_livingsituation_peopleunder18_female",
        "people_under_18_basics_cope_pooled_female_v1",
        "numeric",
        "quant",
        "female",
        (UNDER18_PRIMARY_QID, UNDER18_COPE_QID),
    ),
    SexSpecificDefinition(
        "bin_maritalstatus_currentmaritalstatus__divorced_male",
        "marital_status_divorced_basics_cope_pooled_male_v1",
        "binary",
        "binary",
        "male",
        (MARITAL_PRIMARY_QID, MARITAL_COPE_QID),
    ),
    SexSpecificDefinition(
        "bin_maritalstatus_currentmaritalstatus__divorced_female",
        "marital_status_divorced_basics_cope_pooled_female_v1",
        "binary",
        "binary",
        "female",
        (MARITAL_PRIMARY_QID, MARITAL_COPE_QID),
    ),
    SexSpecificDefinition(
        "bin_maritalstatus_currentmaritalstatus__never_married_male",
        "marital_status_never_married_basics_cope_pooled_male_v1",
        "binary",
        "binary",
        "male",
        (MARITAL_PRIMARY_QID, MARITAL_COPE_QID),
    ),
    SexSpecificDefinition(
        "bin_maritalstatus_currentmaritalstatus__never_married_female",
        "marital_status_never_married_basics_cope_pooled_female_v1",
        "binary",
        "binary",
        "female",
        (MARITAL_PRIMARY_QID, MARITAL_COPE_QID),
    ),
)

DEFINITION_BY_ID = {definition.phenotype_id: definition for definition in DEFINITIONS}
SOURCE_QIDS = frozenset(qid for definition in DEFINITIONS for qid in definition.source_qids)


def phenotype_ids() -> set[str]:
    return set(DEFINITION_BY_ID)


def construction_ids_by_phenotype() -> dict[str, str]:
    return {
        definition.phenotype_id: definition.construction_id
        for definition in DEFINITIONS
    }


def _norm(text: str) -> str:
    return R.norm(text).replace("’", "'")


def _finite_age(age) -> bool:
    try:
        return age is not None and math.isfinite(float(age))
    except (TypeError, ValueError):
        return False


TEXT_ALIASES = {
    GENDER_IDENTITY_QID: {
        "man": "man",
        "woman": "woman",
        "non-binary": "non_binary",
        "non binary": "non_binary",
        "transgender": "transgender",
        "none of these describe me, and i'd like to consider additional options": "additional_options",
        "none of these describe me, and i would like to consider additional options": "additional_options",
    },
    CLOSER_GENDER_QID: {
        "trans man/transgender man/ftm": "trans_man",
        "trans man": "trans_man",
        "transgender man": "trans_man",
        "ftm": "trans_man",
        "trans woman/transgender woman/mtf": "trans_woman",
        "trans woman": "trans_woman",
        "transgender woman": "trans_woman",
        "mtf": "trans_woman",
        "genderqueer": "genderqueer",
        "genderfluid": "genderfluid",
        "gender variant": "gender_variant",
        "two spirit": "two_spirit",
        "two-spirit": "two_spirit",
        "questioning or unsure of your gender identity": "unsure",
        "unsure": "unsure",
        "none of these describe me, and i want to specify": "specified_gender",
        "specified gender": "specified_gender",
    },
    SEXUAL_ORIENTATION_QID: {
        "gay": "gay",
        "lesbian": "lesbian",
        "bisexual": "bisexual",
        "straight": "straight",
        "straight; that is, not gay or lesbian, etc": "straight",
        "none": "none",
        "none of these describe me and i'd like to see additional options": "none",
        "none of these describe me and i would like to see additional options": "none",
    },
    ACTIVE_DUTY_QID: {"yes": "yes", "no": "no"},
    MARITAL_PRIMARY_QID: {
        "married": "married",
        "divorced": "divorced",
        "widowed": "widowed",
        "separated": "separated",
        "never married": "never_married",
        "living with partner": "living_with_partner",
    },
    MARITAL_COPE_QID: {
        "married": "married",
        "divorced": "divorced",
        "widowed": "widowed",
        "separated": "separated",
        "never married": "never_married",
        "living with partner": "living_with_partner",
    },
}

# Stable AoU answer concepts used when present.  Text is still checked for a
# disagreement; releases without a known answer concept use the unique aliases.
CONCEPT_ALIASES = {
    GENDER_IDENTITY_QID: {
        "1585839": "man",
        "1585840": "woman",
        "1585841": "non_binary",
        "1585842": "transgender",
        "1585843": "additional_options",
    },
    CLOSER_GENDER_QID: {
        "1585349": "trans_man",
        "1585350": "trans_woman",
        "1585351": "genderqueer",
        "1585352": "genderfluid",
        "1585353": "gender_variant",
        "701374": "two_spirit",
        "1585354": "unsure",
        "1585355": "specified_gender",
    },
    SEXUAL_ORIENTATION_QID: {
        "1585900": "straight",
        "1585901": "gay",
        "1585902": "lesbian",
        "1585903": "bisexual",
        "1585904": "none",
    },
    ACTIVE_DUTY_QID: {"1585853": "yes", "1585854": "no"},
    MARITAL_PRIMARY_QID: {
        "1585893": "married",
        "1585894": "divorced",
        "1585895": "widowed",
        "1585896": "separated",
        "1585897": "never_married",
        "1585898": "living_with_partner",
    },
}


def _text_alias(qid: str, text: str) -> str | None:
    normalized = _norm(text)
    aliases = TEXT_ALIASES.get(qid, {})
    value = aliases.get(normalized)
    if value is None and qid == SEXUAL_ORIENTATION_QID and normalized.startswith("straight;"):
        value = "straight"
    return value


def _canonical_response(questions: dict, qid: str, participant_id: str):
    """Return (set of canonical substantive answers, age), or None if invalid."""
    question = questions.get(qid)
    if not question:
        return None
    response = question.get("responses", {}).get(participant_id)
    if response is None:
        return None
    age, answers = response
    if not _finite_age(age) or not answers:
        return None
    sidecar = question.get("response_sidecars", {}).get(participant_id, {})
    concept_ids = tuple(sidecar.get("answer_concept_ids", ()))
    if concept_ids and len(concept_ids) != len(answers):
        return None

    canonical = set()
    for index, answer in enumerate(answers):
        if R.is_missing(answer):
            return None
        text_value = _text_alias(qid, answer)
        concept_id = concept_ids[index].strip() if concept_ids else ""
        concept_value = CONCEPT_ALIASES.get(qid, {}).get(concept_id)
        if concept_value is not None and text_value is not None and concept_value != text_value:
            return None
        value = concept_value if concept_value is not None else text_value
        if value is None:
            return None
        canonical.add(value)
    return canonical, float(age)


def _numeric_response(questions: dict, qid: str, participant_id: str):
    question = questions.get(qid)
    if not question:
        return None
    response = question.get("responses", {}).get(participant_id)
    if response is None:
        return None
    age, answers = response
    substantive = [answer for answer in answers if not R.is_missing(answer)]
    if len(substantive) != 1 or not _finite_age(age):
        return None
    try:
        value = float(substantive[0])
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value) or value < 1.0 or value > 20.0 or value != int(value):
        return None
    return value, float(age)


def _participant_ids(questions: dict, *qids: str) -> set[str]:
    return {
        participant_id
        for qid in qids
        for participant_id in questions.get(qid, {}).get("responses", {})
    }


def _base_meta(definition: SexSpecificDefinition, *, item: str, question: str, answer: str):
    return {
        "question_concept_id": "|".join(definition.source_qids),
        "item_concept": item,
        "question": question,
        "answer": answer,
        "ordinal_rule": "",
        "covar_mode": "agepc",
        "sex_filter": definition.sex_filter,
        "construction_id": definition.construction_id,
    }


def _build_transgender(questions: dict, definition: SexSpecificDefinition, target: str):
    values = {}
    for participant_id in _participant_ids(questions, GENDER_IDENTITY_QID, CLOSER_GENDER_QID):
        parent = _canonical_response(questions, GENDER_IDENTITY_QID, participant_id)
        followup = _canonical_response(questions, CLOSER_GENDER_QID, participant_id)
        if parent is None and followup is None:
            continue
        parent_answers = parent[0] if parent else set()
        followup_answers = followup[0] if followup else set()
        ages = [got[1] for got in (parent, followup) if got is not None]
        is_case = "transgender" in parent_answers or target in followup_answers
        values[participant_id] = (float(is_case), statistics.fmean(ages))
    return values


def _build_checkbox_binary(questions: dict, qid: str, case_answers: set[str]):
    values = {}
    for participant_id in _participant_ids(questions, qid):
        got = _canonical_response(questions, qid, participant_id)
        if got is None:
            continue
        answers, age = got
        values[participant_id] = (float(bool(answers & case_answers)), age)
    return values


def _pooled_categorical(questions: dict, primary_qid: str, cope_qid: str):
    values = {}
    for participant_id in _participant_ids(questions, primary_qid, cope_qid):
        primary = _canonical_response(questions, primary_qid, participant_id)
        if primary is not None and len(primary[0]) == 1:
            values[participant_id] = (next(iter(primary[0])), primary[1], 0.0)
            continue
        cope = _canonical_response(questions, cope_qid, participant_id)
        if cope is not None and len(cope[0]) == 1:
            values[participant_id] = (next(iter(cope[0])), cope[1], 1.0)
    return values


def _pooled_numeric(questions: dict, primary_qid: str, cope_qid: str):
    values = {}
    for participant_id in _participant_ids(questions, primary_qid, cope_qid):
        primary = _numeric_response(questions, primary_qid, participant_id)
        if primary is not None:
            values[participant_id] = (primary[0], primary[1], 0.0)
            continue
        cope = _numeric_response(questions, cope_qid, participant_id)
        if cope is not None:
            values[participant_id] = (cope[0], cope[1], 1.0)
    return values


def build_sex_specific_phenotypes(questions: dict):
    """Yield all 12 collision-reserved sex-stratified phenotypes."""
    for phenotype_id, target in (
        ("bin_gender_transgender_expanded_male", "trans_woman"),
        ("bin_gender_transgender_expanded_female", "trans_man"),
    ):
        definition = DEFINITION_BY_ID[phenotype_id]
        values = _build_transgender(questions, definition, target)
        yield phenotype_id, "binary", "binary", values, _base_meta(
            definition,
            item="gender_transgender_expanded",
            question="Transgender identity, expanded with sex-specific closer-gender description",
            answer=(
                "Transgender OR Trans woman/Transgender Woman/MTF"
                if definition.sex_filter == "male"
                else "Transgender OR Trans man/Transgender Man/FTM"
            ),
        )

    for phenotype_id, case_answers in (
        ("bin_thebasics_sexualorientation__gay_male", {"gay"}),
        ("bin_thebasics_sexualorientation__gay_or_lesbian_female", {"gay", "lesbian"}),
    ):
        definition = DEFINITION_BY_ID[phenotype_id]
        values = _build_checkbox_binary(questions, SEXUAL_ORIENTATION_QID, case_answers)
        yield phenotype_id, "binary", "binary", values, _base_meta(
            definition,
            item="thebasics_sexualorientation",
            question="Sexual orientation, sex-stratified",
            answer=" OR ".join(sorted(case_answers)),
        )

    active_values = _build_checkbox_binary(questions, ACTIVE_DUTY_QID, {"yes"})
    for phenotype_id in (
        "bin_activeduty_activedutyservestatus__yes_male",
        "bin_activeduty_activedutyservestatus__yes_female",
    ):
        definition = DEFINITION_BY_ID[phenotype_id]
        yield phenotype_id, "binary", "binary", active_values, _base_meta(
            definition,
            item="activeduty_activedutyservestatus",
            question="Ever served on active duty, sex-stratified",
            answer="Yes; No is control; missing/PNA excluded",
        )

    under18 = _pooled_numeric(questions, UNDER18_PRIMARY_QID, UNDER18_COPE_QID)
    under18_values = {pid: (value, age) for pid, (value, age, _source) in under18.items()}
    under18_source = {pid: source for pid, (_value, _age, source) in under18.items()}
    for phenotype_id in (
        "num_livingsituation_peopleunder18_male",
        "num_livingsituation_peopleunder18_female",
    ):
        definition = DEFINITION_BY_ID[phenotype_id]
        meta = _base_meta(
            definition,
            item="livingsituation_peopleunder18",
            question="People living at home under age 18, Basics-primary/COPE-fill-in",
            answer="Integer 1..20; missing/out-of-range excluded",
        )
        meta.update({
            "extra_covariates": {"from_cope": under18_source},
            "extra_covariates_label": "from_cope",
        })
        yield phenotype_id, "numeric", "quant", under18_values, meta

    marital = _pooled_categorical(questions, MARITAL_PRIMARY_QID, MARITAL_COPE_QID)
    marital_source = {pid: source for pid, (_answer, _age, source) in marital.items()}
    for status in ("divorced", "never_married"):
        for sex_filter in ("male", "female"):
            phenotype_id = (
                f"bin_maritalstatus_currentmaritalstatus__{status}_{sex_filter}"
            )
            definition = DEFINITION_BY_ID[phenotype_id]
            values = {
                pid: (float(answer == status), age)
                for pid, (answer, age, _source) in marital.items()
            }
            meta = _base_meta(
                definition,
                item="maritalstatus_currentmaritalstatus",
                question="Current marital status, Basics-primary/COPE-fill-in and sex-stratified",
                answer=f"{status.replace('_', ' ').title()} vs all other substantive statuses",
            )
            meta.update({
                "extra_covariates": {"from_cope": marital_source},
                "extra_covariates_label": "from_cope",
            })
            yield phenotype_id, "binary", "binary", values, meta


def write_construction_qc(
    path: Path | str,
    *,
    phenotype_id: str,
    values: Mapping[str, tuple[float, float]],
    metadata: Mapping[str, object],
    codebook_fingerprints: Mapping[str, str],
    git_state: Mapping[str, object],
) -> Path:
    definition = DEFINITION_BY_ID[phenotype_id]
    if metadata.get("construction_id") != definition.construction_id:
        raise ValueError(f"Construction mismatch for {phenotype_id}")
    raw = [float(value) for value, _age in values.values() if math.isfinite(float(value))]
    counts = Counter(raw)
    payload = {
        "phenotype_id": phenotype_id,
        "construction_id": definition.construction_id,
        "trait_type": definition.trait_type,
        "kind": definition.kind,
        "sex_filter": definition.sex_filter,
        "covar_mode": "agepc",
        "question_concept_id": metadata.get("question_concept_id", ""),
        "item_concept": metadata.get("item_concept", ""),
        "answer_and_control_definition": metadata.get("answer", ""),
        "missing_policy": "No response, skip, PNA, don't know, invalid, and unrecognized answers are missing",
        "raw_moments": {
            "minimum": min(raw) if raw else None,
            "maximum": max(raw) if raw else None,
            "mean": statistics.fmean(raw) if raw else None,
            "standard_deviation": statistics.pstdev(raw) if raw else None,
        },
        "raw_histogram": {format(value, ".17g"): count for value, count in sorted(counts.items())},
        "prefilter_scored_n": len(raw),
        "final_filtered_n": None,
        "codebook_fingerprints_sha256": dict(codebook_fingerprints),
        "pipeline_git": dict(git_state),
    }
    qc_path = Path(path)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return qc_path


def finalize_construction_qc(path: Path | str | None, final_filtered_n: int) -> None:
    if not path:
        return
    qc_path = Path(path)
    if not qc_path.is_file():
        return
    payload = json.loads(qc_path.read_text())
    payload["final_filtered_n"] = int(final_filtered_n)
    qc_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
