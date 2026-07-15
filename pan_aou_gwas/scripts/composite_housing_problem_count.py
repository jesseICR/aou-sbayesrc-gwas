#!/usr/bin/env python3
"""Approved complete-case count of housing problems selected in SDOH AHC-2."""

from __future__ import annotations

from typing import Mapping

from approved_composites import (
    CONTRADICTION,
    INVALID,
    MISSING,
    SCORED,
    CompositeDefinition,
    CompositeScore,
    normalize_answer,
    register_composite,
    resolve_answer,
    response_provenance,
    selected_response,
)


QUESTION_CONCEPT_ID = "40192402"
ITEM_CONCEPT = "ahc_2"
SOURCE_SURVEY = "Social Determinants of Health"
PHENOTYPE_ID = "comp_housing_problem_count"
CONSTRUCTION_ID = "housing_problem_count_complete_case_v1"

# Stable codebook option, OMOP answer concept ID, and display label.  The OMOP
# concept is the primary match; normalized labels/codes are checked fallbacks
# only when the answer-concept sidecar is absent.
PROBLEM_OPTIONS = (
    ("SDOH_44", "40192460", "Bug infestation"),
    ("SDOH_45", "40192479", "Mold"),
    ("SDOH_46", "40192393", "Lead paint or pipes"),
    ("SDOH_47", "40192434", "Inadequate heat"),
    ("SDOH_48", "40192495", "Oven or stove not working"),
    ("SDOH_49", "40192468", "No or not working smoke detector"),
    ("SDOH_50", "40192444", "Water leaks"),
)
NONE_OPTION = ("SDOH_51", "40192392", "None of the above")

_NONE = "none_of_the_above"
_NON_SUBSTANTIVE = "non_substantive"

ANSWER_CONCEPT_VALUES = {
    concept_id: answer_code for answer_code, concept_id, _label in PROBLEM_OPTIONS
}
ANSWER_CONCEPT_VALUES[NONE_OPTION[1]] = _NONE
# These standard AoU concepts are filtered by latest-valid ingest in production,
# but remain explicit here so manually supplied/audited records cannot become 0.
ANSWER_CONCEPT_VALUES.update(
    {
        "903079": _NON_SUBSTANTIVE,  # PMI: Prefer Not To Answer
        "1332892": _NON_SUBSTANTIVE,  # Prefer not to answer
        "903096": _NON_SUBSTANTIVE,  # PMI: Skip
    }
)


def _answer_aliases() -> dict[str, str]:
    aliases: dict[str, str] = {}
    for answer_code, _concept_id, label in PROBLEM_OPTIONS:
        aliases[normalize_answer(answer_code)] = answer_code
        aliases[normalize_answer(label)] = answer_code
    aliases[normalize_answer(NONE_OPTION[0])] = _NONE
    aliases[normalize_answer(NONE_OPTION[2])] = _NONE
    for label in (
        "Prefer not to answer",
        "PMI: Prefer Not To Answer",
        "Don't know",
        "Not sure",
        "Skip",
        "PMI: Skip",
    ):
        aliases[normalize_answer(label)] = _NON_SUBSTANTIVE
    return aliases


ANSWER_TEXT_VALUES = _answer_aliases()
PROBLEM_CODES = frozenset(answer_code for answer_code, _concept_id, _label in PROBLEM_OPTIONS)


def score_housing_problem_count(
    questions: Mapping[str, dict], iid: str
) -> CompositeScore:
    """Count seven selected problem types for one valid AHC-2 response event."""
    response = selected_response(questions, QUESTION_CONCEPT_ID, iid)
    provenance = response_provenance(response)
    source_provenance = (provenance,) if provenance is not None else ()
    if response is None:
        return CompositeScore(None, 0, MISSING, None, source_provenance)
    if not response.answers:
        return CompositeScore(None, 0, MISSING, response.age, source_provenance)

    selected: set[str] = set()
    for observation in response.answers:
        value, status = resolve_answer(
            observation, ANSWER_CONCEPT_VALUES, ANSWER_TEXT_VALUES
        )
        if status == INVALID:
            # A present but unknown concept or a known concept whose label maps
            # elsewhere is an auditable invalid record.  With no concept, an
            # unrecognized label supplies no valid checkbox selection.
            invalid_status = INVALID if observation.concept_id else MISSING
            return CompositeScore(
                None, 0, invalid_status, response.age, source_provenance
            )
        if value == _NON_SUBSTANTIVE:
            return CompositeScore(None, 0, MISSING, response.age, source_provenance)
        selected.add(str(value))

    selected_problems = selected & PROBLEM_CODES
    if _NONE in selected and selected_problems:
        return CompositeScore(None, 7, CONTRADICTION, response.age, source_provenance)
    if _NONE in selected:
        return CompositeScore(0, 7, SCORED, response.age, source_provenance)
    if selected_problems:
        return CompositeScore(
            len(selected_problems), 7, SCORED, response.age, source_provenance
        )
    return CompositeScore(None, 0, MISSING, response.age, source_provenance)


_ITEM_MAPPINGS = tuple(
    {
        "question_concept_id": QUESTION_CONCEPT_ID,
        "item_concept": ITEM_CONCEPT,
        "answer_code": answer_code,
        "answer_concept_id": answer_concept_id,
        "answer": label,
        "score": 1,
    }
    for answer_code, answer_concept_id, label in PROBLEM_OPTIONS
) + (
    {
        "question_concept_id": QUESTION_CONCEPT_ID,
        "item_concept": ITEM_CONCEPT,
        "answer_code": NONE_OPTION[0],
        "answer_concept_id": NONE_OPTION[1],
        "answer": NONE_OPTION[2],
        "score": 0,
        "zero_anchor": True,
    },
)

_DEFINITION_ANSWER_MAPPING = {
    f"{answer_concept_id} / {answer_code} / {label}": 1
    for answer_code, answer_concept_id, label in PROBLEM_OPTIONS
}
_DEFINITION_ANSWER_MAPPING.update(
    {
        f"{NONE_OPTION[1]} / {NONE_OPTION[0]} / {NONE_OPTION[2]}": 0,
        "None plus any scored problem": "invalid contradiction",
        "absent or unrecognized response": "missing",
    }
)


HOUSING_PROBLEM_COUNT_DEFINITION = register_composite(
    CompositeDefinition(
        phenotype_id=PHENOTYPE_ID,
        construction_id=CONSTRUCTION_ID,
        trait_type="composite",
        kind="quant",
        covar_mode="full",
        description=(
            "Pan-AoU-derived count of seven housing-condition problems selected "
            "in the SDOH AHC-2 checkbox question."
        ),
        source_surveys=(SOURCE_SURVEY,),
        source_qids=(QUESTION_CONCEPT_ID,),
        item_mappings=_ITEM_MAPPINGS,
        answer_mapping=_DEFINITION_ANSWER_MAPPING,
        missing_policy=(
            "One substantive checkbox response event is required. None alone is "
            "zero; no recognized option and non-substantive answers are missing; "
            "unknown concepts, concept/text conflicts, and None-plus-problem sets are "
            "invalid and unscored; no prorating."
        ),
        valid_range=(0, 7),
        scorer=score_housing_problem_count,
        interpretation=(
            "Higher values indicate a larger number of housing-condition problem "
            "types, not greater severity of any one problem."
        ),
        limitations=(
            "This is a pan-AoU-derived count, not a validated AHC total severity scale.",
            "The separate number-of-moves item and other housing questions are excluded.",
        ),
    )
)

# A compact conventional alias for callers that import each registration module.
DEFINITION = HOUSING_PROBLEM_COUNT_DEFINITION
