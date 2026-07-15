#!/usr/bin/env python3
"""Approved complete-case count of four lifetime psychotic experiences."""

from __future__ import annotations

from typing import Mapping

from approved_composites import (
    INCOMPLETE,
    INVALID,
    MISSING,
    SCORED,
    CompositeDefinition,
    CompositeScore,
    finite_mean,
    register_composite,
    resolve_answer,
    response_provenance,
    selected_response,
)


PHENOTYPE_ID = "psych_psychotic_experiences_count"
CONSTRUCTION_ID = "psychotic_experiences_count_four_item_complete_case_v1"
SOURCE_SURVEY = "Behavioral Health and Personality"

# The visual item is intentionally included.  The pre-existing three-item
# psych_psychotic_experiences_any construction is outside this definition and
# remains unchanged.
PSYCHOTIC_EXPERIENCE_ITEMS = (
    {
        "qid": "1703885",
        "item_concept": "mhqukb_49",
        "experience_type": "visual experience",
        "question": (
            "Did you ever see something that wasn't really there that other people "
            "could not see?"
        ),
    },
    {
        "qid": "1703901",
        "item_concept": "cidi5_21",
        "experience_type": "auditory experience",
        "question": "Did you ever hear voices that other people did not hear?",
    },
    {
        "qid": "1703915",
        "item_concept": "cidi5_22",
        "experience_type": "special signs or signals",
        "question": (
            "Did you ever think some force was trying to communicate directly with "
            "you by sending special signs or signals?"
        ),
    },
    {
        "qid": "1703871",
        "item_concept": "cidi5_23",
        "experience_type": "persecutory belief",
        "question": (
            "Did you ever believe there was a plot to harm you or to have people "
            "follow you, but your friends or family did not think this was true?"
        ),
    },
)

PSYCHOTIC_EXPERIENCE_QIDS = tuple(
    str(item["qid"]) for item in PSYCHOTIC_EXPERIENCE_ITEMS
)

# These Yes/No answer concepts are shared by all four live survey items.
_NON_SUBSTANTIVE = "missing"
ANSWER_CONCEPT_VALUES: Mapping[str, int | str] = {
    "1703969": 1,  # Yes
    "1703958": 0,  # No
    "903079": _NON_SUBSTANTIVE,  # PMI: Prefer Not To Answer
    "1332892": _NON_SUBSTANTIVE,  # Prefer not to answer
    "903096": _NON_SUBSTANTIVE,  # PMI: Skip
}

ANSWER_TEXT_VALUES: Mapping[str, int | str] = {
    "yes": 1,
    "no": 0,
    "prefer not to answer": _NON_SUBSTANTIVE,
    "pmi: prefer not to answer": _NON_SUBSTANTIVE,
    "pmi_prefernottoanswer": _NON_SUBSTANTIVE,
    "pmi: skip": _NON_SUBSTANTIVE,
    "skip": _NON_SUBSTANTIVE,
    "don't know": _NON_SUBSTANTIVE,
    "dont know": _NON_SUBSTANTIVE,
    "not sure": _NON_SUBSTANTIVE,
}


def score_psychotic_experiences_count(
    questions: Mapping[str, dict], iid: str
) -> CompositeScore:
    """Count four valid Yes/No lifetime-experience endorsements.

    The approved primary score is complete case.  A known answer concept is
    resolved first and normalized text is used only when no concept is present.
    Any unknown concept, concept/text disagreement, or multiple-answer radio
    response invalidates the composite.  Non-substantive and absent components
    remain missing and are never interpreted as No.
    """
    values: list[int] = []
    ages: list[float | None] = []
    provenance = []
    has_invalid = False
    has_missing = False

    for qid in PSYCHOTIC_EXPERIENCE_QIDS:
        response = selected_response(questions, qid, iid)
        component_provenance = response_provenance(response)
        if component_provenance is not None:
            provenance.append(component_provenance)

        if response is None or not response.answers:
            has_missing = True
            continue
        if len(response.answers) != 1:
            has_invalid = True
            continue

        value, status = resolve_answer(
            response.answers[0], ANSWER_CONCEPT_VALUES, ANSWER_TEXT_VALUES
        )
        if status == INVALID:
            has_invalid = True
            continue
        if value == _NON_SUBSTANTIVE:
            has_missing = True
            continue
        if not isinstance(value, int) or value not in (0, 1):
            has_invalid = True
            continue

        values.append(value)
        ages.append(response.age)

    observed_count = len(values)
    if has_invalid:
        status = INVALID
        raw_score = None
    elif observed_count == len(PSYCHOTIC_EXPERIENCE_QIDS) and not has_missing:
        status = SCORED
        raw_score = sum(values)
    elif observed_count:
        status = INCOMPLETE
        raw_score = None
    else:
        status = MISSING
        raw_score = None

    return CompositeScore(
        raw_score=raw_score,
        observed_component_count=observed_count,
        status=status,
        age=finite_mean(ages),
        source_provenance=tuple(provenance),
    )


PSYCHOTIC_EXPERIENCES_COUNT_DEFINITION = register_composite(
    CompositeDefinition(
        phenotype_id=PHENOTYPE_ID,
        construction_id=CONSTRUCTION_ID,
        trait_type="derived_psych",
        kind="quant",
        covar_mode="full",
        description=(
            "Complete-case count of four distinct types of lifetime unusual or "
            "psychotic experiences"
        ),
        source_surveys=(SOURCE_SURVEY,),
        source_qids=PSYCHOTIC_EXPERIENCE_QIDS,
        item_mappings=PSYCHOTIC_EXPERIENCE_ITEMS,
        answer_mapping={
            "1703958 / No": 0,
            "1703969 / Yes": 1,
            "Prefer not to answer, skip, or don't know": "missing",
        },
        missing_policy=(
            "Complete case: require one valid substantive Yes/No response for all "
            "four items; do not infer No from absence or calculate a partial count."
        ),
        valid_range=(0, 4),
        scorer=score_psychotic_experiences_count,
        sensitive_topics=("mental_health",),
        interpretation=(
            "Higher values indicate more distinct types of lifetime psychotic "
            "experiences, not episode frequency, recency, distress, impairment, "
            "diagnosis, or treatment."
        ),
        limitations=(
            "This is a pan-AoU-derived equal-weight count, not a validated severity scale.",
            "Conditional frequency, age, distress, professional-contact, and treatment "
            "follow-up questions are excluded.",
            "The existing three-item psych_psychotic_experiences_any phenotype is "
            "retained unchanged and is not emitted by this construction.",
        ),
    )
)

DEFINITION = PSYCHOTIC_EXPERIENCES_COUNT_DEFINITION
