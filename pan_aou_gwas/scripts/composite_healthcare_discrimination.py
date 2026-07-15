#!/usr/bin/env python3
"""Approved complete-case medical-setting discrimination composite."""

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


PHENOTYPE_ID = "comp_healthcare_discrimination"
CONSTRUCTION_ID = "healthcare_discrimination_complete_case_v1"
SOURCE_SURVEY = "Social Determinants of Health"

# The order is the approved sdoh_dms_1..sdoh_dms_7 order.  The question text
# is retained here so construction QC does not depend on mutable prompt text in
# a particular survey extract.
HEALTHCARE_DISCRIMINATION_ITEMS = (
    {
        "qid": "40192497",
        "item_concept": "sdoh_dms_1",
        "question": (
            "How often are you treated with less courtesy than other people when you go "
            "to a doctor's office or other health care provider?"
        ),
    },
    {
        "qid": "40192425",
        "item_concept": "sdoh_dms_2",
        "question": (
            "How often are you treated with less respect than other people when you go "
            "to a doctor's office or other health care provider?"
        ),
    },
    {
        "qid": "40192503",
        "item_concept": "sdoh_dms_3",
        "question": (
            "How often do you receive poorer service than others when you go to a "
            "doctor's office or other health care provider?"
        ),
    },
    {
        "qid": "40192505",
        "item_concept": "sdoh_dms_4",
        "question": (
            "How often does a doctor or nurse act as if he or she thinks you are not "
            "smart when you go to a doctor's office or other health care provider?"
        ),
    },
    {
        "qid": "40192423",
        "item_concept": "sdoh_dms_5",
        "question": (
            "How often does a doctor or nurse act as if he or she is afraid of you when "
            "you go to a doctor's office or other health care provider?"
        ),
    },
    {
        "qid": "40192394",
        "item_concept": "sdoh_dms_6",
        "question": (
            "How often do you feel like a doctor or nurse is not listening to what you "
            "were saying when you go to a doctor's office or other health care provider?"
        ),
    },
    {
        "qid": "40192383",
        "item_concept": "sdoh_dms_7",
        "question": (
            "How often does a doctor or nurse act as if he or she is better than you when "
            "you go to a doctor's office or other health care provider?"
        ),
    },
)

HEALTHCARE_DISCRIMINATION_QIDS = tuple(
    str(item["qid"]) for item in HEALTHCARE_DISCRIMINATION_ITEMS
)

# Answer concept IDs were verified in the local pan-AoU survey extract.  The
# metadata codebook labels these same concepts SDOH_20, SDOH_21, SDOH_22,
# SDOH_18, and SDOH_40, respectively.
ANSWER_CONCEPT_VALUES: Mapping[str, int | str] = {
    "40192465": 0,  # Never (SDOH_20)
    "40192481": 1,  # Rarely (SDOH_21)
    "40192429": 2,  # Sometimes (SDOH_22)
    "40192382": 3,  # Most of the time (SDOH_18)
    "40192515": 4,  # Always (SDOH_40)
    "903079": "missing",  # PMI: Prefer Not To Answer
    "903096": "missing",  # PMI: Skip
    "1332892": "missing",  # Prefer not to answer
}

ANSWER_TEXT_VALUES: Mapping[str, int | str] = {
    "never": 0,
    "rarely": 1,
    "sometimes": 2,
    "most of the time": 3,
    "always": 4,
    "pmi: skip": "missing",
    "skip": "missing",
    "prefer not to answer": "missing",
    "pmi: prefer not to answer": "missing",
    "pmi_prefernottoanswer": "missing",
    "don't know": "missing",
    "dont know": "missing",
    "not sure": "missing",
}


def score_healthcare_discrimination(
    questions: Mapping[str, dict], iid: str
) -> CompositeScore:
    """Score the seven-item 0..28 sum for one participant.

    All seven items must contain exactly one valid substantive response.  Known
    answer concepts take priority over text; text is only a fallback when the
    answer concept is absent.  Ages are averaged across finite component-event
    ages, independently of whether an age is available for every valid item.
    """
    values: list[int] = []
    ages: list[float | None] = []
    provenance = []
    has_invalid = False
    has_missing = False

    for qid in HEALTHCARE_DISCRIMINATION_QIDS:
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
        if value == "missing":
            has_missing = True
            continue
        if not isinstance(value, int) or value < 0 or value > 4:
            has_invalid = True
            continue

        values.append(value)
        ages.append(response.age)

    observed_count = len(values)
    if has_invalid:
        status = INVALID
        raw_score = None
    elif observed_count == len(HEALTHCARE_DISCRIMINATION_QIDS) and not has_missing:
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


HEALTHCARE_DISCRIMINATION_DEFINITION = register_composite(
    CompositeDefinition(
        phenotype_id=PHENOTYPE_ID,
        construction_id=CONSTRUCTION_ID,
        trait_type="composite",
        kind="quant",
        covar_mode="full",
        description=(
            "Pan-AoU-derived unweighted sum of seven medical-setting discrimination "
            "frequency items"
        ),
        source_surveys=(SOURCE_SURVEY,),
        source_qids=HEALTHCARE_DISCRIMINATION_QIDS,
        item_mappings=HEALTHCARE_DISCRIMINATION_ITEMS,
        answer_mapping={
            "40192465 / Never": 0,
            "40192481 / Rarely": 1,
            "40192429 / Sometimes": 2,
            "40192382 / Most of the time": 3,
            "40192515 / Always": 4,
            "903079 / PMI: Prefer Not To Answer": "missing",
            "903096 / PMI: Skip": "missing",
            "1332892 / Prefer not to answer": "missing",
        },
        missing_policy=(
            "Complete case: require one valid substantive 0..4 response for all seven "
            "items; do not prorate or infer zero from absence."
        ),
        valid_range=(0, 28),
        scorer=score_healthcare_discrimination,
        interpretation=(
            "Higher scores indicate more frequent and/or more broadly experienced "
            "discrimination in medical settings."
        ),
        limitations=(
            "This is a pan-AoU-derived equal-weight sum, not an officially validated "
            "named scale.",
        ),
    )
)

# A uniform name is convenient for registry-loading code while the explicit
# name above remains self-documenting for direct imports and tests.
DEFINITION = HEALTHCARE_DISCRIMINATION_DEFINITION
