#!/usr/bin/env python3
"""Complete-case count of the six pooled ACS disability domains."""

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


PHENOTYPE_ID = "comp_disability_count"
CONSTRUCTION_ID = "disability_count_basics_life_functioning_pooled_complete_case_v1"
SOURCE_SURVEYS = ("The Basics", "Life Functioning: Basics Survey Companion")

DISABILITY_ITEMS = (
    (
        "903573",
        "disability_deaf",
        "Are you deaf or do you have serious difficulty hearing?",
        "903587",
        "903503",
        "903596",
    ),
    (
        "903574",
        "disability_blind",
        "Are you blind or do you have serious difficulty seeing, even when wearing glasses?",
        "903504",
        "903597",
        "903598",
    ),
    (
        "903575",
        "disability_difficultyconcentrating",
        (
            "Because of a physical, mental, or emotional condition, do you have serious "
            "difficulty concentrating, remembering or making decisions?"
        ),
        "903599",
        "903600",
        "903601",
    ),
    (
        "903576",
        "disability_walkingclimbing",
        "Do you have serious difficulty walking or climbing stairs?",
        "903602",
        "903603",
        "903604",
    ),
    (
        "903577",
        "disability_dressingbathing",
        "Do you have difficulty dressing or bathing?",
        "903605",
        "903606",
        "903607",
    ),
    (
        "903578",
        "disability_errandsalone",
        (
            "Because of a physical, mental, or emotional condition, do you have difficulty "
            "doing errands alone such as visiting doctor's office or shopping?"
        ),
        "903608",
        "903609",
        "903610",
    ),
)

# AoU answer concepts are item-specific for these otherwise identical Yes/No
# scales.  The two generic PMI concepts are included because older releases can
# encode non-substantive responses with them instead of the item-specific PNA.
ANSWER_CONCEPT_VALUES = {
    qid: {
        yes_concept: 1,
        no_concept: 0,
        pna_concept: None,
        "903079": None,  # PMI_PreferNotToAnswer
        "903096": None,  # PMI_Skip
    }
    for qid, _item, _question, yes_concept, no_concept, pna_concept in DISABILITY_ITEMS
}

ANSWER_TEXT_VALUES = {
    "yes": 1,
    "no": 0,
    "prefer not to answer": None,
    "pmi: prefer not to answer": None,
    "don't know": None,
    "dont know": None,
    "not sure": None,
    "skip": None,
    "pmi: skip": None,
    "skipped": None,
}


def score_disability_count(questions: Mapping[str, dict], iid: str) -> CompositeScore:
    """Score one participant after latest-valid responses have been selected.

    Basics and Life Functioning use the same six qids, so pooling occurs before
    this function in the shared latest-response ingest.  If duplicate rows from
    both surveys remain at the selected event, identical substantive answers
    define one domain; conflicting answers make the composite invalid.
    """
    component_values: list[int] = []
    component_ages: list[float | None] = []
    provenance = []
    invalid = False
    has_response = False

    for qid, _item, _question, _yes, _no, _pna in DISABILITY_ITEMS:
        response = selected_response(questions, qid, iid)
        if response is None:
            continue

        has_response = True
        item_provenance = response_provenance(response)
        if item_provenance is not None:
            provenance.append(item_provenance)

        substantive: list[int] = []
        non_substantive = False
        item_invalid = False
        for observation in response.answers:
            value, status = resolve_answer(
                observation,
                ANSWER_CONCEPT_VALUES[qid],
                ANSWER_TEXT_VALUES,
            )
            if status == INVALID:
                item_invalid = True
            elif value is None:
                non_substantive = True
            else:
                substantive.append(int(value))

        distinct_values = set(substantive)
        if item_invalid or len(distinct_values) > 1 or (non_substantive and substantive):
            invalid = True
            continue
        if non_substantive or not substantive:
            continue

        # Multiple identical observations can arise from the pooled sources at
        # one timestamp.  They are one domain, never extra count components.
        component_values.append(next(iter(distinct_values)))
        component_ages.append(response.age)

    observed = len(component_values)
    age = finite_mean(component_ages)
    if invalid:
        return CompositeScore(None, observed, INVALID, age, tuple(provenance))
    if observed == len(DISABILITY_ITEMS):
        return CompositeScore(sum(component_values), observed, SCORED, age, tuple(provenance))
    status = INCOMPLETE if has_response else MISSING
    return CompositeScore(None, observed, status, age, tuple(provenance))


DISABILITY_COUNT_DEFINITION = register_composite(
    CompositeDefinition(
        phenotype_id=PHENOTYPE_ID,
        construction_id=CONSTRUCTION_ID,
        trait_type="composite",
        description=(
            "Complete-case count of six ACS disability domains endorsed across pooled "
            "The Basics and Life Functioning responses."
        ),
        source_surveys=SOURCE_SURVEYS,
        source_qids=tuple(item[0] for item in DISABILITY_ITEMS),
        item_mappings=tuple(
            {
                "item_concept": item,
                "question_concept_id": qid,
                "question": question,
                "answer_concept_ids": {
                    "Yes": yes_concept,
                    "No": no_concept,
                    "Prefer not to answer": pna_concept,
                },
            }
            for qid, item, question, yes_concept, no_concept, pna_concept in DISABILITY_ITEMS
        ),
        answer_mapping={"No": 0, "Yes": 1, "Prefer not to answer": "missing"},
        missing_policy=(
            "Complete case: require one valid Yes/No answer for every one of the six qids; "
            "do not infer No or prorate missing, non-substantive, or invalid domains."
        ),
        valid_range=(0, 6),
        scorer=score_disability_count,
        sensitive_topics=("disability",),
        interpretation=(
            "Higher values indicate more disability domains endorsed; this pan-AoU-derived "
            "count is not symptom severity or an official ACS severity scale."
        ),
        limitations=(
            "Equal-weight domain breadth, not symptom severity.",
            "Requires complete valid responses to all six domains.",
        ),
        kind="quant",
        covar_mode="full",
    )
)

# Uniform alias used by the shared registration loader.
DEFINITION = DISABILITY_COUNT_DEFINITION
