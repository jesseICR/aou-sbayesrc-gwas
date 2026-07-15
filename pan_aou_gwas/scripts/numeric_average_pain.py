"""Numeric Overall Health average-pain phenotype.

The source field is represented as a slider in the survey codebook and remains
eligible for its legacy one-vs-rest binaries.  This focused builder adds the
requested numeric companion without changing those existing phenotype IDs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


PHENOTYPE_ID = "num_overallhealth_averagepain7days"
CONSTRUCTION_ID = "overallhealth_pain_slider_numeric_v1"
QUESTION_CONCEPT_ID = "1585747"
ITEM_CONCEPT = "overallhealth_averagepain7days"
QUESTION = "In the past 7 days, how would you rate your pain on average?"


@dataclass(frozen=True)
class NumericPhenotypeDefinition:
    phenotype_id: str
    construction_id: str
    source_qids: tuple[str, ...]


DEFINITION = NumericPhenotypeDefinition(
    phenotype_id=PHENOTYPE_ID,
    construction_id=CONSTRUCTION_ID,
    source_qids=(QUESTION_CONCEPT_ID,),
)


def _pain_value(answers: tuple[str, ...]) -> float | None:
    """Return a valid integer slider response as a float, otherwise missing."""
    if len(answers) != 1:
        return None
    try:
        value = float(answers[0])
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value) or not value.is_integer() or not 0.0 <= value <= 10.0:
        return None
    return value


def build_numeric_average_pain(questions: dict):
    """Yield the 0--10 average-pain slider through the quantitative GWAS path."""
    question = questions.get(QUESTION_CONCEPT_ID)
    if not question:
        return

    values = {}
    for iid, (age, answers) in question.get("responses", {}).items():
        value = _pain_value(answers)
        if value is not None:
            values[iid] = (value, age)

    yield PHENOTYPE_ID, "numeric", "quant", values, {
        "question_concept_id": QUESTION_CONCEPT_ID,
        "item_concept": ITEM_CONCEPT,
        "question": question.get("question") or QUESTION,
        "answer": "integer slider value 0 (no pain) through 10 (worst pain imaginable)",
        "ordinal_rule": "slider_numeric_0_10",
        "covar_mode": "full",
        "construction_id": CONSTRUCTION_ID,
    }
