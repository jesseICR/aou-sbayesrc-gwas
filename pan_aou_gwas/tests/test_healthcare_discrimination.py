import math
import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from approved_composites import INCOMPLETE, INVALID, MISSING, SCORED  # noqa: E402
import composite_healthcare_discrimination as H  # noqa: E402


LABELS = ("Never", "Rarely", "Sometimes", "Most of the time", "Always")
CONCEPTS = ("40192465", "40192481", "40192429", "40192382", "40192515")


def response_table(
    answers=None,
    *,
    ages=None,
    concepts=None,
    participant="participant",
):
    answers = list(answers if answers is not None else ["Never"] * 7)
    ages = list(ages if ages is not None else [40.0] * 7)
    concepts = list(concepts if concepts is not None else ["40192465"] * 7)
    questions = {}
    for index, qid in enumerate(H.HEALTHCARE_DISCRIMINATION_QIDS):
        if answers[index] is None:
            continue
        answer = answers[index]
        if isinstance(answer, tuple):
            answer_tuple = answer
        else:
            answer_tuple = (answer,)
        concept = concepts[index]
        if isinstance(concept, tuple):
            concept_tuple = concept
        else:
            concept_tuple = (concept,)
        questions[qid] = {
            "question": H.HEALTHCARE_DISCRIMINATION_ITEMS[index]["question"],
            "responses": {participant: (ages[index], answer_tuple)},
            "response_sidecars": {
                participant: {
                    "answer_concept_ids": concept_tuple,
                    "source_surveys": (H.SOURCE_SURVEY,) * len(answer_tuple),
                    "response_timestamps": (f"2025-01-{index + 1:02d}",)
                    * len(answer_tuple),
                }
            },
        }
    return questions


class HealthcareDiscriminationTest(unittest.TestCase):
    def test_registered_definition_metadata(self):
        definition = H.HEALTHCARE_DISCRIMINATION_DEFINITION
        self.assertEqual(definition.phenotype_id, "comp_healthcare_discrimination")
        self.assertEqual(
            definition.construction_id,
            "healthcare_discrimination_complete_case_v1",
        )
        self.assertEqual((definition.trait_type, definition.kind), ("composite", "quant"))
        self.assertEqual(definition.covar_mode, "full")
        self.assertEqual(definition.valid_range, (0, 28))
        self.assertEqual(len(definition.source_qids), 7)
        self.assertEqual(
            tuple(mapping["item_concept"] for mapping in definition.item_mappings),
            tuple(f"sdoh_dms_{index}" for index in range(1, 8)),
        )

    def test_valid_endpoints(self):
        low = H.score_healthcare_discrimination(response_table(), "participant")
        self.assertEqual(
            (low.raw_score, low.observed_component_count, low.status),
            (0, 7, SCORED),
        )

        high = H.score_healthcare_discrimination(
            response_table(["Always"] * 7, concepts=["40192515"] * 7),
            "participant",
        )
        self.assertEqual(
            (high.raw_score, high.observed_component_count, high.status),
            (28, 7, SCORED),
        )

    def test_all_intermediate_values_and_score(self):
        values = [0, 1, 2, 3, 4, 3, 2]
        score = H.score_healthcare_discrimination(
            response_table(
                [LABELS[value] for value in values],
                concepts=[CONCEPTS[value] for value in values],
            ),
            "participant",
        )
        self.assertEqual(score.raw_score, sum(values))
        self.assertEqual(score.status, SCORED)

    def test_text_fallback_and_normalization(self):
        score = H.score_healthcare_discrimination(
            response_table(
                ["  NEVER ", "Rarely", "Sometimes", "Most   of the time", "Always", "Rarely", "Never"],
                concepts=[""] * 7,
            ),
            "participant",
        )
        self.assertEqual((score.raw_score, score.status), (11, SCORED))

    def test_concept_is_authoritative_when_text_is_not_an_alias(self):
        score = H.score_healthcare_discrimination(
            response_table(
                ["release-specific label"] * 7,
                concepts=["40192429"] * 7,
            ),
            "participant",
        )
        self.assertEqual((score.raw_score, score.status), (14, SCORED))

    def test_unknown_concept_is_invalid_even_with_known_text(self):
        concepts = ["40192465"] * 7
        concepts[2] = "999999999"
        score = H.score_healthcare_discrimination(
            response_table(["Never"] * 7, concepts=concepts), "participant"
        )
        self.assertEqual(
            (score.raw_score, score.observed_component_count, score.status),
            (None, 6, INVALID),
        )

    def test_concept_text_conflict_is_invalid(self):
        answers = ["Never"] * 7
        concepts = ["40192465"] * 7
        answers[4] = "Always"
        score = H.score_healthcare_discrimination(
            response_table(answers, concepts=concepts), "participant"
        )
        self.assertEqual(
            (score.raw_score, score.observed_component_count, score.status),
            (None, 6, INVALID),
        )

    def test_missing_component_is_incomplete_not_zero(self):
        answers = ["Never"] * 7
        answers[-1] = None
        score = H.score_healthcare_discrimination(response_table(answers), "participant")
        self.assertEqual(
            (score.raw_score, score.observed_component_count, score.status),
            (None, 6, INCOMPLETE),
        )

    def test_non_substantive_response_is_missing_component(self):
        answers = ["Never"] * 7
        concepts = ["40192465"] * 7
        answers[1] = "PMI: Skip"
        concepts[1] = "903096"
        score = H.score_healthcare_discrimination(
            response_table(answers, concepts=concepts), "participant"
        )
        self.assertEqual(
            (score.raw_score, score.observed_component_count, score.status),
            (None, 6, INCOMPLETE),
        )

    def test_prefer_not_to_answer_text_fallback_is_missing(self):
        answers = ["Never"] * 7
        concepts = ["40192465"] * 7
        answers[3] = "Prefer not to answer"
        concepts[3] = ""
        score = H.score_healthcare_discrimination(
            response_table(answers, concepts=concepts), "participant"
        )
        self.assertEqual(
            (score.raw_score, score.observed_component_count, score.status),
            (None, 6, INCOMPLETE),
        )

    def test_no_responses_is_missing(self):
        score = H.score_healthcare_discrimination({}, "participant")
        self.assertEqual(
            (score.raw_score, score.observed_component_count, score.status, score.age),
            (None, 0, MISSING, None),
        )

    def test_multiple_radio_answers_are_invalid(self):
        answers = ["Never"] * 7
        concepts = ["40192465"] * 7
        answers[0] = ("Never", "Always")
        concepts[0] = ("40192465", "40192515")
        score = H.score_healthcare_discrimination(
            response_table(answers, concepts=concepts), "participant"
        )
        self.assertEqual(score.status, INVALID)
        self.assertEqual(score.observed_component_count, 6)

    def test_mean_finite_age_and_provenance(self):
        ages = [20.0, 21.0, float("nan"), 23.0, 24.0, float("inf"), 26.0]
        score = H.score_healthcare_discrimination(
            response_table(ages=ages), "participant"
        )
        self.assertTrue(math.isclose(score.age, 22.8))
        self.assertEqual(len(score.source_provenance), 7)
        first = score.source_provenance[0]
        self.assertEqual(first.qid, "40192497")
        self.assertEqual(first.source_survey, H.SOURCE_SURVEY)
        self.assertEqual(first.response_timestamp, "2025-01-01")
        self.assertEqual(first.answer_concept_ids, ("40192465",))


if __name__ == "__main__":
    unittest.main()
