import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import approved_composites as AC  # noqa: E402
import composite_housing_problem_count as H  # noqa: E402


def response_table(observations=None, age=47.25, include_response=True):
    responses = {}
    sidecars = {}
    if include_response:
        observations = tuple(observations or ())
        responses["participant"] = (age, tuple(text for text, _concept_id in observations))
        sidecars["participant"] = {
            "answer_concept_ids": tuple(concept_id for _text, concept_id in observations),
            "source_surveys": tuple(H.SOURCE_SURVEY for _ in observations),
            "response_timestamps": tuple("2024-01-02T03:04:05Z" for _ in observations),
        }
    return {
        H.QUESTION_CONCEPT_ID: {
            "question": "Think about the place you live.",
            "responses": responses,
            "response_sidecars": sidecars,
        }
    }


def score(observations=None, **kwargs):
    return H.score_housing_problem_count(response_table(observations, **kwargs), "participant")


class HousingProblemCountTest(unittest.TestCase):
    def test_none_only_is_observed_zero(self):
        result = score((("None of the above", "40192392"),))
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (0, 7, AC.SCORED),
        )

    def test_all_problems_and_intermediate_score(self):
        all_problems = tuple((label, concept_id) for _code, concept_id, label in H.PROBLEM_OPTIONS)
        maximum = score(all_problems)
        self.assertEqual((maximum.raw_score, maximum.status), (7, AC.SCORED))

        intermediate = score(
            (
                ("Bug infestation", "40192460"),
                ("Inadequate heat", "40192434"),
                ("Water leaks", "40192444"),
            )
        )
        self.assertEqual(
            (intermediate.raw_score, intermediate.observed_component_count), (3, 7)
        )

    def test_text_fallback_and_duplicate_option(self):
        result = score((("  mold ", ""), ("SDOH_45", "")))
        self.assertEqual((result.raw_score, result.status), (1, AC.SCORED))

    def test_none_plus_problem_is_contradiction(self):
        result = score(
            (("None of the above", "40192392"), ("Mold", "40192479"))
        )
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 7, AC.CONTRADICTION),
        )

    def test_absent_response_and_no_selection_are_missing(self):
        absent = score(include_response=False)
        empty = score(())
        self.assertEqual((absent.raw_score, absent.status), (None, AC.MISSING))
        self.assertEqual((empty.raw_score, empty.status), (None, AC.MISSING))
        self.assertEqual(empty.observed_component_count, 0)

    def test_non_substantive_answers_are_missing(self):
        for observation in (
            ("PMI: Prefer Not To Answer", "903079"),
            ("Prefer not to answer", ""),
            ("PMI: Skip", "903096"),
        ):
            with self.subTest(observation=observation):
                result = score((observation,))
                self.assertEqual((result.raw_score, result.status), (None, AC.MISSING))

    def test_unknown_concept_is_invalid_even_with_known_text(self):
        result = score((("Mold", "99999999"),))
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 0, AC.INVALID),
        )

    def test_unrecognized_text_without_concept_is_missing(self):
        result = score((("Unmapped checkbox label", ""),))
        self.assertEqual((result.raw_score, result.status), (None, AC.MISSING))

    def test_concept_text_conflict_is_invalid(self):
        result = score((("Mold", "40192460"),))  # concept denotes Bug infestation
        self.assertEqual((result.raw_score, result.status), (None, AC.INVALID))

    def test_response_event_age_and_provenance_are_retained(self):
        result = score((("Water leaks", "40192444"),), age=61.5)
        self.assertEqual(result.age, 61.5)
        self.assertEqual(len(result.source_provenance), 1)
        provenance = result.source_provenance[0]
        self.assertEqual(provenance.qid, H.QUESTION_CONCEPT_ID)
        self.assertEqual(provenance.source_survey, H.SOURCE_SURVEY)
        self.assertEqual(provenance.answer_concept_ids, ("40192444",))
        self.assertEqual(provenance.response_timestamp, "2024-01-02T03:04:05Z")

    def test_registered_definition_is_quantitative_complete_case_composite(self):
        definition = H.DEFINITION
        self.assertEqual(definition.phenotype_id, "comp_housing_problem_count")
        self.assertEqual(
            definition.construction_id, "housing_problem_count_complete_case_v1"
        )
        self.assertEqual((definition.trait_type, definition.kind), ("composite", "quant"))
        self.assertEqual(definition.valid_range, (0, 7))
        self.assertEqual(definition.source_qids, ("40192402",))


if __name__ == "__main__":
    unittest.main()
