import math
import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import approved_composites as AC  # noqa: E402
import composite_psychotic_experiences_count as P  # noqa: E402


def response_table(
    answers=None,
    *,
    concepts=None,
    ages=None,
    participant="participant",
):
    answers = list(answers if answers is not None else ["No"] * 4)
    concepts = list(concepts if concepts is not None else ["1703958"] * 4)
    ages = list(ages if ages is not None else [40.0] * 4)
    questions = {}
    for index, qid in enumerate(P.PSYCHOTIC_EXPERIENCE_QIDS):
        if answers[index] is None:
            continue
        answer = answers[index]
        answer_tuple = answer if isinstance(answer, tuple) else (answer,)
        concept = concepts[index]
        concept_tuple = concept if isinstance(concept, tuple) else (concept,)
        questions[qid] = {
            "question": P.PSYCHOTIC_EXPERIENCE_ITEMS[index]["question"],
            "responses": {participant: (ages[index], answer_tuple)},
            "response_sidecars": {
                participant: {
                    "answer_concept_ids": concept_tuple,
                    "source_surveys": (P.SOURCE_SURVEY,) * len(answer_tuple),
                    "response_timestamps": (
                        f"2025-02-{index + 1:02d}T12:00:00Z",
                    )
                    * len(answer_tuple),
                }
            },
        }
    return questions


def score(answers=None, **kwargs):
    return P.score_psychotic_experiences_count(
        response_table(answers, **kwargs), "participant"
    )


class PsychoticExperiencesCountTest(unittest.TestCase):
    def test_registered_definition_metadata_and_sensitivity(self):
        definition = P.DEFINITION
        self.assertEqual(definition.phenotype_id, "psych_psychotic_experiences_count")
        self.assertEqual(
            definition.construction_id,
            "psychotic_experiences_count_four_item_complete_case_v1",
        )
        self.assertEqual((definition.trait_type, definition.kind), ("derived_psych", "quant"))
        self.assertEqual(definition.covar_mode, "full")
        self.assertEqual(definition.valid_range, (0, 4))
        self.assertEqual(definition.sensitive_topics, ("mental_health",))
        self.assertEqual(
            definition.source_qids,
            ("1703885", "1703901", "1703915", "1703871"),
        )
        self.assertEqual(
            tuple(item["item_concept"] for item in definition.item_mappings),
            ("mhqukb_49", "cidi5_21", "cidi5_22", "cidi5_23"),
        )

    def test_valid_endpoints(self):
        minimum = score()
        maximum = score(["Yes"] * 4, concepts=["1703969"] * 4)
        self.assertEqual(
            (minimum.raw_score, minimum.observed_component_count, minimum.status),
            (0, 4, AC.SCORED),
        )
        self.assertEqual(
            (maximum.raw_score, maximum.observed_component_count, maximum.status),
            (4, 4, AC.SCORED),
        )

    def test_every_intermediate_score(self):
        for expected in range(1, 4):
            with self.subTest(expected=expected):
                answers = ["Yes"] * expected + ["No"] * (4 - expected)
                concepts = ["1703969"] * expected + ["1703958"] * (4 - expected)
                result = score(answers, concepts=concepts)
                self.assertEqual(
                    (result.raw_score, result.observed_component_count, result.status),
                    (expected, 4, AC.SCORED),
                )

    def test_visual_only_scores_one_and_corrected_any_check_is_true(self):
        result = score(
            ["Yes", "No", "No", "No"],
            concepts=["1703969", "1703958", "1703958", "1703958"],
        )
        self.assertEqual((result.raw_score, result.status), (1, AC.SCORED))

        # This is intentionally an internal consistency check only.  This
        # module must not register or emit a replacement binary phenotype.
        corrected_any = int(result.raw_score >= 1)
        self.assertEqual(corrected_any, 1)

    def test_corrected_any_check_equals_count_threshold_for_all_levels(self):
        for expected in range(5):
            answers = ["Yes"] * expected + ["No"] * (4 - expected)
            concepts = ["1703969"] * expected + ["1703958"] * (4 - expected)
            result = score(answers, concepts=concepts)
            corrected_any = int(result.raw_score >= 1)
            self.assertEqual(corrected_any, int(expected >= 1))

    def test_text_fallback_is_normalized(self):
        result = score(
            ["  YES ", "no", "Yes", " No  "],
            concepts=[""] * 4,
        )
        self.assertEqual((result.raw_score, result.status), (2, AC.SCORED))

    def test_known_concept_is_authoritative_for_unmapped_release_label(self):
        result = score(
            ["release-specific label"] * 4,
            concepts=["1703969", "1703958", "1703969", "1703958"],
        )
        self.assertEqual((result.raw_score, result.status), (2, AC.SCORED))

    def test_missing_component_is_incomplete_not_no(self):
        answers = ["No"] * 4
        answers[2] = None
        result = score(answers)
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 3, AC.INCOMPLETE),
        )

    def test_no_responses_is_missing(self):
        result = P.score_psychotic_experiences_count({}, "participant")
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status, result.age),
            (None, 0, AC.MISSING, None),
        )

    def test_non_substantive_concept_and_text_are_missing_components(self):
        cases = (
            ("PMI: Prefer Not To Answer", "903079"),
            ("Prefer not to answer", "1332892"),
            ("Don't know", ""),
            ("PMI: Skip", "903096"),
        )
        for answer, concept in cases:
            with self.subTest(answer=answer, concept=concept):
                answers = ["No"] * 4
                concepts = ["1703958"] * 4
                answers[1] = answer
                concepts[1] = concept
                result = score(answers, concepts=concepts)
                self.assertEqual(
                    (result.raw_score, result.observed_component_count, result.status),
                    (None, 3, AC.INCOMPLETE),
                )

    def test_unknown_concept_is_invalid_even_with_known_text(self):
        concepts = ["1703958"] * 4
        concepts[0] = "999999999"
        result = score(["No"] * 4, concepts=concepts)
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 3, AC.INVALID),
        )

    def test_unrecognized_text_without_concept_is_invalid(self):
        answers = ["No"] * 4
        concepts = ["1703958"] * 4
        answers[3] = "Unrecognized answer"
        concepts[3] = ""
        result = score(answers, concepts=concepts)
        self.assertEqual((result.raw_score, result.status), (None, AC.INVALID))

    def test_concept_text_conflict_is_invalid(self):
        answers = ["No"] * 4
        concepts = ["1703958"] * 4
        answers[0] = "Yes"
        result = score(answers, concepts=concepts)
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 3, AC.INVALID),
        )

    def test_multiple_radio_answers_are_invalid(self):
        answers = ["No"] * 4
        concepts = ["1703958"] * 4
        answers[0] = ("No", "Yes")
        concepts[0] = ("1703958", "1703969")
        result = score(answers, concepts=concepts)
        self.assertEqual((result.raw_score, result.status), (None, AC.INVALID))

    def test_mean_finite_age_and_provenance(self):
        result = score(ages=[20.0, float("nan"), 24.0, float("inf")])
        self.assertTrue(math.isclose(result.age, 22.0))
        self.assertEqual(len(result.source_provenance), 4)
        first = result.source_provenance[0]
        self.assertEqual(first.qid, "1703885")
        self.assertEqual(first.source_survey, P.SOURCE_SURVEY)
        self.assertEqual(first.answer_concept_ids, ("1703958",))
        self.assertEqual(first.response_timestamp, "2025-02-01T12:00:00Z")


if __name__ == "__main__":
    unittest.main()
