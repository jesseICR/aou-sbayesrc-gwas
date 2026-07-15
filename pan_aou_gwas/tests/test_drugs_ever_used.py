import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import approved_composites as AC  # noqa: E402
import composite_drugs_ever_used as D  # noqa: E402


IID = "participant"


def response_table(observations=None, *, age=48.5, include_response=True):
    responses = {}
    sidecars = {}
    if include_response:
        observations = tuple(observations or ())
        responses[IID] = (age, tuple(text for text, _concept_id in observations))
        sidecars[IID] = {
            "answer_concept_ids": tuple(
                concept_id for _text, concept_id in observations
            ),
            "source_surveys": tuple(D.SOURCE_SURVEY for _ in observations),
            "response_timestamps": tuple(
                "2025-03-04T05:06:07Z" for _ in observations
            ),
        }
    return {
        D.QUESTION_CONCEPT_ID: {
            "question": "Which of these drugs have you ever used?",
            "responses": responses,
            "response_sidecars": sidecars,
        }
    }


def score(observations=None, **kwargs):
    return D.score_drugs_ever_used(response_table(observations, **kwargs), IID)


class DrugsEverUsedTest(unittest.TestCase):
    def test_registered_quantitative_definition_and_limit_metadata(self):
        definition = D.DEFINITION
        self.assertEqual(definition.phenotype_id, "num_drugs_ever_used")
        self.assertEqual(
            definition.construction_id,
            "drugs_ever_used_nine_class_checkbox_v1",
        )
        self.assertEqual((definition.trait_type, definition.kind), ("composite", "quant"))
        self.assertEqual(definition.covar_mode, "full")
        self.assertEqual(definition.valid_range, (0, 9))
        self.assertEqual(definition.source_qids, ("1585636",))
        self.assertEqual(len(D.DRUG_CLASS_OPTIONS), 9)
        self.assertIn(definition, AC.registered_composites())
        limitations = " ".join(definition.limitations).lower()
        for required in ("breadth", "severity", "frequency", "unequal", "self-reported", "other"):
            self.assertIn(required, limitations)

    def test_none_only_is_zero_with_response_event_age_and_provenance(self):
        result = score(((D.NONE_OPTION[2], D.NONE_OPTION[1]),), age=51.25)
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (0, 9, AC.SCORED),
        )
        self.assertEqual(result.age, 51.25)
        self.assertEqual(len(result.source_provenance), 1)
        provenance = result.source_provenance[0]
        self.assertEqual(provenance.qid, "1585636")
        self.assertEqual(provenance.source_survey, "Lifestyle")
        self.assertEqual(provenance.response_timestamp, "2025-03-04T05:06:07Z")
        self.assertEqual(provenance.answer_concept_ids, ("1585648",))

    def test_all_valid_endpoints_and_intermediate_scores(self):
        all_classes = tuple(
            (label, concept_id)
            for _answer_code, concept_id, label in D.DRUG_CLASS_OPTIONS
        )
        maximum = score(all_classes)
        self.assertEqual(
            (maximum.raw_score, maximum.observed_component_count, maximum.status),
            (9, 9, AC.SCORED),
        )

        for expected in range(1, 9):
            with self.subTest(expected=expected):
                result = score(all_classes[:expected])
                self.assertEqual(
                    (result.raw_score, result.observed_component_count, result.status),
                    (expected, 9, AC.SCORED),
                )

    def test_duplicate_class_is_counted_once_and_text_fallback_is_checked(self):
        result = score(
            (
                ("Marijuana", ""),
                ("WhichDrugsUsed_MarijuanaUse", ""),
                ("Cocaine", ""),
            )
        )
        self.assertEqual((result.raw_score, result.status), (2, AC.SCORED))

    def test_none_plus_scored_class_is_contradiction(self):
        result = score(
            (
                (D.NONE_OPTION[2], D.NONE_OPTION[1]),
                (D.DRUG_CLASS_OPTIONS[0][2], D.DRUG_CLASS_OPTIONS[0][1]),
            )
        )
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 9, AC.CONTRADICTION),
        )

    def test_absent_response_and_no_selection_are_missing(self):
        absent = score(include_response=False)
        empty = score(())
        self.assertEqual(
            (absent.raw_score, absent.observed_component_count, absent.status),
            (None, 0, AC.MISSING),
        )
        self.assertEqual(
            (empty.raw_score, empty.observed_component_count, empty.status),
            (None, 0, AC.MISSING),
        )

    def test_other_is_ignored_but_other_only_is_missing(self):
        other_only = score(((D.OTHER_OPTION[2], D.OTHER_OPTION[1]),))
        self.assertEqual(
            (other_only.raw_score, other_only.observed_component_count, other_only.status),
            (None, 0, AC.MISSING),
        )

        scored_plus_other = score(
            (
                (D.OTHER_OPTION[2], D.OTHER_OPTION[1]),
                (D.DRUG_CLASS_OPTIONS[0][2], D.DRUG_CLASS_OPTIONS[0][1]),
            )
        )
        self.assertEqual(
            (scored_plus_other.raw_score, scored_plus_other.status),
            (1, AC.SCORED),
        )

        none_plus_other = score(
            (
                (D.NONE_OPTION[2], D.NONE_OPTION[1]),
                (D.OTHER_OPTION[2], D.OTHER_OPTION[1]),
            )
        )
        self.assertEqual(
            (none_plus_other.raw_score, none_plus_other.status),
            (0, AC.SCORED),
        )

    def test_prefer_not_to_answer_and_skips_are_missing(self):
        for observation in (
            ("PMI: Prefer Not To Answer", "903079"),
            ("Prefer not to answer", ""),
            ("PMI: Skip", "903096"),
            ("Don't know", ""),
        ):
            with self.subTest(observation=observation):
                result = score((observation,))
                self.assertEqual(
                    (result.raw_score, result.observed_component_count, result.status),
                    (None, 0, AC.MISSING),
                )

    def test_unknown_concept_is_invalid_even_with_known_text(self):
        result = score((("Marijuana", "99999999"),))
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 0, AC.INVALID),
        )

    def test_concept_text_conflict_is_invalid(self):
        marijuana_concept = D.DRUG_CLASS_OPTIONS[0][1]
        cocaine_label = D.DRUG_CLASS_OPTIONS[1][2]
        result = score(((cocaine_label, marijuana_concept),))
        self.assertEqual(
            (result.raw_score, result.observed_component_count, result.status),
            (None, 0, AC.INVALID),
        )

    def test_unrecognized_text_without_concept_is_missing(self):
        result = score((("Unmapped checkbox option", ""),))
        self.assertEqual((result.raw_score, result.status), (None, AC.MISSING))


if __name__ == "__main__":
    unittest.main()
