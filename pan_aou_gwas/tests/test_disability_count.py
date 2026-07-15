import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import approved_composites as AC  # noqa: E402
import composite_disability_count as D  # noqa: E402


IID = "participant"


def response_table(
    values=None,
    *,
    ages=None,
    surveys=None,
    use_concepts=True,
):
    values = values or {qid: 0 for qid, *_rest in D.DISABILITY_ITEMS}
    ages = ages or {}
    surveys = surveys or {}
    questions = {}
    for qid, _item, _question, yes_concept, no_concept, _pna_concept in D.DISABILITY_ITEMS:
        if qid not in values:
            continue
        value = values[qid]
        if value == 1:
            text, concept = "Yes", yes_concept
        elif value == 0:
            text, concept = "No", no_concept
        else:
            if isinstance(value, tuple):
                text, concept = value
            else:
                text, concept = str(value), ""
        questions[qid] = {
            "responses": {IID: (ages.get(qid, 40.0), (text,))},
            "response_sidecars": {
                IID: {
                    "answer_concept_ids": (concept if use_concepts else "",),
                    "source_surveys": (surveys.get(qid, "The Basics"),),
                    "response_timestamps": ("2025-01-01",),
                }
            },
        }
    return questions


class DisabilityCountTest(unittest.TestCase):
    def test_definition_is_registered_quantitative_complete_case_composite(self):
        self.assertEqual(D.DEFINITION.phenotype_id, "comp_disability_count")
        self.assertEqual(
            D.DEFINITION.construction_id,
            "disability_count_basics_life_functioning_pooled_complete_case_v1",
        )
        self.assertEqual((D.DEFINITION.trait_type, D.DEFINITION.kind), ("composite", "quant"))
        self.assertEqual(D.DEFINITION.covar_mode, "full")
        self.assertEqual(D.DEFINITION.valid_range, (0, 6))
        self.assertIn(D.DEFINITION, AC.registered_composites())

    def test_endpoints_intermediate_and_mean_finite_age(self):
        all_no = D.score_disability_count(response_table(), IID)
        self.assertEqual((all_no.raw_score, all_no.observed_component_count), (0, 6))
        self.assertEqual(all_no.status, AC.SCORED)

        all_yes_values = {qid: 1 for qid, *_rest in D.DISABILITY_ITEMS}
        all_yes = D.score_disability_count(response_table(all_yes_values), IID)
        self.assertEqual((all_yes.raw_score, all_yes.observed_component_count), (6, 6))

        intermediate = {
            qid: index % 2 for index, (qid, *_rest) in enumerate(D.DISABILITY_ITEMS)
        }
        ages = {
            qid: age
            for age, (qid, *_rest) in zip((20, 30, 40, 50, 60, 70), D.DISABILITY_ITEMS)
        }
        score = D.score_disability_count(response_table(intermediate, ages=ages), IID)
        self.assertEqual(score.raw_score, 3)
        self.assertEqual(score.age, 45.0)

        ages[D.DISABILITY_ITEMS[-1][0]] = float("nan")
        score = D.score_disability_count(response_table(intermediate, ages=ages), IID)
        self.assertEqual(score.raw_score, 3)
        self.assertEqual(score.age, 40.0)

    def test_each_missing_domain_makes_complete_case_score_missing(self):
        complete = {qid: 1 for qid, *_rest in D.DISABILITY_ITEMS}
        for qid, *_rest in D.DISABILITY_ITEMS:
            with self.subTest(qid=qid):
                values = dict(complete)
                values.pop(qid)
                score = D.score_disability_count(response_table(values), IID)
                self.assertIsNone(score.raw_score)
                self.assertEqual(score.observed_component_count, 5)
                self.assertEqual(score.status, AC.INCOMPLETE)

        no_response = D.score_disability_count({}, IID)
        self.assertEqual(no_response.status, AC.MISSING)
        self.assertEqual(no_response.observed_component_count, 0)

    def test_non_substantive_answers_are_missing_not_zero(self):
        qid, _item, _question, _yes, _no, pna = D.DISABILITY_ITEMS[2]
        questions = response_table()
        questions[qid]["responses"][IID] = (42.0, ("Prefer not to answer",))
        questions[qid]["response_sidecars"][IID]["answer_concept_ids"] = (pna,)
        score = D.score_disability_count(questions, IID)
        self.assertEqual(score.status, AC.INCOMPLETE)
        self.assertEqual(score.observed_component_count, 5)
        self.assertIsNone(score.raw_score)

        questions[qid]["responses"][IID] = (42.0, ("Don't know",))
        questions[qid]["response_sidecars"][IID]["answer_concept_ids"] = ("",)
        score = D.score_disability_count(questions, IID)
        self.assertEqual(score.status, AC.INCOMPLETE)

    def test_unknown_concept_and_concept_text_conflict_are_invalid(self):
        qid, _item, _question, yes, _no, _pna = D.DISABILITY_ITEMS[0]
        questions = response_table()
        questions[qid]["response_sidecars"][IID]["answer_concept_ids"] = ("999999999",)
        score = D.score_disability_count(questions, IID)
        self.assertEqual(score.status, AC.INVALID)
        self.assertEqual(score.observed_component_count, 5)

        questions = response_table()
        questions[qid]["responses"][IID] = (40.0, ("No",))
        questions[qid]["response_sidecars"][IID]["answer_concept_ids"] = (yes,)
        score = D.score_disability_count(questions, IID)
        self.assertEqual(score.status, AC.INVALID)
        self.assertIsNone(score.raw_score)

    def test_normalized_text_fallback_is_checked_per_qid(self):
        values = {
            qid: ("  YES  " if index < 2 else "No")
            for index, (qid, *_rest) in enumerate(D.DISABILITY_ITEMS)
        }
        score = D.score_disability_count(
            response_table(values, use_concepts=False), IID
        )
        self.assertEqual(score.status, AC.SCORED)
        self.assertEqual(score.raw_score, 2)

        # A real concept belonging to another disability qid is not accepted
        # merely because the accompanying text says Yes.
        first_qid, *_first = D.DISABILITY_ITEMS[0]
        other_yes = D.DISABILITY_ITEMS[1][3]
        questions = response_table()
        questions[first_qid]["responses"][IID] = (40.0, ("Yes",))
        questions[first_qid]["response_sidecars"][IID]["answer_concept_ids"] = (other_yes,)
        self.assertEqual(D.score_disability_count(questions, IID).status, AC.INVALID)

    def test_cross_survey_pooling_and_duplicate_domain_prevention(self):
        values = {qid: 1 for qid, *_rest in D.DISABILITY_ITEMS}
        surveys = {
            qid: ("The Basics" if index < 3 else "Life Functioning")
            for index, (qid, *_rest) in enumerate(D.DISABILITY_ITEMS)
        }
        score = D.score_disability_count(
            response_table(values, surveys=surveys), IID
        )
        self.assertEqual((score.raw_score, score.observed_component_count), (6, 6))
        self.assertEqual(
            {source.source_survey for source in score.source_provenance},
            {"The Basics", "Life Functioning"},
        )

        qid, _item, _question, yes, _no, _pna = D.DISABILITY_ITEMS[0]
        questions = response_table(values)
        questions[qid]["responses"][IID] = (40.0, ("Yes", "Yes"))
        questions[qid]["response_sidecars"][IID] = {
            "answer_concept_ids": (yes, yes),
            "source_surveys": ("The Basics", "Life Functioning"),
            "response_timestamps": ("2025-01-01", "2025-01-01"),
        }
        duplicate = D.score_disability_count(questions, IID)
        self.assertEqual((duplicate.raw_score, duplicate.observed_component_count), (6, 6))
        self.assertEqual(
            duplicate.source_provenance[0].source_survey,
            "Life Functioning|The Basics",
        )

    def test_conflicting_duplicate_domain_is_invalid(self):
        qid, _item, _question, yes, no, _pna = D.DISABILITY_ITEMS[0]
        questions = response_table()
        questions[qid]["responses"][IID] = (40.0, ("Yes", "No"))
        questions[qid]["response_sidecars"][IID] = {
            "answer_concept_ids": (yes, no),
            "source_surveys": ("The Basics", "Life Functioning"),
            "response_timestamps": ("2025-01-01", "2025-01-01"),
        }
        score = D.score_disability_count(questions, IID)
        self.assertEqual(score.status, AC.INVALID)
        self.assertEqual(score.observed_component_count, 5)


if __name__ == "__main__":
    unittest.main()
