import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


QIDS = {
    "mhqukb_49": "1703885",
    "cidi5_21": "1703901",
    "cidi5_22": "1703915",
    "cidi5_23": "1703871",
    "mhqukb_54": "1703890",
}
SCREENERS = ("mhqukb_49", "cidi5_21", "cidi5_22", "cidi5_23")
DISTRESS_ANSWERS = (
    "Not distressing at all. It was a positive experience",
    "Not distressing, a neutral experience",
    "A bit distressing",
    "Quite distressing",
    "Very distressing",
)


def question_table(rows):
    questions = {}
    for item, qid in QIDS.items():
        responses = {
            iid: (age, answers if isinstance(answers, tuple) else (answers,))
            for iid, participant in rows.items()
            if (entry := participant.get(item)) is not None
            for age, answers in (entry,)
        }
        questions[qid] = {"question": item, "responses": responses}
    return questions


def build_distress(rows):
    built = P.build_derived_psych_phenotypes(
        question_table(rows),
        item_labels={},
        qid_by_item=QIDS,
    )
    return next(row for row in built if row[0] == "ord_mhqukb_54_pop")


def screen_answers(yes_item=None, age=40.0):
    return {
        item: (age, "Yes" if item == yes_item else "No")
        for item in SCREENERS
    }


class PsychosisDistressPopulationTest(unittest.TestCase):
    def test_screener_negative_floor_and_all_followup_levels(self):
        rows = {"floor": screen_answers(age=30.0)}
        for expected, answer in enumerate(DISTRESS_ANSWERS):
            rows[f"level_{expected}"] = {
                **screen_answers("mhqukb_49", age=40.0 + expected),
                "mhqukb_54": (50.0 + expected, answer),
            }

        pheno_id, trait_type, kind, values, metadata = build_distress(rows)

        self.assertEqual(pheno_id, "ord_mhqukb_54_pop")
        self.assertEqual((trait_type, kind), ("ordinal", "quant"))
        self.assertEqual(values["floor"], (-1.0, 30.0))
        for expected in range(5):
            self.assertEqual(
                values[f"level_{expected}"],
                (float(expected), 50.0 + expected),
            )
        self.assertEqual(
            metadata["construction_id"], "psychosis_distress_population_zero_v1"
        )
        self.assertEqual(metadata["ordinal_rule"], "distress_0_4_population_floor")
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertIn("mental_health", metadata["sensitive_topics"])
        self.assertIn("psychosis", metadata["sensitive_topics"])
        self.assertEqual(
            metadata["question_concept_id"],
            "1703885|1703901|1703915|1703871|1703890",
        )
        self.assertEqual(
            metadata["item_concept"],
            "mhqukb_49|cidi5_21|cidi5_22|cidi5_23|mhqukb_54",
        )

    def test_each_of_four_screeners_can_open_the_followup_gate(self):
        rows = {
            item: {
                **screen_answers(item, age=45.0),
                "mhqukb_54": (46.0, "A bit distressing"),
            }
            for item in SCREENERS
        }

        values = build_distress(rows)[3]

        for item in SCREENERS:
            self.assertEqual(values[item], (2.0, 46.0))

    def test_incomplete_negative_screen_does_not_earn_floor(self):
        rows = {
            "one_missing": {
                item: (40.0, "No")
                for item in SCREENERS[:-1]
            },
            "one_prefer_not": {
                **screen_answers(age=41.0),
                "cidi5_23": (41.0, "Prefer not to answer"),
            },
            "valid_floor": screen_answers(age=42.0),
        }

        values = build_distress(rows)[3]

        self.assertEqual(values, {"valid_floor": (-1.0, 42.0)})

    def test_endorser_requires_one_valid_distress_answer(self):
        rows = {
            "no_followup": screen_answers("cidi5_21"),
            "prefer_not": {
                **screen_answers("cidi5_21"),
                "mhqukb_54": (40.0, "Prefer not to answer"),
            },
            "unknown": {
                **screen_answers("cidi5_21"),
                "mhqukb_54": (40.0, "Unmapped distress response"),
            },
            "multiple": {
                **screen_answers("cidi5_21"),
                "mhqukb_54": (
                    40.0,
                    ("A bit distressing", "Very distressing"),
                ),
            },
            "valid": {
                **screen_answers("cidi5_21"),
                "mhqukb_54": (43.0, "Very distressing"),
            },
        }

        values = build_distress(rows)[3]

        self.assertEqual(values, {"valid": (4.0, 43.0)})

    def test_exact_filtered_run_routes_to_derived_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {"ord_mhqukb_54_pop"},
                exact={"ord_mhqukb_54_pop"},
            )
        )


if __name__ == "__main__":
    unittest.main()
