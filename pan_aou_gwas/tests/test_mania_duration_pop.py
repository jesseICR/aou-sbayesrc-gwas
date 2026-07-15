import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


QIDS = {
    "mhqukb_43": "1703923",
    "mhqukb_44": "1703928",
    "mhqukb_46": "1703879",
}

DURATION_VALUES = {
    "Less than 24 hours": 0.5,
    "At least a day, but less than four days": 2.5,
    "At least four days in a row but less than a week": 5.5,
    "A week or more": 10.0,
}


def question_table(rows):
    questions = {}
    for item, qid in QIDS.items():
        responses = {
            iid: (float(age), tuple(answers))
            for iid, participant in rows.items()
            if (entry := participant.get(item)) is not None
            for age, answers in (entry,)
        }
        questions[qid] = {"question": item, "responses": responses}
    return questions


def build_duration(rows):
    built = P.build_derived_psych_phenotypes(
        question_table(rows),
        item_labels={},
        qid_by_item=QIDS,
    )
    return next(row for row in built if row[0] == "ord_mhqukb_46_pop")


class ManiaDurationPopulationTest(unittest.TestCase):
    def test_floor_all_duration_bands_and_metadata(self):
        rows = {
            "floor": {
                "mhqukb_43": (40, ("No",)),
                "mhqukb_44": (41, ("No",)),
            }
        }
        for i, answer in enumerate(DURATION_VALUES):
            rows[f"duration_{i}"] = {
                "mhqukb_43": (50 + i, ("Yes",)),
                "mhqukb_44": (50 + i, ("No",)),
                "mhqukb_46": (60 + i, (answer,)),
            }

        pheno_id, trait_type, kind, values, metadata = build_duration(rows)

        self.assertEqual(pheno_id, "ord_mhqukb_46_pop")
        self.assertEqual((trait_type, kind), ("ordinal", "quant"))
        self.assertEqual(values["floor"], (-1.0, 40.0))
        for i, expected in enumerate(DURATION_VALUES.values()):
            self.assertEqual(values[f"duration_{i}"], (expected, 60.0 + i))
        self.assertEqual(
            metadata["construction_id"], "mania_followup_population_zero_v1"
        )
        self.assertEqual(
            metadata["ordinal_rule"],
            "episode_duration_days_midpoint_population_zero",
        )
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertIn("mental_health", metadata["sensitive_topics"])
        self.assertEqual(
            metadata["question_concept_id"], "1703923|1703928|1703879"
        )

    def test_either_screener_yes_opens_the_followup(self):
        rows = {
            "euphoric": {
                "mhqukb_43": (40, ("Yes",)),
                "mhqukb_44": (40, ("No",)),
                "mhqukb_46": (40, ("Less than 24 hours",)),
            },
            "irritable": {
                "mhqukb_43": (41, ("No",)),
                "mhqukb_44": (41, ("Yes",)),
                "mhqukb_46": (41, ("A week or more",)),
            },
        }

        values = build_duration(rows)[3]

        self.assertEqual(values["euphoric"][0], 0.5)
        self.assertEqual(values["irritable"][0], 10.0)

    def test_incomplete_gate_or_invalid_followup_remains_missing(self):
        rows = {
            "one_no_one_missing": {
                "mhqukb_43": (40, ("No",)),
            },
            "endorser_no_followup": {
                "mhqukb_43": (41, ("Yes",)),
                "mhqukb_44": (41, ("No",)),
            },
            "endorser_prefer_not": {
                "mhqukb_43": (42, ("Yes",)),
                "mhqukb_44": (42, ("No",)),
                "mhqukb_46": (42, ("Prefer not to answer",)),
            },
            "endorser_dont_know": {
                "mhqukb_43": (43, ("No",)),
                "mhqukb_44": (43, ("Yes",)),
                "mhqukb_46": (43, ("Don't know",)),
            },
            "endorser_unknown": {
                "mhqukb_43": (44, ("No",)),
                "mhqukb_44": (44, ("Yes",)),
                "mhqukb_46": (44, ("About two weeks",)),
            },
            "valid_floor": {
                "mhqukb_43": (45, ("No",)),
                "mhqukb_44": (45, ("No",)),
            },
        }

        values = build_duration(rows)[3]

        self.assertEqual(values, {"valid_floor": (-1.0, 45.0)})

    def test_exact_id_routes_to_the_derived_psych_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {"ord_mhqukb_46_pop"},
                exact={"ord_mhqukb_46_pop"},
            )
        )


if __name__ == "__main__":
    unittest.main()
