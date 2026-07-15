import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


QIDS = {
    "mhqukb_5": "dep5",
    "mhqukb_6": "dep6",
    "mhqukb_21": "dep21",
    "mhqukb_22": "dep22",
    "mhqukb_24": "dep24",
}


def question_table(rows):
    questions = {}
    for item, qid in QIDS.items():
        responses = {
            iid: (age, tuple(answers))
            for iid, participant in rows.items()
            if (entry := participant.get(item)) is not None
            for age, answers in (entry,)
        }
        questions[qid] = {"question": item, "responses": responses}
    return questions


def build_population_traits(rows):
    return {
        built[0]: built
        for built in P.build_population_gated_phenotypes(
            question_table(rows),
            item_labels={},
            qid_by_item=QIDS,
        )
    }


class DepressionImpairmentPopulationTest(unittest.TestCase):
    def test_floor_and_every_valid_impairment_level(self):
        rows = {
            "floor": {
                "mhqukb_5": (40.0, ("No",)),
                "mhqukb_6": (41.0, ("No",)),
            },
            "not_at_all": {
                "mhqukb_5": (42.0, ("Yes",)),
                "mhqukb_6": (42.0, ("No",)),
                "mhqukb_22": (43.0, ("Not at all",)),
            },
            "a_little": {
                "mhqukb_5": (44.0, ("No",)),
                "mhqukb_6": (44.0, ("Yes",)),
                "mhqukb_22": (45.0, ("A little",)),
            },
            "somewhat": {
                "mhqukb_5": (46.0, ("Yes",)),
                "mhqukb_22": (47.0, ("Somewhat",)),
            },
            "a_lot": {
                "mhqukb_6": (48.0, ("Yes",)),
                "mhqukb_22": (49.0, ("A lot",)),
            },
        }

        pheno_id, trait_type, kind, values, metadata = build_population_traits(rows)[
            "ord_mhqukb_22_pop"
        ]

        self.assertEqual(pheno_id, "ord_mhqukb_22_pop")
        self.assertEqual((trait_type, kind), ("ordinal", "quant"))
        self.assertEqual(
            values,
            {
                "floor": (-1.0, 40.0),
                "not_at_all": (0.0, 43.0),
                "a_little": (1.0, 45.0),
                "somewhat": (2.0, 47.0),
                "a_lot": (3.0, 49.0),
            },
        )
        self.assertEqual(
            metadata["construction_id"],
            "mhq_depression_followup_population_zero_v1",
        )
        self.assertEqual(metadata["ordinal_rule"], "impact_0_3_population_zero")
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertEqual(metadata["question_concept_id"], "dep5|dep6|dep22")

    def test_missing_ambiguous_and_non_substantive_answers_are_omitted(self):
        rows = {
            "one_no_one_missing": {
                "mhqukb_5": (40.0, ("No",)),
            },
            "positive_without_followup": {
                "mhqukb_5": (41.0, ("Yes",)),
                "mhqukb_6": (41.0, ("No",)),
            },
            "positive_prefer_not": {
                "mhqukb_5": (42.0, ("Yes",)),
                "mhqukb_6": (42.0, ("No",)),
                "mhqukb_22": (42.0, ("Prefer not to answer",)),
            },
            "positive_unknown": {
                "mhqukb_5": (43.0, ("Yes",)),
                "mhqukb_6": (43.0, ("No",)),
                "mhqukb_22": (43.0, ("Unknown impact level",)),
            },
            "valid_floor": {
                "mhqukb_5": (44.0, ("No",)),
                "mhqukb_6": (44.0, ("No",)),
                "mhqukb_22": (45.0, ("A lot",)),
            },
        }

        values = build_population_traits(rows)["ord_mhqukb_22_pop"][3]

        self.assertEqual(values, {"valid_floor": (-1.0, 44.0)})

    def test_existing_depression_population_traits_keep_zero_floor(self):
        rows = {
            "floor": {
                "mhqukb_5": (40.0, ("No",)),
                "mhqukb_6": (40.0, ("No",)),
            }
        }

        built = build_population_traits(rows)

        self.assertEqual(built["ord_mhqukb_21_pop"][3]["floor"], (0.0, 40.0))
        self.assertEqual(built["ord_mhqukb_24_pop"][3]["floor"], (0.0, 40.0))
        self.assertEqual(built["ord_mhqukb_22_pop"][3]["floor"], (-1.0, 40.0))

    def test_exact_filtered_run_routes_to_population_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {"ord_mhqukb_22_pop"},
                prefixes=("num_", "ord_"),
            )
        )


if __name__ == "__main__":
    unittest.main()
