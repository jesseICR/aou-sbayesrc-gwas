import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


QIDS = {
    "mhqukb_43": "1703923",
    "mhqukb_44": "1703928",
    "mhqukb_45": "1703906",
}

SYMPTOMS = (
    "I was more active than usual",
    "I was more talkative than usual",
    "I needed less sleep than usual",
    "I was more creative or had more ideas than usual",
    "I was more restless than usual",
    "I was more confident than usual",
    "My thoughts were racing",
    "I was easily distracted",
)


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


def build_count(rows):
    built = P.build_derived_psych_phenotypes(
        question_table(rows),
        item_labels={},
        qid_by_item=QIDS,
    )
    return next(row for row in built if row[0] == "ord_mania_symptom_count_pop")


class ManiaSymptomCountPopulationTest(unittest.TestCase):
    def test_every_valid_count_and_screener_negative_floor(self):
        rows = {
            "floor": {
                "mhqukb_43": (40.0, ("No",)),
                "mhqukb_44": (41.0, ("No",)),
            }
        }
        for count in range(9):
            selected = ("None of the above",) if count == 0 else SYMPTOMS[:count]
            rows[f"count_{count}"] = {
                "mhqukb_43": (50.0 + count, ("Yes",)),
                "mhqukb_44": (50.0 + count, ("No",)),
                "mhqukb_45": (50.0 + count, selected),
            }

        pheno_id, trait_type, kind, values, metadata = build_count(rows)

        self.assertEqual(pheno_id, "ord_mania_symptom_count_pop")
        self.assertEqual((trait_type, kind), ("derived_psych", "quant"))
        self.assertEqual(values["floor"], (-1.0, 40.0))
        for count in range(9):
            self.assertEqual(values[f"count_{count}"][0], float(count))
        self.assertEqual(
            metadata["construction_id"], "mania_symptom_count_population_zero_v1"
        )
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertIn("mental_health", metadata["sensitive_topics"])
        self.assertEqual(
            metadata["question_concept_id"], "1703923|1703928|1703906"
        )

    def test_either_screener_yes_is_an_endorser_and_symptoms_are_distinct(self):
        rows = {
            "irritable_only": {
                "mhqukb_43": (45.0, ("No",)),
                "mhqukb_44": (45.0, ("Yes",)),
                "mhqukb_45": (
                    45.0,
                    (SYMPTOMS[0], SYMPTOMS[0], "None of the above"),
                ),
            },
            "euphoric_only": {
                "mhqukb_43": (46.0, ("Yes",)),
                "mhqukb_44": (46.0, ("No",)),
                "mhqukb_45": (46.0, (SYMPTOMS[0], SYMPTOMS[1])),
            },
        }

        values = build_count(rows)[3]

        self.assertEqual(values["irritable_only"], (1.0, 45.0))
        self.assertEqual(values["euphoric_only"], (2.0, 46.0))

    def test_missing_or_unmappable_followup_does_not_become_zero(self):
        rows = {
            "one_no_one_missing": {
                "mhqukb_43": (40.0, ("No",)),
            },
            "yes_no_followup": {
                "mhqukb_43": (41.0, ("Yes",)),
                "mhqukb_44": (41.0, ("No",)),
            },
            "yes_prefer_not": {
                "mhqukb_43": (42.0, ("Yes",)),
                "mhqukb_44": (42.0, ("No",)),
                "mhqukb_45": (42.0, ("Prefer not to answer",)),
            },
            "yes_skip_mixed_with_symptom": {
                "mhqukb_43": (43.0, ("Yes",)),
                "mhqukb_44": (43.0, ("No",)),
                "mhqukb_45": (43.0, ("Skip", SYMPTOMS[0])),
            },
            "yes_unknown": {
                "mhqukb_43": (44.0, ("Yes",)),
                "mhqukb_44": (44.0, ("No",)),
                "mhqukb_45": (44.0, ("Unmapped symptom option",)),
            },
            "valid_control": {
                "mhqukb_43": (45.0, ("No",)),
                "mhqukb_44": (45.0, ("No",)),
            },
        }

        values = build_count(rows)[3]

        self.assertEqual(values, {"valid_control": (-1.0, 45.0)})

    def test_exact_filtered_run_routes_to_derived_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {"ord_mania_symptom_count_pop"},
                exact={"ord_mania_symptom_count_pop"},
            )
        )


if __name__ == "__main__":
    unittest.main()
