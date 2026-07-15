import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


EDS_QIDS = P.EDS_BREADTH_SOURCE_QIDS
FOLLOWUP_QID = P.EDS_ATTRIBUTION_QID
GROUNDS = (
    "Your ancestry or national origins",
    "Your gender",
    "Your race",
    "Your age",
    "Your religion",
    "Your height",
    "Your weight",
    "Some other aspect of your physical appearance",
    "Your sexual orientation",
    "Your education or income level",
)


def question_table(rows):
    questions = {}
    for index, qid in enumerate(EDS_QIDS):
        responses = {
            iid: (float(age), tuple(answers))
            for iid, participant in rows.items()
            if (entry := participant.get(f"eds_{index + 1}")) is not None
            for age, answers in (entry,)
        }
        questions[qid] = {"question": f"sdoh_eds_{index + 1}", "responses": responses}
    questions[FOLLOWUP_QID] = {
        "question": "sdoh_eds_follow_up_1",
        "responses": {
            iid: (float(age), tuple(answers))
            for iid, participant in rows.items()
            if (entry := participant.get("followup")) is not None
            for age, answers in (entry,)
        },
    }
    return questions


def all_never(age=40.0):
    return {f"eds_{index}": (age + index, ("Never",)) for index in range(1, 10)}


def reporter(followup, *, age=50.0):
    return {
        "eds_1": (age, ("A few times a year",)),
        "followup": (age + 5.0, followup),
    }


def build(rows):
    built = P.build_population_gated_phenotypes(
        question_table(rows),
        item_labels={},
        qid_by_item={},
    )
    return next(row for row in built if row[0] == P.EDS_BREADTH_PHENO_ID)


class SdohEdsDiscriminationBreadthPopulationTest(unittest.TestCase):
    def test_floor_and_every_valid_endpoint(self):
        rows = {
            "floor": all_never(),
            "one": reporter((GROUNDS[0],)),
            "ten": reporter(GROUNDS),
        }

        pheno_id, trait_type, kind, values, metadata = build(rows)

        self.assertEqual(pheno_id, "ord_sdoh_eds_discrimination_breadth_pop")
        self.assertEqual((trait_type, kind), ("derived", "quant"))
        self.assertEqual(values["floor"], (-1.0, 41.0))
        self.assertEqual(values["one"], (1.0, 55.0))
        self.assertEqual(values["ten"], (10.0, 55.0))
        self.assertEqual(
            metadata["construction_id"],
            "eds_attribution_breadth_population_zero_v1",
        )
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertEqual(
            metadata["question_concept_id"],
            "|".join((*EDS_QIDS, FOLLOWUP_QID)),
        )

    def test_distinct_ground_count_ignores_known_other(self):
        rows = {
            "duplicate": reporter((GROUNDS[0], GROUNDS[0], GROUNDS[1])),
            "with_other": reporter((GROUNDS[0], "Other (specify)")),
            "code_aliases": reporter(("SDOH_29", "SDOH_38", "SDOH_38")),
        }

        values = build(rows)[3]

        self.assertEqual(values["duplicate"][0], 2.0)
        self.assertEqual(values["with_other"][0], 1.0)
        self.assertEqual(values["code_aliases"][0], 2.0)

    def test_incomplete_gate_and_uncodeable_followups_are_missing(self):
        rows = {
            "incomplete_all_never": {
                "eds_1": (40, ("Never",)),
            },
            "unknown_gate_answer": {
                **all_never(),
                "eds_9": (49, ("Occasionally-ish",)),
            },
            "reporter_no_followup": {
                "eds_1": (50, ("A few times a month",)),
            },
            "reporter_other_only": reporter(("Other (specify)",)),
            "reporter_skip": reporter(("Skip",)),
            "reporter_mixed_skip": reporter((GROUNDS[0], "Prefer not to answer")),
            "reporter_unknown": reporter(("An unmapped attribution",)),
            "valid_floor": all_never(60),
        }

        values = build(rows)[3]

        self.assertEqual(values, {"valid_floor": (-1.0, 61.0)})

    def test_any_valid_positive_eds_item_establishes_reporter(self):
        rows = {
            "partial_positive": {
                "eds_4": (47, ("At least once a week",)),
                "followup": (48, (GROUNDS[2],)),
            },
            "positive_plus_missing": {
                "eds_1": (49, ("Almost every day",)),
                "eds_2": (49, ("Prefer not to answer",)),
                "followup": (50, (GROUNDS[3],)),
            },
        }

        values = build(rows)[3]

        self.assertEqual(values["partial_positive"], (1.0, 48.0))
        self.assertEqual(values["positive_plus_missing"], (1.0, 50.0))

    def test_exact_filtered_run_routes_to_population_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {P.EDS_BREADTH_PHENO_ID},
                prefixes=("num_", "ord_"),
            )
        )


if __name__ == "__main__":
    unittest.main()
