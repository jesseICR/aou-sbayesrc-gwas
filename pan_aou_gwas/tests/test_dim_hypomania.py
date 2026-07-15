import sys
import unittest
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


QIDS = {
    "mhqukb_43": "1703923",
    "mhqukb_44": "1703928",
    "mhqukb_45": "1703906",
    "mhqukb_46": "1703879",
    "mhqukb_47": "1703894",
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
            iid: (float(age), answers if isinstance(answers, tuple) else (answers,))
            for iid, participant in rows.items()
            if (entry := participant.get(item)) is not None
            for age, answers in (entry,)
        }
        questions[qid] = {"question": item, "responses": responses}
    return questions


def endorser(age, symptom_answers, duration, impairment, irritable=False):
    return {
        "mhqukb_43": (age, "No" if irritable else "Yes"),
        "mhqukb_44": (age, "Yes" if irritable else "No"),
        "mhqukb_45": (age, symptom_answers),
        "mhqukb_46": (age, duration),
        "mhqukb_47": (age, impairment),
    }


def build_dimension(rows):
    built = P.build_derived_psych_phenotypes(
        question_table(rows),
        item_labels={},
        qid_by_item=QIDS,
    )
    return next(row for row in built if row[0] == "dim_hypomania")


class HypomaniaDimensionTest(unittest.TestCase):
    def test_floor_equal_weight_z_sum_shift_and_metadata(self):
        rows = {
            "floor": {
                "mhqukb_43": (30.0, "No"),
                "mhqukb_44": (31.0, "No"),
            },
            "low": endorser(
                40.0,
                ("None of the above",),
                "Less than 24 hours",
                "No problems",
            ),
            "middle": endorser(
                41.0,
                SYMPTOMS[:3],
                "At least four days in a row but less than a week",
                "No problems",
                irritable=True,
            ),
            "high": endorser(
                42.0,
                SYMPTOMS,
                "A week or more",
                (
                    "Needed treatment or caused problems with work, relationships, "
                    "finances, the law or other aspects of life"
                ),
            ),
        }

        pheno_id, trait_type, kind, values, metadata = build_dimension(rows)

        components = np.asarray(
            [
                (0.0, 0.5, 0.0),
                (3.0, 5.5, 0.0),
                (8.0, 10.0, 1.0),
            ]
        )
        expected = ((components - components.mean(axis=0)) / components.std(axis=0)).sum(axis=1)
        expected -= expected.min()

        self.assertEqual(pheno_id, "dim_hypomania")
        self.assertEqual((trait_type, kind), ("derived_psych", "quant"))
        self.assertEqual(values["floor"], (-1.0, 30.0))
        np.testing.assert_allclose(
            [values[iid][0] for iid in ("low", "middle", "high")],
            expected,
        )
        self.assertEqual(values["low"][0], 0.0)
        self.assertTrue(all(values[iid][0] >= 0 for iid in ("low", "middle", "high")))
        self.assertEqual(
            metadata["construction_id"], "dim_hypomania_severity_hurdle_v1"
        )
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertEqual(metadata["component_standardization"], "complete-case endorsers; population SD (ddof=0)")
        self.assertEqual(metadata["n_complete_case_endorsers"], 3)
        self.assertEqual(metadata["degenerate_sd_components"], "")
        self.assertIn("mental_health", metadata["sensitive_topics"])
        self.assertEqual(
            metadata["question_concept_id"],
            "1703923|1703928|1703906|1703879|1703894",
        )

    def test_degenerate_component_contributes_zero_without_nan(self):
        rows = {
            "low": endorser(
                50.0,
                ("None of the above",),
                "Less than 24 hours",
                "No problems",
            ),
            "high": endorser(
                51.0,
                SYMPTOMS[:2],
                "A week or more",
                "No problems",
            ),
        }

        values = build_dimension(rows)[3]
        metadata = build_dimension(rows)[4]

        self.assertEqual(values["low"][0], 0.0)
        self.assertGreater(values["high"][0], 0.0)
        self.assertTrue(all(np.isfinite(value) for value, _age in values.values()))
        self.assertEqual(metadata["degenerate_sd_components"], "impairment")
        self.assertEqual(
            metadata["degenerate_sd_policy"],
            "constant/nonfinite component contributes z=0",
        )

    def test_incomplete_or_invalid_endorser_is_missing(self):
        base = {
            "mhqukb_43": (45.0, "Yes"),
            "mhqukb_44": (45.0, "No"),
        }
        rows = {
            "ambiguous_negative": {"mhqukb_43": (40.0, "No")},
            "missing_checklist": {
                **base,
                "mhqukb_46": (45.0, "Less than 24 hours"),
                "mhqukb_47": (45.0, "No problems"),
            },
            "unknown_checklist": {
                **endorser(
                    46.0,
                    ("Unknown symptom",),
                    "Less than 24 hours",
                    "No problems",
                ),
            },
            "refused_checklist_mixed": {
                **endorser(
                    47.0,
                    (SYMPTOMS[0], "Prefer not to answer"),
                    "Less than 24 hours",
                    "No problems",
                ),
            },
            "unknown_duration": {
                **endorser(
                    48.0,
                    SYMPTOMS[:1],
                    "About two weeks",
                    "No problems",
                ),
            },
            "unknown_impairment": {
                **endorser(
                    49.0,
                    SYMPTOMS[:1],
                    "Less than 24 hours",
                    "A few problems",
                ),
            },
            "valid_floor": {
                "mhqukb_43": (52.0, "No"),
                "mhqukb_44": (52.0, "No"),
            },
            "valid_explicit_zero": endorser(
                53.0,
                ("None of the above",),
                "Less than 24 hours",
                "No problems",
            ),
        }

        values = build_dimension(rows)[3]

        self.assertEqual(set(values), {"valid_floor", "valid_explicit_zero"})
        self.assertEqual(values["valid_floor"], (-1.0, 52.0))
        self.assertEqual(values["valid_explicit_zero"], (0.0, 53.0))

    def test_exact_filtered_run_routes_to_derived_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {"dim_hypomania"},
                exact={"dim_hypomania"},
            )
        )


if __name__ == "__main__":
    unittest.main()
