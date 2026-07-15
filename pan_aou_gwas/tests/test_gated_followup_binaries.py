import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import gated_followup_binaries as G  # noqa: E402
import pan_aou_gwas as P  # noqa: E402


def question(question, responses):
    return {
        "question": question,
        "responses": {
            iid: (age, tuple(answers))
            for iid, (age, answers) in responses.items()
        },
    }


def built_by_id(questions):
    return {
        result[0]: result
        for result in G.build_gated_followup_binary_phenotypes(questions)
    }


class GatedFollowupBinariesTest(unittest.TestCase):
    def test_sexuality_subtype_uses_expanded_parent_controls(self):
        questions = {
            G.SEXUALITY_SPEC.parent_qid: question(
                "Which of the following best represents how you think of yourself?",
                {
                    "case": (30.0, ("Sexual Orientation: None",)),
                    "other_followup": (31.0, ("Sexual Orientation: None",)),
                    "straight": (32.0, ("Sexual Orientation: Straight",)),
                    "bisexual": (33.0, ("Sexual Orientation: Bisexual",)),
                    "gay": (34.0, ("Sexual Orientation: Gay",)),
                    "lesbian": (35.0, ("Sexual Orientation: Lesbian",)),
                    "gate_no_followup": (36.0, ("Sexual Orientation: None",)),
                    "parent_skip": (37.0, ("Skip",)),
                },
            ),
            G.SEXUALITY_SPEC.followup_qid: question(
                "Are any of these a closer description of how you think of yourself?",
                {
                    "case": (40.0, ("Sexuality Closer Description: Queer",)),
                    "other_followup": (
                        41.0,
                        ("Sexuality Closer Description: Asexual",),
                    ),
                    "followup_prefer_not": (42.0, ("Prefer Not To Answer",)),
                    "followup_dont_know": (
                        43.0,
                        ("Sexuality Closer Description: Dont Know",),
                    ),
                },
            ),
        }

        phenotype_id = (
            "bin_genderidentity_sexualitycloserdescriptio__queer"
        )
        pheno_id, trait_type, kind, values, metadata = built_by_id(questions)[
            phenotype_id
        ]

        self.assertEqual(pheno_id, phenotype_id)
        self.assertEqual((trait_type, kind), ("binary", "binary"))
        self.assertEqual(
            values,
            {
                "case": (1.0, 40.0),
                "other_followup": (0.0, 41.0),
                "straight": (0.0, 32.0),
                "bisexual": (0.0, 33.0),
                "gay": (0.0, 34.0),
                "lesbian": (0.0, 35.0),
            },
        )
        self.assertNotIn("gate_no_followup", values)
        self.assertNotIn("parent_skip", values)
        self.assertNotIn("followup_prefer_not", values)
        self.assertNotIn("followup_dont_know", values)
        self.assertEqual(
            metadata["construction_id"],
            "sexuality_closer_description_expanded_parent_controls_v3",
        )
        self.assertEqual(metadata["question_concept_id"], "1585899|1585357")

    def test_current_living_subtype_uses_own_and_rent_as_controls(self):
        questions = {
            G.CURRENT_LIVING_SPEC.parent_qid: question(
                "Do you own or rent the place where you live?",
                {
                    "case": (50.0, ("Current Home Own: Other Arrangement",)),
                    "other_followup": (
                        51.0,
                        ("Current Home Own: Other Arrangement",),
                    ),
                    "own": (52.0, ("Current Home Own: Own",)),
                    "rent": (53.0, ("Current Home Own: Rent",)),
                    "gate_no_followup": (
                        54.0,
                        ("Current Home Own: Other Arrangement",),
                    ),
                    "parent_pna": (55.0, ("Prefer Not To Answer",)),
                },
            ),
            G.CURRENT_LIVING_SPEC.followup_qid: question(
                "Where are you currently living?",
                {
                    "case": (60.0, ("Current Living: Shelter",)),
                    "other_followup": (61.0, ("Current Living: Family",)),
                    "followup_skip": (62.0, ("Skip",)),
                },
            ),
        }

        phenotype_id = "bin_livingsituation_currentliving__shelter"
        _pheno_id, _trait_type, _kind, values, metadata = built_by_id(questions)[
            phenotype_id
        ]

        self.assertEqual(
            values,
            {
                "case": (1.0, 60.0),
                "other_followup": (0.0, 61.0),
                "own": (0.0, 52.0),
                "rent": (0.0, 53.0),
            },
        )
        self.assertNotIn("gate_no_followup", values)
        self.assertNotIn("parent_pna", values)
        self.assertNotIn("followup_skip", values)
        self.assertEqual(
            metadata["construction_id"],
            "current_living_expanded_parent_controls_v2",
        )
        self.assertEqual(metadata["question_concept_id"], "1585370|1585402")

    def test_multiselect_closer_description_retains_each_selected_case(self):
        questions = {
            G.SEXUALITY_SPEC.parent_qid: question(
                "Which of the following best represents how you think of yourself?",
                {
                    "multi": (29.0, ("Sexual Orientation: None",)),
                    "straight": (30.0, ("Sexual Orientation: Straight",)),
                },
            ),
            G.SEXUALITY_SPEC.followup_qid: question(
                "Are any of these a closer description of how you think of yourself?",
                {
                    "multi": (
                        31.0,
                        (
                            "Sexuality Closer Description: Queer",
                            "Sexuality Closer Description: Not Figured Out",
                        ),
                    ),
                },
            ),
        }

        built = built_by_id(questions)
        queer = built[
            "bin_genderidentity_sexualitycloserdescriptio__queer"
        ][3]
        not_figured_out = built[
            "bin_genderidentity_sexualitycloserdescriptio__not_figured_out"
        ][3]

        self.assertEqual(queer["multi"], (1.0, 31.0))
        self.assertEqual(not_figured_out["multi"], (1.0, 31.0))
        self.assertEqual(queer["straight"], (0.0, 30.0))
        self.assertEqual(not_figured_out["straight"], (0.0, 30.0))

    def test_exact_routing_and_all_dynamic_subtype_ids_force_matrix_rebuild(self):
        for prefix in G.PHENOTYPE_PREFIXES:
            phenotype_id = f"{prefix}test_option"
            self.assertTrue(P.wants_phenotype_source(
                {phenotype_id}, prefixes=G.PHENOTYPE_PREFIXES
            ))
            self.assertTrue(P.should_rebuild_existing_phenotype(phenotype_id))
        self.assertFalse(P.should_rebuild_existing_phenotype("bin_unrelated__yes"))
        self.assertEqual(G.FOLLOWUP_QIDS, frozenset({"1585357", "1585402"}))

    def test_versioned_construction_invalidates_stale_gwas_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            glm = root / "test.glm.linear"
            sumstats = root / "test.sumstats.tsv.gz"
            params = root / "test.gwas.params.tsv"
            glm.write_text("old glm\n")
            sumstats.write_text("old sumstats\n")
            row = {
                "pheno_id": (
                    "bin_genderidentity_sexualitycloserdescriptio__queer"
                ),
                "pheno_name": "test_resid",
                "trait_type": "binary",
                "kind": "binary",
                "n": 1000,
                "n_cases": 100,
                "n_controls": 900,
                "covar_mode": "full",
                "sex_filter": "all",
                "extra_covariates": "",
                "construction_id": "old_followup_only_controls_v1",
                "glm": str(glm),
                "sumstats": str(sumstats),
                "gwas_params": str(params),
            }
            P.write_gwas_params(params, row)
            row["construction_id"] = G.SEXUALITY_SPEC.construction_id
            self.assertFalse(P.gwas_params_match(row))
            self.assertFalse(P.output_complete(row))


if __name__ == "__main__":
    unittest.main()
