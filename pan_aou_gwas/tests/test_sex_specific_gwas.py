import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402
import sex_specific_gwas as S  # noqa: E402


def question(responses, concept_ids=None):
    out = {
        "question": "test question",
        "responses": {
            participant_id: (age, tuple(answers))
            for participant_id, (age, answers) in responses.items()
        },
        "response_sidecars": {},
    }
    for participant_id, ids in (concept_ids or {}).items():
        out["response_sidecars"][participant_id] = {
            "answer_concept_ids": tuple(ids),
            "source_surveys": tuple("The Basics" for _ in ids),
            "response_timestamps": tuple("2025-01-01" for _ in ids),
        }
    return out


def built_by_id(questions):
    return {
        result[0]: result
        for result in S.build_sex_specific_phenotypes(questions)
    }


class SexSpecificGwasTest(unittest.TestCase):
    def test_registry_has_twelve_unique_collision_reserved_ids(self):
        self.assertEqual(len(S.DEFINITIONS), 12)
        self.assertEqual(len(S.phenotype_ids()), 12)
        self.assertTrue(all(
            definition.phenotype_id.endswith(("_male", "_female"))
            for definition in S.DEFINITIONS
        ))
        self.assertTrue(all(definition.sex_filter in {"male", "female"} for definition in S.DEFINITIONS))

    def test_transgender_definitions_add_requested_closer_gender_case(self):
        questions = {
            S.GENDER_IDENTITY_QID: question({
                "generic": (30.0, ("Gender Identity: Transgender",)),
                "man": (31.0, ("Gender Identity: Man",)),
                "additional": (
                    32.0,
                    ("Gender Identity: None of these describe me, and I'd like to consider additional options",),
                ),
                "pna": (33.0, ("PMI: Prefer Not To Answer",)),
            }),
            S.CLOSER_GENDER_QID: question({
                "trans_woman": (40.0, ("Closer Gender Description: Trans woman",)),
                "trans_man": (41.0, ("Closer Gender Description: Trans man",)),
                "additional": (42.0, ("Closer Gender Description: Genderqueer",)),
            }),
        }
        built = built_by_id(questions)
        male = built["bin_gender_transgender_expanded_male"]
        female = built["bin_gender_transgender_expanded_female"]

        self.assertEqual(male[3]["generic"][0], 1.0)
        self.assertEqual(female[3]["generic"][0], 1.0)
        self.assertEqual(male[3]["trans_woman"][0], 1.0)
        self.assertEqual(female[3]["trans_woman"][0], 0.0)
        self.assertEqual(female[3]["trans_man"][0], 1.0)
        self.assertEqual(male[3]["trans_man"][0], 0.0)
        self.assertEqual(male[3]["additional"], (0.0, 37.0))
        self.assertNotIn("pna", male[3])
        self.assertEqual(male[4]["sex_filter"], "male")
        self.assertEqual(female[4]["sex_filter"], "female")
        self.assertEqual(male[4]["covar_mode"], "agepc")

    def test_concept_text_disagreement_is_missing(self):
        questions = {
            S.SEXUAL_ORIENTATION_QID: question(
                {
                    "gay": (30.0, ("Sexual Orientation: Gay",)),
                    "conflict": (31.0, ("Sexual Orientation: Lesbian",)),
                },
                {
                    "gay": ("1585901",),
                    "conflict": ("1585901",),
                },
            )
        }
        male = built_by_id(questions)["bin_thebasics_sexualorientation__gay_male"]
        self.assertEqual(male[3], {"gay": (1.0, 30.0)})

    def test_gay_and_gay_or_lesbian_controls_use_all_other_substantive_answers(self):
        questions = {
            S.SEXUAL_ORIENTATION_QID: question({
                "gay": (30.0, ("Sexual Orientation: Gay",)),
                "lesbian": (31.0, ("Sexual Orientation: Lesbian",)),
                "straight": (32.0, ("Sexual Orientation: Straight",)),
                "bisexual": (33.0, ("Sexual Orientation: Bisexual",)),
                "missing": (34.0, ("Skip",)),
            })
        }
        built = built_by_id(questions)
        male = built["bin_thebasics_sexualorientation__gay_male"][3]
        female = built["bin_thebasics_sexualorientation__gay_or_lesbian_female"][3]
        self.assertEqual(male["gay"][0], 1.0)
        self.assertEqual(male["lesbian"][0], 0.0)
        self.assertEqual(female["gay"][0], 1.0)
        self.assertEqual(female["lesbian"][0], 1.0)
        self.assertEqual(female["straight"][0], 0.0)
        self.assertEqual(female["bisexual"][0], 0.0)
        self.assertNotIn("missing", male)

    def test_active_duty_yes_no_and_missing(self):
        questions = {
            S.ACTIVE_DUTY_QID: question({
                "yes": (30.0, ("Active Duty: Yes",)),
                "no": (31.0, ("Active Duty: No",)),
                "pna": (32.0, ("Prefer not to answer",)),
            })
        }
        built = built_by_id(questions)
        values = built["bin_activeduty_activedutyservestatus__yes_female"][3]
        self.assertEqual(values, {"yes": (1.0, 30.0), "no": (0.0, 31.0)})

    def test_marital_status_is_basics_primary_with_cope_fill_in(self):
        questions = {
            S.MARITAL_PRIMARY_QID: question({
                "primary_divorced": (40.0, ("Current Marital Status: Divorced",)),
                "primary_wins": (41.0, ("Current Marital Status: Married",)),
                "primary_pna": (42.0, ("Prefer not to answer",)),
            }),
            S.MARITAL_COPE_QID: question({
                "cope_never": (50.0, ("Never married",)),
                "primary_wins": (51.0, ("Divorced",)),
                "primary_pna": (52.0, ("Divorced",)),
            }),
        }
        built = built_by_id(questions)
        divorced = built["bin_maritalstatus_currentmaritalstatus__divorced_male"]
        never = built["bin_maritalstatus_currentmaritalstatus__never_married_female"]
        self.assertEqual(divorced[3]["primary_divorced"], (1.0, 40.0))
        self.assertEqual(divorced[3]["primary_wins"], (0.0, 41.0))
        self.assertEqual(divorced[3]["primary_pna"], (1.0, 52.0))
        self.assertEqual(never[3]["cope_never"], (1.0, 50.0))
        self.assertEqual(divorced[4]["extra_covariates"]["from_cope"]["primary_pna"], 1.0)

    def test_under18_numeric_is_pooled_range_checked_and_quantitative(self):
        questions = {
            S.UNDER18_PRIMARY_QID: question({
                "primary": (40.0, ("3",)),
                "primary_wins": (41.0, ("2",)),
                "zero": (42.0, ("0",)),
                "fraction": (43.0, ("2.5",)),
            }),
            S.UNDER18_COPE_QID: question({
                "cope": (50.0, ("4",)),
                "primary_wins": (51.0, ("9",)),
            }),
        }
        built = built_by_id(questions)
        male = built["num_livingsituation_peopleunder18_male"]
        self.assertEqual(male[1:3], ("numeric", "quant"))
        self.assertEqual(male[3], {
            "primary": (3.0, 40.0),
            "primary_wins": (2.0, 41.0),
            "cope": (4.0, 50.0),
        })
        self.assertEqual(male[4]["extra_covariates"]["from_cope"]["cope"], 1.0)

    def test_exact_routing_rebuild_and_prior_manifest_collision(self):
        for phenotype_id in S.phenotype_ids():
            self.assertTrue(P.wants_phenotype_source({phenotype_id}, exact=S.phenotype_ids()))
            self.assertTrue(P.should_rebuild_existing_phenotype(phenotype_id))

        definition = S.DEFINITIONS[0]
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata"
            metadata.mkdir()
            path = metadata / "phenotype_manifest.tsv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=("pheno_id", "construction_id"), delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "pheno_id": definition.phenotype_id,
                    "construction_id": "different_definition_v1",
                })
            with self.assertRaisesRegex(RuntimeError, "Reserved phenotype ID collision"):
                P.validate_reserved_phenotype_ids(Path(tmp), S.DEFINITIONS)

    def test_construction_qc_records_sex_and_final_filtered_n(self):
        definition = S.DEFINITIONS[0]
        metadata = {
            "construction_id": definition.construction_id,
            "question_concept_id": "1585838|1585348",
            "item_concept": "gender_transgender_expanded",
            "answer": "test definition",
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "qc.json"
            S.write_construction_qc(
                path,
                phenotype_id=definition.phenotype_id,
                values={"a": (0.0, 30.0), "b": (1.0, 31.0)},
                metadata=metadata,
                codebook_fingerprints={"codebook": "sha"},
                git_state={"revision": "abc", "dirty": True},
            )
            S.finalize_construction_qc(path, 2)
            payload = json.loads(path.read_text())
            self.assertEqual(payload["sex_filter"], "male")
            self.assertEqual(payload["final_filtered_n"], 2)


if __name__ == "__main__":
    unittest.main()
