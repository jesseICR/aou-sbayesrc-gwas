import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import numeric_average_pain as N  # noqa: E402
import pan_aou_gwas as P  # noqa: E402


class NumericAveragePainTest(unittest.TestCase):
    def test_valid_slider_endpoints_and_intermediate_values(self):
        questions = {
            N.QUESTION_CONCEPT_ID: {
                "question": N.QUESTION,
                "responses": {
                    "zero": (40.0, ("0",)),
                    "middle": (41.0, ("5",)),
                    "ten": (42.0, ("10",)),
                },
            }
        }
        built = list(N.build_numeric_average_pain(questions))
        pheno_id, trait_type, kind, values, metadata = built[0]
        self.assertEqual(
            (pheno_id, trait_type, kind),
            (N.PHENOTYPE_ID, "numeric", "quant"),
        )
        self.assertEqual(
            values,
            {
                "zero": (0.0, 40.0),
                "middle": (5.0, 41.0),
                "ten": (10.0, 42.0),
            },
        )
        self.assertEqual(metadata["construction_id"], N.CONSTRUCTION_ID)
        self.assertEqual(metadata["covar_mode"], "full")

    def test_invalid_noninteger_out_of_range_and_artifact_are_missing(self):
        questions = {
            N.QUESTION_CONCEPT_ID: {
                "question": N.QUESTION,
                "responses": {
                    "negative": (40.0, ("-1",)),
                    "high": (40.0, ("11",)),
                    "fraction": (40.0, ("2.5",)),
                    "nan": (40.0, ("nan",)),
                    "artifact": (40.0, ("No matching concept",)),
                    "prefer_not": (40.0, ("PMI: Prefer Not To Answer",)),
                    "multiple": (40.0, ("1", "2")),
                },
            }
        }
        built = list(N.build_numeric_average_pain(questions))
        self.assertEqual(built[0][3], {})

    def test_no_source_question_emits_nothing(self):
        self.assertEqual(list(N.build_numeric_average_pain({})), [])

    def test_reserved_id_collision_uses_construction_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            metadata = outdir / "metadata"
            metadata.mkdir()
            manifest = metadata / "phenotype_manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("pheno_id", "construction_id"))
                writer.writerow((N.PHENOTYPE_ID, "unrelated_definition_v1"))
            with self.assertRaisesRegex(RuntimeError, "Reserved phenotype ID collision"):
                P.validate_reserved_phenotype_ids(outdir, [N.DEFINITION])

    def test_exact_routing_selects_new_numeric_phenotype(self):
        self.assertTrue(P.wants_phenotype_source(
            {N.PHENOTYPE_ID}, exact=(N.PHENOTYPE_ID,)
        ))
        self.assertFalse(P.wants_phenotype_source(
            {"num_other"}, exact=(N.PHENOTYPE_ID,)
        ))


if __name__ == "__main__":
    unittest.main()
