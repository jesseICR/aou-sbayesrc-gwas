import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import new_quantitative_gwas as NQ  # noqa: E402
import pan_aou_gwas as P  # noqa: E402


class NewQuantitativeRegistryTest(unittest.TestCase):
    def test_registry_reserves_exactly_eight_unique_ids(self):
        expected = {
            "num_overallhealth_averagepain7days",
            "ord_mania_symptom_count_pop",
            "mhq_pcl_symptom_sum_pop",
            "ord_mhqukb_22_pop",
            "ord_mhqukb_46_pop",
            "ord_mhqukb_54_pop",
            "ord_sdoh_eds_discrimination_breadth_pop",
            "dim_hypomania",
        }
        self.assertEqual(NQ.phenotype_ids(), expected)
        self.assertEqual(len(NQ.DEFINITIONS), len(expected))
        self.assertEqual(
            len({definition.construction_id for definition in NQ.DEFINITIONS}),
            len(expected),
        )
        self.assertTrue(expected <= P.REBUILD_EXISTING_PHENO_IDS)

    def test_manifest_collision_requires_matching_construction(self):
        definition = NQ.DEFINITIONS[0]
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            metadata = outdir / "metadata"
            metadata.mkdir()
            manifest = metadata / "phenotype_manifest.tsv"
            manifest.write_text(
                "pheno_id\tconstruction_id\n"
                f"{definition.phenotype_id}\tan_existing_different_trait_v1\n"
            )
            with self.assertRaisesRegex(RuntimeError, "Reserved phenotype ID collision"):
                P.validate_reserved_phenotype_ids(outdir, NQ.DEFINITIONS)
            manifest.write_text(
                "pheno_id\tconstruction_id\n"
                f"{definition.phenotype_id}\t{definition.construction_id}\n"
            )
            P.validate_reserved_phenotype_ids(outdir, NQ.DEFINITIONS)

    def test_construction_qc_records_raw_distribution_and_final_n(self):
        definition = NQ.DEFINITIONS[1]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "qc.json"
            NQ.write_construction_qc(
                path,
                phenotype_id=definition.phenotype_id,
                trait_type=definition.trait_type,
                kind="quant",
                values={"a": (-1.0, 40.0), "b": (0.0, 41.0), "c": (0.0, 42.0)},
                metadata={
                    "construction_id": definition.construction_id,
                    "question_concept_id": "q1|q2",
                    "item_concept": "item1|item2",
                    "question": "test",
                    "answer": "-1 floor then 0..8",
                    "ordinal_rule": "test_rule",
                    "covar_mode": "full",
                },
                codebook_fingerprints={"book": "sha"},
                git_state={"revision": "abc", "dirty": True},
            )
            payload = json.loads(path.read_text())
            self.assertEqual(payload["raw_histogram"], {"-1": 1, "0": 2})
            self.assertEqual(payload["raw_moments"]["minimum"], -1.0)
            self.assertEqual(payload["raw_moments"]["maximum"], 0.0)
            self.assertIsNone(payload["final_filtered_n"])
            NQ.finalize_construction_qc(path, 2)
            self.assertEqual(json.loads(path.read_text())["final_filtered_n"], 2)


if __name__ == "__main__":
    unittest.main()
