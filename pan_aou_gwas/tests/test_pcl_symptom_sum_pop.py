import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


TRAUMA_CODES = [*(f"mhqukb_{i}" for i in range(34, 43)), "cidi5_4", "cidi5_5"]
PCL_CODES = [f"pcl_{i}" for i in range(1, 6)]


def make_questions(records):
    qid_by_item = {code: f"qid_{code}" for code in [*TRAUMA_CODES, *PCL_CODES]}
    questions = {
        qid: {"question": code, "responses": {}}
        for code, qid in qid_by_item.items()
    }
    for iid, code, answer, age in records:
        questions[qid_by_item[code]]["responses"][iid] = (float(age), (answer,))
    item_labels = {code: code for code in qid_by_item}
    return questions, item_labels, qid_by_item


class PclSymptomSumPopulationTest(unittest.TestCase):
    def test_complete_case_hurdle_scoring_and_metadata(self):
        records = []

        # A complete all-Never trauma battery earns the distinct -1 floor.
        for code in TRAUMA_CODES:
            records.append(("no_trauma", code, "Never", 40))

        # One lifetime trauma endorsement is sufficient to open the gate; the
        # other trauma items need not be observed once an endorsement exists.
        records.append(("zero", "mhqukb_34", "Yes, within the last 12 months", 50))
        records.append(("mid", "cidi5_4", "Yes, but not in the last 12 months", 60))
        records.append(("maximum", "cidi5_5", "Yes, within the last 12 months", 70))
        for code, answer in zip(
            PCL_CODES,
            ["Not at all", "A little bit", "Moderately", "Quite a bit", "Extremely"],
        ):
            records.append(("mid", code, answer, 61))
        for code in PCL_CODES:
            records.append(("zero", code, "Not at all", 51))
            records.append(("maximum", code, "Extremely", 71))

        # Endorsers missing a PCL component are not assigned a partial sum.
        records.append(("partial_pcl", "mhqukb_34", "Yes, within the last 12 months", 80))
        for code in PCL_CODES[:-1]:
            records.append(("partial_pcl", code, "Moderately", 81))

        # Missing/invalid trauma data must not be interpreted as no trauma.
        for code in TRAUMA_CODES[:-1]:
            records.append(("partial_never", code, "Never", 90))
        for code in TRAUMA_CODES:
            answer = "Prefer not to answer" if code == "mhqukb_40" else "Never"
            records.append(("nonsubstantive", code, answer, 100))

        questions, item_labels, qid_by_item = make_questions(records)
        built = {
            pheno_id: (trait_type, kind, values, metadata)
            for pheno_id, trait_type, kind, values, metadata in P.build_derived_psych_phenotypes(
                questions, item_labels, qid_by_item
            )
        }
        trait_type, kind, values, metadata = built["mhq_pcl_symptom_sum_pop"]

        self.assertEqual((trait_type, kind), ("derived_psych", "quant"))
        self.assertEqual(values["no_trauma"], (-1.0, 40.0))
        self.assertEqual(values["zero"], (0.0, 51.0))
        self.assertEqual(values["mid"], (10.0, 61.0))
        self.assertEqual(values["maximum"], (20.0, 71.0))
        self.assertNotIn("partial_pcl", values)
        self.assertNotIn("partial_never", values)
        self.assertNotIn("nonsubstantive", values)
        self.assertEqual(metadata["construction_id"], "pcl_symptom_sum_population_zero_v1")
        self.assertEqual(metadata["ordinal_rule"], "pcl_intensity_0_4_population_floor")
        self.assertEqual(metadata["covar_mode"], "full")
        self.assertIn("mental_health", metadata["sensitive_topics"])
        self.assertEqual(set(metadata["item_concept"].split("|")), set(TRAUMA_CODES + PCL_CODES))

    def test_exact_id_routes_to_the_derived_psych_builder(self):
        self.assertTrue(
            P.wants_phenotype_source(
                {"mhq_pcl_symptom_sum_pop"}, prefixes=("psych_", "mhq_", "bin_mania_")
            )
        )


if __name__ == "__main__":
    unittest.main()
