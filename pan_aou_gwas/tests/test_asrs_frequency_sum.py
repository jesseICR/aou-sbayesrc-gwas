import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import pan_aou_gwas as P  # noqa: E402


def response_table(answer_by_item, age=40.0):
    questions = {}
    for item_code, qid in P.ASRS_FREQUENCY_ITEMS:
        answer = answer_by_item.get(item_code)
        if answer is not None:
            questions[qid] = {
                "question": item_code,
                "responses": {"participant": (age, (answer,))},
            }
    return questions


class AsrsFrequencySumTest(unittest.TestCase):
    def test_score_endpoints_and_direction(self):
        self.assertEqual(P.asrs_complete_case_score([0.0] * 6), 0.0)
        self.assertEqual(P.asrs_complete_case_score([4.0] * 6), 24.0)
        self.assertEqual(P.asrs_complete_case_score([0, 1, 2, 3, 4, 4]), 14.0)

    def test_incomplete_score_is_missing(self):
        self.assertIsNone(P.asrs_complete_case_score([0.0] * 5))
        self.assertIsNone(P.asrs_complete_case_score([0, 1, 2, None, 3, 4]))

    def test_builder_emits_only_new_complete_case_phenotype(self):
        answers = {
            "asrs_1": "Never",
            "asrs_2": "Rarely",
            "asrs_3": "Sometimes",
            "asrs_4": "Often",
            "asrs_5": "Very often",
            "asrs_6": "Very often",
        }
        built = list(P.build_asrs_frequency_sum_phenotype(response_table(answers)))
        self.assertEqual(len(built), 1)
        pheno_id, trait_type, kind, values, metadata = built[0]
        self.assertEqual(pheno_id, "comp_asrs_adhd_0_24")
        self.assertEqual((trait_type, kind), ("composite", "quant"))
        self.assertEqual(values["participant"], (14.0, 40.0))
        self.assertEqual(
            metadata["construction_id"], "asrs_frequency_sum_complete_case_v1"
        )

    def test_builder_drops_each_missing_response_path(self):
        complete = {item_code: "Often" for item_code, _qid in P.ASRS_FREQUENCY_ITEMS}
        for item_code, _qid in P.ASRS_FREQUENCY_ITEMS:
            with self.subTest(item=item_code):
                answers = dict(complete)
                answers.pop(item_code)
                built = list(P.build_asrs_frequency_sum_phenotype(response_table(answers)))
                self.assertEqual(built[0][3], {})

    def test_skip_is_not_scored_and_qc_has_all_raw_levels(self):
        answers = {item_code: "Often" for item_code, _qid in P.ASRS_FREQUENCY_ITEMS}
        answers["asrs_3"] = "Skip"
        with tempfile.TemporaryDirectory() as tmp:
            qc_path = Path(tmp) / "asrs_qc.tsv"
            built = list(P.build_asrs_frequency_sum_phenotype(
                response_table(answers), qc_path
            ))
            self.assertEqual(built[0][3], {})
            with qc_path.open() as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        raw_levels = {
            int(row["metric"])
            for row in rows
            if row["section"] == "raw_score"
        }
        self.assertEqual(raw_levels, set(range(25)))


if __name__ == "__main__":
    unittest.main()
