import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import approved_composites as AC  # noqa: E402
import gen_composite_docs as G  # noqa: E402
import pan_aou_gwas as P  # noqa: E402


class ApprovedCompositeInfrastructureTest(unittest.TestCase):
    def test_ingest_sidecars_align_without_changing_response_tuple(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "survey.csv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow([
                    "person_id", "survey", "question_concept_id", "question",
                    "answer_concept_id", "answer", "survey_datetime", "age_at_survey",
                ])
                writer.writerows([
                    ["p1", "Survey A", "q1", "Question", "11", "Yes", "2024-01-01", "40"],
                    ["p1", "Survey A", "q1", "Question", "903079", "PMI: Prefer Not To Answer", "2025-01-01", "41"],
                    ["p2", "Survey B", "q2", "Checkbox", "21", "Alpha", "2024-02-01", "50"],
                    ["p2", "Survey B", "q2", "Checkbox", "22", "Beta", "2024-02-01", "50"],
                    ["p3", "Survey A", "q1", "Question", "903079", "PMI: Prefer Not To Answer", "2025-01-01", "60"],
                ])
            questions = P.build_latest_responses(
                [path], {"p1", "p2", "p3"}, retain_latest_missing_qids={"q1"}
            )

        self.assertEqual(questions["q1"]["responses"]["p1"], (40.0, ("Yes",)))
        self.assertEqual(
            questions["q1"]["response_sidecars"]["p1"],
            {
                "answer_concept_ids": ("11",),
                "source_surveys": ("Survey A",),
                "response_timestamps": ("2024-01-01",),
            },
        )
        response = AC.selected_response(questions, "q2", "p2")
        self.assertEqual(tuple(answer.text for answer in response.answers), ("Alpha", "Beta"))
        self.assertEqual(tuple(answer.concept_id for answer in response.answers), ("21", "22"))
        self.assertEqual(tuple(answer.source_survey for answer in response.answers), ("Survey B", "Survey B"))
        self.assertEqual(
            questions["q1"]["responses"]["p3"],
            (60.0, ("PMI: Prefer Not To Answer",)),
        )

    def test_concept_first_resolution_and_checked_text_fallback(self):
        concepts = {"1": 1, "0": 0}
        aliases = {"yes": 1, "no": 0, "ambiguous": (0, 1)}
        self.assertEqual(
            AC.resolve_answer(AC.AnswerObservation("Yes", "1"), concepts, aliases),
            (1, AC.SCORED),
        )
        self.assertEqual(
            AC.resolve_answer(AC.AnswerObservation("Yes", "0"), concepts, aliases)[1],
            AC.INVALID,
        )
        self.assertEqual(
            AC.resolve_answer(AC.AnswerObservation("Yes"), concepts, aliases),
            (1, AC.SCORED),
        )
        self.assertEqual(
            AC.resolve_answer(AC.AnswerObservation("Yes", "unknown"), concepts, aliases)[1],
            AC.INVALID,
        )
        self.assertEqual(
            AC.resolve_answer(AC.AnswerObservation("ambiguous"), concepts, aliases)[1],
            AC.INVALID,
        )

    def test_reserved_manifest_collision_fails_loudly(self):
        definition = AC.CompositeDefinition(
            phenotype_id="reserved_id",
            construction_id="expected_v1",
            trait_type="composite",
            description="test",
            source_surveys=("Survey",),
            source_qids=("qid",),
            item_mappings=({"item_concept": "item", "question_concept_id": "qid"},),
            answer_mapping={"No": 0, "Yes": 1},
            missing_policy="complete case",
            valid_range=(0, 1),
            scorer=lambda _questions, _iid: AC.CompositeScore(None, 0, AC.MISSING, None),
        )
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            metadata = outdir / "metadata"
            metadata.mkdir()
            manifest = metadata / "phenotype_manifest.tsv"
            manifest.write_text("pheno_id\tconstruction_id\nreserved_id\tother_v1\n")
            with self.assertRaisesRegex(RuntimeError, "Reserved phenotype ID collision"):
                P.validate_reserved_phenotype_ids(outdir, [definition])
            manifest.write_text("pheno_id\tconstruction_id\nreserved_id\texpected_v1\n")
            P.validate_reserved_phenotype_ids(outdir, [definition])

    def test_construction_qc_includes_required_distributions(self):
        definition = AC.CompositeDefinition(
            phenotype_id="qc_id",
            construction_id="qc_v1",
            trait_type="composite",
            description="test",
            source_surveys=("Survey",),
            source_qids=("qid",),
            item_mappings=({"item_concept": "item", "question_concept_id": "qid"},),
            answer_mapping={"No": 0, "Yes": 1},
            missing_policy="complete case",
            valid_range=(0, 1),
            scorer=lambda _questions, _iid: AC.CompositeScore(None, 0, AC.MISSING, None),
        )
        rows = [
            ("p0", AC.CompositeScore(0, 1, AC.SCORED, 40.0)),
            ("p1", AC.CompositeScore(1, 1, AC.SCORED, 41.0)),
            ("bad", AC.CompositeScore(None, 0, AC.INVALID, None)),
            ("contra", AC.CompositeScore(None, 1, AC.CONTRADICTION, 42.0)),
        ]
        payload = AC._construction_qc(definition, rows, {"book": "sha"}, {"revision": "abc", "dirty": True})
        self.assertEqual(payload["raw_histogram"], {"0": 1, "1": 1})
        self.assertEqual(payload["invalid_count"], 1)
        self.assertEqual(payload["contradiction_count"], 1)
        self.assertIn("observed_component_distribution", payload)
        self.assertIsNone(payload["final_filtered_n"])
        self.assertEqual(payload["codebook_fingerprints_sha256"], {"book": "sha"})

    def test_registry_contains_exactly_the_five_approved_quantitative_ids(self):
        expected = {
            "comp_healthcare_discrimination",
            "comp_disability_count",
            "comp_housing_problem_count",
            "psych_psychotic_experiences_count",
            "num_drugs_ever_used",
        }
        self.assertEqual(AC.approved_phenotype_ids(), expected)
        definitions = {definition.phenotype_id: definition for definition in AC.registered_composites()}
        self.assertTrue(all(definition.kind == "quant" for definition in definitions.values()))
        self.assertEqual(definitions["psych_psychotic_experiences_count"].trait_type, "derived_psych")
        self.assertTrue(
            all(
                definition.trait_type == "composite"
                for pheno_id, definition in definitions.items()
                if pheno_id != "psych_psychotic_experiences_count"
            )
        )
        self.assertTrue(P.wants_phenotype_source({"num_drugs_ever_used"}, exact=expected))
        self.assertFalse(P.wants_phenotype_source({"height_cm"}, exact=expected))

    def test_registered_builder_writes_internal_qc_and_construction_metadata(self):
        answers = (
            "Bug infestation",
            "Mold",
        )
        questions = {
            "40192402": {
                "question": "housing",
                "responses": {"participant": (50.0, answers)},
                "response_sidecars": {
                    "participant": {
                        "answer_concept_ids": ("40192460", "40192479"),
                        "source_surveys": (
                            "Social Determinants of Health",
                            "Social Determinants of Health",
                        ),
                        "response_timestamps": ("2025-01-01", "2025-01-01"),
                    }
                },
            }
        }
        with tempfile.TemporaryDirectory() as tmp:
            qc_dir = Path(tmp) / "qc"
            built = list(AC.build_registered_composite_phenotypes(
                questions,
                phenotype_ids={"comp_housing_problem_count"},
                qc_dir=qc_dir,
                codebook_fingerprints={"book": "sha"},
                git_state={"revision": "abc", "dirty": True},
            ))
            pheno_id, trait_type, kind, values, metadata = built[0]
            self.assertEqual((pheno_id, trait_type, kind), (
                "comp_housing_problem_count", "composite", "quant"
            ))
            self.assertEqual(values, {"participant": (2.0, 50.0)})
            participant_qc = qc_dir / "comp_housing_problem_count.participants.tsv"
            construction_qc = Path(metadata["_approved_composite_qc_path"])
            self.assertTrue(participant_qc.is_file())
            self.assertTrue(construction_qc.is_file())
            AC.finalize_construction_qc(construction_qc, 1)
            payload = __import__("json").loads(construction_qc.read_text())
            self.assertEqual(payload["construction_id"], "housing_problem_count_complete_case_v1")
            self.assertEqual(payload["raw_histogram"]["2"], 1)
            self.assertEqual(payload["final_filtered_n"], 1)

    def test_generated_approved_documentation_matches_checked_in_section(self):
        generated = "\n".join(G.approved_composite_lines()).strip()
        checked_in = (SCRIPT_DIR.parent / "metadata" / "COMPOSITE_SCORES.md").read_text()
        self.assertIn(generated, checked_in)
        for definition in AC.registered_composites():
            self.assertEqual(generated.count(f"#### `{definition.phenotype_id}`"), 1)
            self.assertIn(f"`{definition.construction_id}`", generated)


if __name__ == "__main__":
    unittest.main()
