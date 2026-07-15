"""Registry for the eight approved ordinal/quantitative survey GWAS.

The registry reserves output identity independently of the builders.  A prior
manifest may reuse one of these IDs only when its construction ID is the same;
otherwise the pipeline stops before phenotype construction.
"""

from __future__ import annotations

import json
import math
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class NewQuantitativeDefinition:
    phenotype_id: str
    construction_id: str
    trait_type: str


DEFINITIONS = (
    NewQuantitativeDefinition(
        "num_overallhealth_averagepain7days",
        "overallhealth_pain_slider_numeric_v1",
        "numeric",
    ),
    NewQuantitativeDefinition(
        "ord_mania_symptom_count_pop",
        "mania_symptom_count_population_zero_v1",
        "derived_psych",
    ),
    NewQuantitativeDefinition(
        "mhq_pcl_symptom_sum_pop",
        "pcl_symptom_sum_population_zero_v1",
        "derived_psych",
    ),
    NewQuantitativeDefinition(
        "ord_mhqukb_22_pop",
        "mhq_depression_followup_population_zero_v1",
        "ordinal",
    ),
    NewQuantitativeDefinition(
        "ord_mhqukb_46_pop",
        "mania_followup_population_zero_v1",
        "ordinal",
    ),
    NewQuantitativeDefinition(
        "ord_mhqukb_54_pop",
        "psychosis_distress_population_zero_v1",
        "ordinal",
    ),
    NewQuantitativeDefinition(
        "ord_sdoh_eds_discrimination_breadth_pop",
        "eds_attribution_breadth_population_zero_v1",
        "derived",
    ),
    NewQuantitativeDefinition(
        "dim_hypomania",
        "dim_hypomania_severity_hurdle_v1",
        "derived_psych",
    ),
)

# Exact survey inputs needed to construct the eight definitions.  This is also
# the response-ingest allowlist used by memory-bounded targeted reruns.
SOURCE_QIDS = frozenset({
    # Overall Health average-pain slider.
    "1585747",
    # Depression screen and impairment.
    "1704045", "1703998", "1704022",
    # Mania screen, symptom checklist, duration, and impairment.
    "1703923", "1703928", "1703906", "1703879", "1703894",
    # Four psychotic-experience screeners and distress.
    "1703885", "1703901", "1703915", "1703871", "1703890",
    # Eleven lifetime-trauma gates and five primary PCL items.
    "1703988", "1704021", "1703976", "1704055", "1704035", "1703994",
    "1704011", "1704029", "1704020", "1704056", "1703989",
    "1704051", "1704009", "1704019", "1703997", "1703971",
    # Nine SDOH EDS items and attribution checkbox.
    "40192466", "40192489", "40192416", "40192490", "40192380",
    "40192395", "40192496", "40192519", "40192451", "40192428",
})


def phenotype_ids() -> set[str]:
    return {definition.phenotype_id for definition in DEFINITIONS}


def construction_ids_by_phenotype() -> dict[str, str]:
    return {
        definition.phenotype_id: definition.construction_id
        for definition in DEFINITIONS
    }


def write_construction_qc(
    path: Path | str,
    *,
    phenotype_id: str,
    trait_type: str,
    kind: str,
    values: Mapping[str, tuple[float, float]],
    metadata: Mapping[str, object],
    codebook_fingerprints: Mapping[str, str],
    git_state: Mapping[str, object],
) -> Path:
    """Write auditable prefilter raw-score provenance for one reserved trait."""
    expected = construction_ids_by_phenotype().get(phenotype_id)
    if expected is None:
        raise ValueError(f"Unregistered new phenotype: {phenotype_id}")
    if metadata.get("construction_id") != expected:
        raise ValueError(f"Construction mismatch for {phenotype_id}")

    raw = [
        float(score)
        for score, _age in values.values()
        if math.isfinite(float(score))
    ]
    counts = Counter(raw)
    payload = {
        "phenotype_id": phenotype_id,
        "construction_id": expected,
        "trait_type": trait_type,
        "kind": kind,
        "covar_mode": metadata.get("covar_mode", "full"),
        "question_concept_id": metadata.get("question_concept_id", ""),
        "item_concept": metadata.get("item_concept", ""),
        "question": metadata.get("question", ""),
        "answer_mapping_and_missing_policy": metadata.get("answer", ""),
        "ordinal_rule": metadata.get("ordinal_rule", ""),
        "sensitive_topics": metadata.get("sensitive_topics", ""),
        "additional_provenance": {
            key: metadata[key]
            for key in (
                "component_standardization",
                "component_names",
                "component_means",
                "component_sds",
                "degenerate_sd_policy",
                "degenerate_sd_components",
                "n_complete_case_endorsers",
            )
            if key in metadata
        },
        "raw_moments": {
            "minimum": min(raw) if raw else None,
            "maximum": max(raw) if raw else None,
            "mean": statistics.fmean(raw) if raw else None,
            "standard_deviation": statistics.pstdev(raw) if raw else None,
        },
        "raw_histogram": {
            format(value, ".17g"): count
            for value, count in sorted(counts.items())
        },
        "prefilter_scored_n": len(raw),
        "final_filtered_n": None,
        "codebook_fingerprints_sha256": dict(codebook_fingerprints),
        "pipeline_git": dict(git_state),
    }
    qc_path = Path(path)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return qc_path


def finalize_construction_qc(path: Path | str | None, final_filtered_n: int) -> None:
    if not path:
        return
    qc_path = Path(path)
    if not qc_path.is_file():
        return
    payload = json.loads(qc_path.read_text())
    payload["final_filtered_n"] = int(final_filtered_n)
    qc_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
