#!/usr/bin/env python3
"""Shared interface for explicitly approved pan-AoU composite phenotypes.

The ordinary survey builders intentionally keep their compact ``(age, answers)``
response tuples.  Approved composites additionally consume aligned provenance
sidecars populated by :func:`pan_aou_gwas.build_latest_responses` so answer
concept IDs, source surveys, and response-event timestamps can be audited.
"""

from __future__ import annotations

import hashlib
import json
import math
import statistics
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Sequence


SCORED = "scored"
MISSING = "missing"
INCOMPLETE = "incomplete"
INVALID = "invalid"
CONTRADICTION = "contradiction"


def normalize_answer(value: object) -> str:
    """Normalize answer labels without erasing meaningful words."""
    return " ".join(str(value or "").strip().lower().split())


@dataclass(frozen=True)
class AnswerObservation:
    text: str
    concept_id: str = ""
    source_survey: str = ""
    response_timestamp: str = ""


@dataclass(frozen=True)
class SelectedResponse:
    qid: str
    age: float | None
    answers: tuple[AnswerObservation, ...]


@dataclass(frozen=True)
class SourceProvenance:
    qid: str
    source_survey: str
    response_timestamp: str
    answer_concept_ids: tuple[str, ...]


@dataclass(frozen=True)
class CompositeScore:
    raw_score: int | None
    observed_component_count: int
    status: str
    age: float | None
    source_provenance: tuple[SourceProvenance, ...] = ()


CompositeScorer = Callable[[Mapping[str, dict], str], CompositeScore]


@dataclass(frozen=True)
class CompositeDefinition:
    phenotype_id: str
    construction_id: str
    trait_type: str
    description: str
    source_surveys: tuple[str, ...]
    source_qids: tuple[str, ...]
    item_mappings: tuple[Mapping[str, object], ...]
    answer_mapping: Mapping[str, object]
    missing_policy: str
    valid_range: tuple[int, int]
    scorer: CompositeScorer = field(repr=False, compare=False)
    sensitive_topics: tuple[str, ...] = ()
    interpretation: str = ""
    limitations: tuple[str, ...] = ()
    kind: str = "quant"
    covar_mode: str = "full"


_REGISTRY: dict[str, CompositeDefinition] = {}


def register_composite(definition: CompositeDefinition) -> CompositeDefinition:
    """Register one approved definition and reject identity collisions."""
    if definition.kind != "quant":
        raise ValueError(f"{definition.phenotype_id}: approved composites must use kind=quant")
    if not definition.phenotype_id or not definition.construction_id:
        raise ValueError("Approved composites require phenotype and construction IDs")
    if definition.phenotype_id in _REGISTRY:
        raise ValueError(f"Duplicate approved composite phenotype ID: {definition.phenotype_id}")
    for existing in _REGISTRY.values():
        if existing.construction_id == definition.construction_id:
            raise ValueError(
                f"Duplicate approved composite construction ID: {definition.construction_id}"
            )
    _REGISTRY[definition.phenotype_id] = definition
    return definition


def registered_composites() -> tuple[CompositeDefinition, ...]:
    return tuple(_REGISTRY[key] for key in sorted(_REGISTRY))


def approved_phenotype_ids() -> frozenset[str]:
    return frozenset(_REGISTRY)


def approved_source_qids() -> frozenset[str]:
    return frozenset(qid for definition in _REGISTRY.values() for qid in definition.source_qids)


def selected_response(questions: Mapping[str, dict], qid: str, iid: str) -> SelectedResponse | None:
    """Adapt a legacy response tuple plus its optional aligned provenance sidecar."""
    question = questions.get(str(qid))
    if not question:
        return None
    response = question.get("responses", {}).get(iid)
    if response is None:
        return None
    age, texts = response
    texts = tuple(str(text) for text in texts)
    sidecar = question.get("response_sidecars", {}).get(iid, {})
    concept_ids = tuple(str(value or "") for value in sidecar.get("answer_concept_ids", ()))
    surveys = tuple(str(value or "") for value in sidecar.get("source_surveys", ()))
    timestamps = tuple(str(value or "") for value in sidecar.get("response_timestamps", ()))
    if concept_ids and len(concept_ids) != len(texts):
        raise ValueError(f"{qid}/{iid}: answer-concept sidecar is not aligned")
    if surveys and len(surveys) != len(texts):
        raise ValueError(f"{qid}/{iid}: source-survey sidecar is not aligned")
    if timestamps and len(timestamps) != len(texts):
        raise ValueError(f"{qid}/{iid}: timestamp sidecar is not aligned")
    observations = tuple(
        AnswerObservation(
            text=text,
            concept_id=concept_ids[index] if concept_ids else "",
            source_survey=surveys[index] if surveys else "",
            response_timestamp=timestamps[index] if timestamps else "",
        )
        for index, text in enumerate(texts)
    )
    try:
        parsed_age = float(age)
    except (TypeError, ValueError):
        parsed_age = None
    if parsed_age is not None and not math.isfinite(parsed_age):
        parsed_age = None
    return SelectedResponse(str(qid), parsed_age, observations)


def response_provenance(response: SelectedResponse | None) -> SourceProvenance | None:
    if response is None:
        return None
    surveys = sorted({answer.source_survey for answer in response.answers if answer.source_survey})
    timestamps = sorted(
        {answer.response_timestamp for answer in response.answers if answer.response_timestamp}
    )
    return SourceProvenance(
        qid=response.qid,
        source_survey="|".join(surveys),
        response_timestamp="|".join(timestamps),
        answer_concept_ids=tuple(answer.concept_id for answer in response.answers),
    )


def finite_mean(values: Sequence[float | None]) -> float | None:
    finite = [float(value) for value in values if value is not None and math.isfinite(value)]
    return sum(finite) / len(finite) if finite else None


def resolve_answer(
    observation: AnswerObservation,
    concept_values: Mapping[str, object],
    text_values: Mapping[str, object],
) -> tuple[object | None, str]:
    """Resolve concept-first, using text only when the concept ID is absent.

    A known concept is authoritative unless its text independently maps to a
    different value.  A present but unknown concept is invalid even if its text
    looks familiar.  Text fallback must be an explicitly unique normalized
    alias; callers express ambiguous aliases as a sequence/set of values.
    """
    concept_id = str(observation.concept_id or "").strip()
    text_key = normalize_answer(observation.text)
    if concept_id:
        if concept_id not in concept_values:
            return None, INVALID
        value = concept_values[concept_id]
        if text_key in text_values:
            text_value = text_values[text_key]
            if isinstance(text_value, (set, frozenset, list, tuple)):
                values = set(text_value)
                if len(values) == 1:
                    text_value = next(iter(values))
                else:
                    return None, INVALID
            if text_value != value:
                return None, INVALID
        return value, SCORED

    if not text_key or text_key not in text_values:
        return None, INVALID
    value = text_values[text_key]
    if isinstance(value, (set, frozenset, list, tuple)):
        values = set(value)
        if len(values) != 1:
            return None, INVALID
        value = next(iter(values))
    return value, SCORED


def source_summary(provenance: Sequence[SourceProvenance]) -> tuple[str, str]:
    surveys = sorted({p.source_survey for p in provenance if p.source_survey})
    events = sorted(
        f"{p.qid}:{p.response_timestamp}"
        for p in provenance
        if p.response_timestamp
    )
    return "|".join(surveys), "|".join(events)


def file_fingerprints(paths: Sequence[Path | str]) -> dict[str, str]:
    fingerprints: dict[str, str] = {}
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file():
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        fingerprints[str(path)] = digest.hexdigest()
    return fingerprints


def pipeline_revision(repo_root: Path | str | None) -> dict[str, object]:
    if repo_root is None:
        return {"revision": "", "dirty": None}
    root = Path(repo_root)
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, capture_output=True, check=True
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=root,
                text=True,
                capture_output=True,
                check=True,
            ).stdout.strip()
        )
        return {"revision": revision, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"revision": "", "dirty": None}


def _write_participant_qc(path: Path, rows: Sequence[tuple[str, CompositeScore]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write(
            "IID\tobserved_component_count\tstatus\traw_score\t"
            "source_surveys\tresponse_events\n"
        )
        for iid, score in rows:
            surveys, events = source_summary(score.source_provenance)
            raw = "" if score.raw_score is None else str(score.raw_score)
            handle.write(
                f"{iid}\t{score.observed_component_count}\t{score.status}\t{raw}\t"
                f"{surveys}\t{events}\n"
            )


def _construction_qc(
    definition: CompositeDefinition,
    rows: Sequence[tuple[str, CompositeScore]],
    fingerprints: Mapping[str, str],
    git_state: Mapping[str, object],
) -> dict[str, object]:
    status_counts = Counter(score.status for _iid, score in rows)
    component_counts = Counter(score.observed_component_count for _iid, score in rows)
    raw_values = [
        score.raw_score
        for _iid, score in rows
        if score.status == SCORED and score.raw_score is not None
    ]
    lo, hi = definition.valid_range
    histogram = Counter(raw_values)
    return {
        "phenotype_id": definition.phenotype_id,
        "construction_id": definition.construction_id,
        "trait_type": definition.trait_type,
        "kind": definition.kind,
        "covar_mode": definition.covar_mode,
        "description": definition.description,
        "source_surveys": list(definition.source_surveys),
        "item_qid_mappings": [dict(mapping) for mapping in definition.item_mappings],
        "answer_mapping": dict(definition.answer_mapping),
        "missing_policy": definition.missing_policy,
        "valid_range": [lo, hi],
        "sensitive_topics": list(definition.sensitive_topics),
        "status_counts": dict(sorted(status_counts.items())),
        "invalid_count": status_counts.get(INVALID, 0),
        "contradiction_count": status_counts.get(CONTRADICTION, 0),
        "raw_moments": {
            "minimum": min(raw_values) if raw_values else None,
            "maximum": max(raw_values) if raw_values else None,
            "mean": statistics.fmean(raw_values) if raw_values else None,
            "standard_deviation": statistics.pstdev(raw_values) if raw_values else None,
        },
        "raw_histogram": {str(value): histogram.get(value, 0) for value in range(lo, hi + 1)},
        "observed_component_distribution": {
            str(value): count for value, count in sorted(component_counts.items())
        },
        "prefilter_scored_n": len(raw_values),
        "final_filtered_n": None,
        "codebook_fingerprints_sha256": dict(fingerprints),
        "pipeline_git": dict(git_state),
    }


def build_registered_composite_phenotypes(
    questions: Mapping[str, dict],
    phenotype_ids: set[str] | frozenset[str] | None = None,
    qc_dir: Path | str | None = None,
    codebook_fingerprints: Mapping[str, str] | None = None,
    git_state: Mapping[str, object] | None = None,
):
    """Yield all requested registered composites through the quantitative path."""
    requested = set(phenotype_ids or ())
    for definition in registered_composites():
        if requested and definition.phenotype_id not in requested:
            continue
        participant_ids = sorted(
            {
                iid
                for qid in definition.source_qids
                for iid in questions.get(qid, {}).get("responses", {})
            }
        )
        scored_rows = [(iid, definition.scorer(questions, iid)) for iid in participant_ids]
        values = {
            iid: (float(score.raw_score), score.age)
            for iid, score in scored_rows
            if score.status == SCORED
            and score.raw_score is not None
            and score.age is not None
            and math.isfinite(score.age)
        }
        lo, hi = definition.valid_range
        if any(raw < lo or raw > hi or raw != int(raw) for raw, _age in values.values()):
            raise ValueError(f"{definition.phenotype_id}: scorer emitted an invalid raw value")

        qc_path = None
        if qc_dir is not None:
            base = Path(qc_dir)
            participant_path = base / f"{definition.phenotype_id}.participants.tsv"
            qc_path = base / f"{definition.phenotype_id}.construction.json"
            _write_participant_qc(participant_path, scored_rows)
            construction = _construction_qc(
                definition,
                scored_rows,
                codebook_fingerprints or {},
                git_state or {"revision": "", "dirty": None},
            )
            qc_path.parent.mkdir(parents=True, exist_ok=True)
            qc_path.write_text(json.dumps(construction, indent=2, sort_keys=True) + "\n")

        yield definition.phenotype_id, definition.trait_type, definition.kind, values, {
            "question_concept_id": "|".join(definition.source_qids),
            "item_concept": "|".join(
                str(mapping.get("item_concept", "")) for mapping in definition.item_mappings
            ),
            "question": definition.description,
            "answer": (
                f"complete-case raw integer score; range {lo}..{hi}; "
                f"construction={definition.construction_id}"
            ),
            "ordinal_rule": "approved_complete_case_composite",
            "covar_mode": definition.covar_mode,
            "construction_id": definition.construction_id,
            "sensitive_topics": "|".join(definition.sensitive_topics),
            "_approved_composite": True,
            "_approved_composite_qc_path": str(qc_path) if qc_path else "",
        }


def finalize_construction_qc(path: Path | str | None, final_filtered_n: int) -> None:
    """Fill the post-covariate/sample-filter N after phenotype preparation."""
    if not path:
        return
    qc_path = Path(path)
    if not qc_path.is_file():
        return
    payload = json.loads(qc_path.read_text())
    payload["final_filtered_n"] = int(final_filtered_n)
    qc_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
