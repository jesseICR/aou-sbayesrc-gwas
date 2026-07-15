"""Population-referenced binaries for gated The Basics follow-up questions.

The ordinary survey one-vs-rest builder uses only respondents to the source
question as its denominator.  That is not the intended denominator for these
two follow-ups: a substantive non-gate response to the parent question proves
that the participant is not a case for any gated follow-up subtype and is
therefore a valid control.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import ordinal_rules as R


@dataclass(frozen=True)
class GatedFollowupSpec:
    parent_qid: str
    followup_qid: str
    item_concept: str
    parent_gate_answers: frozenset[str]
    parent_control_answers: frozenset[str]
    construction_id: str
    description: str

    @property
    def phenotype_prefix(self) -> str:
        return f"bin_{slug(self.item_concept)}__"


SEXUALITY_SPEC = GatedFollowupSpec(
    parent_qid="1585899",
    followup_qid="1585357",
    item_concept="genderidentity_sexualitycloserdescription",
    parent_gate_answers=frozenset({"none"}),
    parent_control_answers=frozenset({"straight", "bisexual", "gay", "lesbian"}),
    construction_id="sexuality_closer_description_expanded_parent_controls_v3",
    description=(
        "Sexuality closer-description subtype; controls include other valid "
        "closer-description options plus substantive non-None sexual-orientation responders"
    ),
)

CURRENT_LIVING_SPEC = GatedFollowupSpec(
    parent_qid="1585370",
    followup_qid="1585402",
    item_concept="livingsituation_currentliving",
    parent_gate_answers=frozenset({"other arrangement"}),
    parent_control_answers=frozenset({"own", "rent"}),
    construction_id="current_living_expanded_parent_controls_v2",
    description=(
        "Current-living subtype; controls include other valid living subtypes "
        "plus substantive Own/Rent parent-question responders"
    ),
)

SPECS = (SEXUALITY_SPEC, CURRENT_LIVING_SPEC)
SOURCE_QIDS = frozenset(
    qid for spec in SPECS for qid in (spec.parent_qid, spec.followup_qid)
)
FOLLOWUP_QIDS = frozenset(spec.followup_qid for spec in SPECS)


def slug(text: str) -> str:
    """Match the pipeline's stable, 40-character phenotype-ID slugging."""
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")[:40] or "x"


def answer_tail(text: str) -> str:
    text = str(text or "").strip()
    if ":" in text:
        _prefix, tail = text.split(":", 1)
        if tail.strip():
            return tail.strip()
    return text


def answer_norm(text: str) -> str:
    return R.norm(answer_tail(text))


def answer_slug(text: str) -> str:
    return slug(answer_tail(text))


PHENOTYPE_PREFIXES = tuple(spec.phenotype_prefix for spec in SPECS)


def is_missing_answer(text: str) -> bool:
    return R.is_missing(text) or R.is_missing(answer_tail(text))


def _single_substantive_response(question: dict, participant_id: str):
    response = question.get("responses", {}).get(participant_id)
    if response is None:
        return None
    age, answers = response
    substantive = [answer for answer in answers if not is_missing_answer(answer)]
    if len(substantive) != 1:
        return None
    return answer_norm(substantive[0]), substantive[0], age


def _substantive_response_set(question: dict, participant_id: str):
    """Return every substantive option selected at the latest response event.

    The sexuality closer-description question is select-all-that-apply.  Its
    option-level binaries therefore must retain participants who selected more
    than one valid description: the selected options are cases and every other
    option is a control.  Current-living responses are normally single-select,
    but the same representation is harmless and keeps the shared builder exact.
    """
    response = question.get("responses", {}).get(participant_id)
    if response is None:
        return None
    age, answers = response
    substantive = [answer for answer in answers if not is_missing_answer(answer)]
    if not substantive:
        return None
    normalized = frozenset(answer_norm(answer) for answer in substantive)
    return normalized, tuple(substantive), age


def _finite_age(age) -> bool:
    try:
        return age is not None and not math.isnan(float(age))
    except (TypeError, ValueError):
        return False


def build_gated_followup_binary_phenotypes(questions, specs=SPECS):
    """Yield repaired one-vs-rest binaries with expanded parent controls.

    A valid follow-up response takes precedence and supplies the case/control
    status and age.  When no valid follow-up exists, only a recognized
    non-gate parent answer supplies a control.  Gate-answer participants without
    a substantive follow-up are excluded because their subtype is unknown.
    """
    for spec in specs:
        parent = questions.get(spec.parent_qid, {})
        followup = questions.get(spec.followup_qid, {})
        if not followup:
            continue

        canonical_by_norm: dict[str, str] = {}
        for participant_id in followup.get("responses", {}):
            got = _substantive_response_set(followup, participant_id)
            if got is None:
                continue
            normalized, raw_answers, _age = got
            for option_norm, raw in zip(
                (answer_norm(answer) for answer in raw_answers), raw_answers
            ):
                canonical_by_norm.setdefault(option_norm, raw)

        participant_ids = set(parent.get("responses", {})) | set(
            followup.get("responses", {})
        )
        for target_norm in sorted(canonical_by_norm):
            target_raw = canonical_by_norm[target_norm]
            values = {}
            for participant_id in participant_ids:
                followup_value = _substantive_response_set(
                    followup, participant_id
                )
                if followup_value is not None:
                    observed_norms, _observed_raw, age = followup_value
                    if _finite_age(age):
                        values[participant_id] = (
                            1.0 if target_norm in observed_norms else 0.0,
                            float(age),
                        )
                    continue

                parent_value = _single_substantive_response(parent, participant_id)
                if parent_value is None:
                    continue
                parent_norm, _parent_raw, age = parent_value
                if parent_norm in spec.parent_control_answers and _finite_age(age):
                    values[participant_id] = (0.0, float(age))
                # Gate responders without a usable follow-up, unexpected parent
                # concepts, and all missing/non-substantive answers stay absent.

            phenotype_id = f"{spec.phenotype_prefix}{answer_slug(target_raw)}"
            yield phenotype_id, "binary", "binary", values, {
                "question_concept_id": f"{spec.parent_qid}|{spec.followup_qid}",
                "item_concept": spec.item_concept,
                "question": followup.get("question", spec.description),
                "answer": target_raw,
                "ordinal_rule": "",
                "covar_mode": "full",
                "construction_id": spec.construction_id,
                "control_definition": (
                    "other substantive follow-up options plus parent non-gate "
                    f"answers: {', '.join(sorted(spec.parent_control_answers))}; "
                    "gate responders without substantive follow-up excluded"
                ),
            }
