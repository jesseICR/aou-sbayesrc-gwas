#!/usr/bin/env python3
"""Approved count of recreational drug classes ever used in Lifestyle."""

from __future__ import annotations

from typing import Mapping

from approved_composites import (
    CONTRADICTION,
    INVALID,
    MISSING,
    SCORED,
    CompositeDefinition,
    CompositeScore,
    normalize_answer,
    register_composite,
    resolve_answer,
    response_provenance,
    selected_response,
)


QUESTION_CONCEPT_ID = "1585636"
ITEM_CONCEPT = "recreationaldruguse_whichdrugsused"
SOURCE_SURVEY = "Lifestyle"
PHENOTYPE_ID = "num_drugs_ever_used"
CONSTRUCTION_ID = "drugs_ever_used_nine_class_checkbox_v1"

# Stable answer code, answer concept ID, and current codebook display label.
# Each distinct selected option contributes one, irrespective of frequency,
# recency, quantity, or the relative risk of the class.
DRUG_CLASS_OPTIONS = (
    (
        "WhichDrugsUsed_MarijuanaUse",
        "1585637",
        "Marijuana (cannabis, pot, grass, hash, weed, etc.)",
    ),
    (
        "WhichDrugsUsed_CocaineUse",
        "1585638",
        "Cocaine (coke, crack, etc.)",
    ),
    (
        "WhichDrugsUsed_HallucinogensUse",
        "1585643",
        "Hallucinogens (LSD, acid, mushrooms, PCP, Special K, ecstasy, etc.)",
    ),
    (
        "WhichDrugsUsed_PrescriptionStimulantsUse",
        "1585639",
        (
            "Prescription stimulants for non-medical reasons (Ritalin, Concerta, "
            "Dexedrine, Adderall, diet pills, etc.)"
        ),
    ),
    (
        "WhichDrugsUsed_SedativesUse",
        "1585642",
        (
            "Sedatives or sleeping pills for non-medical reasons (Valium, Serepax, "
            "Ativan, Xanax, Librium, Rohypnol, GHB, etc.)"
        ),
    ),
    (
        "WhichDrugsUsed_PrescriptionOpioidsUse",
        "1585645",
        (
            "Prescription opioids for non-medical reasons (fentanyl, oxycodone "
            "[OxyContin, Percocet], hydrocodone [Vicodin], methadone, buprenorphine, "
            "etc.)"
        ),
    ),
    (
        "WhichDrugsUsed_MethamphetamineUse",
        "1585640",
        (
            "Other stimulants (methamphetamine, speed, crystal meth, ice, k2/spice, "
            "bath salts, etc.)"
        ),
    ),
    (
        "WhichDrugsUsed_InhalantsUse",
        "1585641",
        "Inhalants (nitrous oxide, glue, gas, paint thinner, etc.)",
    ),
    (
        "WhichDrugsUsed_StreetOpioidsUse",
        "1585644",
        "Street opioids (heroin, opium, etc.)",
    ),
)

NONE_OPTION = (
    "WhichDrugsUsed_NoneOfTheseDrugs",
    "1585648",
    "None of these drugs",
)
OTHER_OPTION = (
    "WhichDrugsUsed_OtherSpecify",
    "1585646",
    "Other (Specify)",
)

_NONE = "none_of_these_drugs"
_OTHER = "other_specify"
_NON_SUBSTANTIVE = "non_substantive"

ANSWER_CONCEPT_VALUES = {
    concept_id: answer_code
    for answer_code, concept_id, _label in DRUG_CLASS_OPTIONS
}
ANSWER_CONCEPT_VALUES.update(
    {
        NONE_OPTION[1]: _NONE,
        OTHER_OPTION[1]: _OTHER,
        # Latest-valid ingest normally prevents these from superseding an
        # earlier substantive event, but explicit handling keeps audited or
        # manually supplied records missing rather than turning them into 0.
        "903079": _NON_SUBSTANTIVE,  # PMI: Prefer Not To Answer
        "1332892": _NON_SUBSTANTIVE,  # Prefer not to answer
        "903096": _NON_SUBSTANTIVE,  # PMI: Skip
    }
)


def _answer_aliases() -> dict[str, str]:
    aliases: dict[str, str] = {}
    for answer_code, _concept_id, label in DRUG_CLASS_OPTIONS:
        aliases[normalize_answer(answer_code)] = answer_code
        aliases[normalize_answer(label)] = answer_code

    # Short labels occur in hand-curated extracts and are unique within this
    # answer set, so they are safe checked fallbacks when concepts are absent.
    short_labels = {
        "Marijuana": "WhichDrugsUsed_MarijuanaUse",
        "Cocaine": "WhichDrugsUsed_CocaineUse",
        "Hallucinogens": "WhichDrugsUsed_HallucinogensUse",
        "Prescription stimulants": "WhichDrugsUsed_PrescriptionStimulantsUse",
        "Sedatives": "WhichDrugsUsed_SedativesUse",
        "Prescription opioids": "WhichDrugsUsed_PrescriptionOpioidsUse",
        "Methamphetamine": "WhichDrugsUsed_MethamphetamineUse",
        "Inhalants": "WhichDrugsUsed_InhalantsUse",
        "Street opioids": "WhichDrugsUsed_StreetOpioidsUse",
    }
    aliases.update(
        {normalize_answer(label): answer_code for label, answer_code in short_labels.items()}
    )
    aliases[normalize_answer(NONE_OPTION[0])] = _NONE
    aliases[normalize_answer(NONE_OPTION[2])] = _NONE
    aliases[normalize_answer(OTHER_OPTION[0])] = _OTHER
    aliases[normalize_answer(OTHER_OPTION[2])] = _OTHER
    aliases[normalize_answer("Other specify")] = _OTHER
    for label in (
        "Prefer not to answer",
        "PMI: Prefer Not To Answer",
        "PMI_PreferNotToAnswer",
        "Don't know",
        "Not sure",
        "Skip",
        "PMI: Skip",
    ):
        aliases[normalize_answer(label)] = _NON_SUBSTANTIVE
    return aliases


ANSWER_TEXT_VALUES = _answer_aliases()
DRUG_CLASS_CODES = frozenset(
    answer_code for answer_code, _concept_id, _label in DRUG_CLASS_OPTIONS
)


def score_drugs_ever_used(questions: Mapping[str, dict], iid: str) -> CompositeScore:
    """Count distinct scored drug-class options in one checkbox response event."""
    response = selected_response(questions, QUESTION_CONCEPT_ID, iid)
    provenance = response_provenance(response)
    source_provenance = (provenance,) if provenance is not None else ()
    if response is None:
        return CompositeScore(None, 0, MISSING, None, source_provenance)
    if not response.answers:
        return CompositeScore(None, 0, MISSING, response.age, source_provenance)

    selected: set[str] = set()
    for observation in response.answers:
        value, status = resolve_answer(
            observation, ANSWER_CONCEPT_VALUES, ANSWER_TEXT_VALUES
        )
        if status == INVALID:
            # Unknown present concepts are invalid even when their label looks
            # familiar.  Without a concept, wholly unrecognized text supplies
            # no recognized checkbox selection and is therefore missing.
            invalid_status = INVALID if observation.concept_id else MISSING
            return CompositeScore(
                None, 0, invalid_status, response.age, source_provenance
            )
        if value == _NON_SUBSTANTIVE:
            return CompositeScore(None, 0, MISSING, response.age, source_provenance)
        selected.add(str(value))

    selected_classes = selected & DRUG_CLASS_CODES
    if _NONE in selected and selected_classes:
        return CompositeScore(None, 9, CONTRADICTION, response.age, source_provenance)
    if _NONE in selected:
        # Other is deliberately ignored, so None plus Other remains the
        # explicit observed-zero response specified by the protocol.
        return CompositeScore(0, 9, SCORED, response.age, source_provenance)
    if selected_classes:
        return CompositeScore(
            len(selected_classes), 9, SCORED, response.age, source_provenance
        )

    # Other-only and empty/unrecognized response events do not establish that
    # none of the nine listed classes was used.
    return CompositeScore(None, 0, MISSING, response.age, source_provenance)


_ITEM_MAPPINGS = tuple(
    {
        "question_concept_id": QUESTION_CONCEPT_ID,
        "item_concept": ITEM_CONCEPT,
        "answer_code": answer_code,
        "answer_concept_id": answer_concept_id,
        "answer": label,
        "score": 1,
    }
    for answer_code, answer_concept_id, label in DRUG_CLASS_OPTIONS
) + (
    {
        "question_concept_id": QUESTION_CONCEPT_ID,
        "item_concept": ITEM_CONCEPT,
        "answer_code": NONE_OPTION[0],
        "answer_concept_id": NONE_OPTION[1],
        "answer": NONE_OPTION[2],
        "score": 0,
        "zero_anchor": True,
    },
    {
        "question_concept_id": QUESTION_CONCEPT_ID,
        "item_concept": ITEM_CONCEPT,
        "answer_code": OTHER_OPTION[0],
        "answer_concept_id": OTHER_OPTION[1],
        "answer": OTHER_OPTION[2],
        "score": "ignored",
    },
)

_DEFINITION_ANSWER_MAPPING = {
    f"{concept_id} / {answer_code} / {label}": 1
    for answer_code, concept_id, label in DRUG_CLASS_OPTIONS
}
_DEFINITION_ANSWER_MAPPING.update(
    {
        f"{NONE_OPTION[1]} / {NONE_OPTION[0]} / {NONE_OPTION[2]}": 0,
        f"{OTHER_OPTION[1]} / {OTHER_OPTION[0]} / {OTHER_OPTION[2]}": "ignored",
        "None plus any scored class": "invalid contradiction",
        "Other only, absent, empty, or non-substantive response": "missing",
    }
)


DRUGS_EVER_USED_DEFINITION = register_composite(
    CompositeDefinition(
        phenotype_id=PHENOTYPE_ID,
        construction_id=CONSTRUCTION_ID,
        trait_type="composite",
        kind="quant",
        covar_mode="full",
        description=(
            "Pan-AoU-derived count of nine recreational drug classes reported ever "
            "used in one Lifestyle checkbox question."
        ),
        source_surveys=(SOURCE_SURVEY,),
        source_qids=(QUESTION_CONCEPT_ID,),
        item_mappings=_ITEM_MAPPINGS,
        answer_mapping=_DEFINITION_ANSWER_MAPPING,
        missing_policy=(
            "One substantive checkbox response event is required. None alone anchors "
            "zero; Other is ignored but Other-only is missing; absence, empty selection, "
            "and non-substantive answers are missing; unknown concepts, concept/text "
            "conflicts, and None-plus-scored-class sets are invalid and unscored."
        ),
        valid_range=(0, 9),
        scorer=score_drugs_ever_used,
        interpretation=(
            "Higher values indicate broader lifetime exposure across the nine listed "
            "drug classes, not severity, frequency, quantity, recency, or dependence."
        ),
        limitations=(
            "Breadth is not severity or frequency: one lifetime trial and frequent use each contribute one.",
            "Survey classes have unequal prevalence and risk but receive equal weight.",
            "Lifetime use is self-reported and may be affected by recall and stigma.",
            "Other (Specify) is excluded, so unlisted-only use is not counted and Other-only is missing.",
        ),
    )
)

DEFINITION = DRUGS_EVER_USED_DEFINITION
