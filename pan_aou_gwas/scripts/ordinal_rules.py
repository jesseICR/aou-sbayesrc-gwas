#!/usr/bin/env python3
"""Ordinal-scale mapping knowledge base for All of Us survey GWAS phenotypes.

This module is the single source of truth for how closed-ended survey answers
become ordinal numeric values. It is answer-text driven, matching how the
pipeline runs on the All of Us Researcher Workbench: OMOP `observation` rows
carry an answer concept whose name tail matches the codebook answer label, so
mapping on the normalised answer label works without a REDCap<->OMOP crosswalk.

There are three layers, applied in priority order per question:

  1. ITEM_OVERRIDES   -- explicit per-question answer->value dicts (keyed by the
                         REDCap `item_concept`). Used for miscellaneous scales
                         that no shared template captures, and to force a
                         specific rule where the generic template would be
                         ambiguous.
  2. TEMPLATES        -- shared answer-phrase families (Likert/frequency/etc.).
                         A question maps to a template iff every non-missing
                         answer label is present in the template (or its local
                         missing set).
  3. NOMINAL_ITEMS /  -- questions that are deliberately NOT ordinal (nominal
     NOMINAL_SIGS        categories); these get one-vs-rest binary GWAS only.

Anything single-select that matches none of the above is emitted to the flag
list for human review, so no arguably-ordinal question is silently dropped.

Direction convention (see specsheet Section "Ordinal mapping philosophy"):
higher value = more symptoms / stronger agreement / more frequent / better
self-rated health / more of the named construct. Signed change scales use
negative=less, 0=same, positive=more.
"""

from __future__ import annotations

import re


# --------------------------------------------------------------------------- #
# Normalisation
# --------------------------------------------------------------------------- #

# Global "answered but non-informative" answers -> always missing for ordinal.
GLOBAL_MISSING = {
    "",
    "skip",
    "pmi: skip",
    "pmi_skip",
    "prefer not to answer",
    "pmi: prefer not to answer",
    "pmi_prefernottoanswer",
    "pmi_prefer_not_to_answer",
    "don't know",
    "dont know",
    "pmi_dontknow",
    "pmi_don'tknow",
    "pmi_dont_know",
    "pmi_don't_know",
    "don't know/not sure",
    "don't know / not sure",
    "not sure",
    "i don't know",
    "i dont know",
    "unknown",
    "no answer",
}


def norm(label: str) -> str:
    """Normalise an answer label for template matching.

    Lowercase, collapse whitespace, drop a trailing period, take the tail after
    the last ':' (handles OMOP 'Category: Answer' concept names), and normalise
    the dont/don't apostrophe.
    """
    text = re.sub(r"<[^>]+>", " ", str(label or ""))
    text = re.sub(r"\s+", " ", text).strip()
    if ":" in text:
        text = text.split(":")[-1].strip()
    text = text.lower().rstrip(".").strip()
    text = text.replace("dont ", "don't ").replace("don t ", "don't ")
    return text


def is_missing(label: str) -> bool:
    return norm(label) in GLOBAL_MISSING


# --------------------------------------------------------------------------- #
# Shared templates: answer label -> ordinal value
# --------------------------------------------------------------------------- #
# Each template: {"map": {norm_label: value}, "local_missing": {norm_label},
#                 "confidence": "high|medium", "desc": ...}

def _t(map_, desc, confidence="high", local_missing=None):
    return {
        "map": {norm(k): v for k, v in map_.items()},
        "local_missing": {norm(x) for x in (local_missing or [])},
        "confidence": confidence,
        "desc": desc,
    }


TEMPLATES: dict[str, dict] = {
    # ---- frequency: symptom/behaviour how-often ------------------------------
    "freq_never_veryoften_0_4": _t(
        {"never": 0, "rarely": 1, "sometimes": 2, "often": 3, "very often": 4},
        "Frequency Never..Very often; higher = more frequent.",
    ),
    "freq_pss_0_4": _t(
        {"never": 0, "almost never": 1, "sometimes": 2, "fairly often": 3, "very often": 4},
        "Perceived-Stress-Scale frequency Never..Very often; higher = more frequent.",
    ),
    "freq_never_often_0_3": _t(
        {"never": 0, "rarely": 1, "sometimes": 2, "often": 3},
        "Frequency Never..Often (no 'very often'); higher = more frequent.",
    ),
    "freq_never_always_0_4": _t(
        {"never": 0, "rarely": 1, "sometimes": 2, "most of the time": 3, "always": 4},
        "Frequency Never..Always; higher = more frequent.",
    ),
    "freq_always_never_0_4": _t(
        {"always": 4, "often": 3, "sometimes": 2, "occasionally": 1, "never": 0},
        "Frequency Always..Never as worded; higher = more frequent.",
    ),
    "freq_never_always_na_0_3": _t(
        {"never": 0, "sometimes": 1, "most of the time": 2, "always": 3},
        "Frequency Never..Always with a Not-applicable option; higher = more frequent.",
        local_missing=["not applicable"],
    ),
    "freq_always_none_0_3": _t(
        {"always": 3, "most of the time": 2, "some of the time": 1, "none of the time": 0},
        "Frequency Always..None of the time; higher = more frequent, as worded.",
    ),
    "time_none_all_0_4": _t(
        {
            "none of the time": 0,
            "a little of the time": 1,
            "some of the time": 2,
            "most of the time": 3,
            "all of the time": 4,
        },
        "Proportion of time None..All; higher = more of the time.",
    ),
    "time_all_none_0_4": _t(
        {
            "all or almost all of the time": 4,
            "most of the time": 3,
            "some of the time": 2,
            "a little of the time": 1,
            "none of the time": 0,
        },
        "Proportion of time All..None as worded; higher = more of the time.",
    ),
    "phq_gad_0_3": _t(
        {
            "not at all": 0,
            "several days": 1,
            "over half the days": 2,
            "more than half the days": 2,
            "more than half of the days": 2,
            "nearly all days": 3,
            "nearly every day": 3,
        },
        "PHQ/GAD symptom frequency Not at all..Nearly every day; higher = more days.",
    ),
    "days_last5_midpoint": _t(
        {
            "none of the days (0 days)": 0.0,
            "a few days (1-2 days)": 1.5,
            "most days (3-4 days)": 3.5,
            "every day": 5.0,
        },
        "Days in last 5 (midpoints); every day = 5.",
        local_missing=["i don't know"],
    ),
    "freq_event_0_5": _t(
        {
            "never": 0,
            "less than once a year": 1,
            "a few times a year": 2,
            "a few times a month": 3,
            "at least once a week": 4,
            "almost everyday": 5,
            "almost every day": 5,
        },
        "Event frequency Never..Almost every day; higher = more frequent.",
    ),
    "freq_event_0_3": _t(
        {
            "never": 0,
            "a few times a month": 1,
            "at least once a week": 2,
            "almost everyday": 3,
            "almost every day": 3,
        },
        "Event frequency Never..Almost every day (4 levels); higher = more frequent.",
    ),
    "freq_covid_contact_0_3": _t(
        {
            "only a few times": 0,
            "1-3 times per month": 1,
            "1-5 times per week": 2,
            "daily": 3,
        },
        "Contact frequency Only a few times..Daily; higher = more frequent.",
    ),
    # ---- agreement -----------------------------------------------------------
    "agree_bfi_1_5": _t(
        {
            "disagree strongly": 1,
            "disagree a little": 2,
            "neutral; no opinion": 3,
            "agree a little": 4,
            "agree strongly": 5,
        },
        "BFI-2 agreement Disagree strongly..Agree strongly; higher = stronger agreement.",
    ),
    "agree_1_4": _t(
        {"strongly disagree": 1, "disagree": 2, "agree": 3, "strongly agree": 4},
        "Agreement without neutral midpoint; higher = stronger agreement.",
    ),
    "agree_1_5": _t(
        {
            "strongly disagree": 1,
            "disagree": 2,
            "neither agree nor disagree": 3,
            "agree": 4,
            "strongly agree": 5,
        },
        "Agreement Strongly disagree..Strongly agree; higher = stronger agreement.",
    ),
    "agree_neutral_1_5": _t(
        {
            "strongly disagree": 1,
            "disagree": 2,
            "neutral (neither agree nor disagree)": 3,
            "agree": 4,
            "strongly agree": 5,
        },
        "Agreement with explicit neutral label; higher = stronger agreement.",
    ),
    "agree_somewhat_1_4": _t(
        {
            "strongly disagree": 1,
            "somewhat disagree": 2,
            "somewhat agree": 3,
            "strongly agree": 4,
        },
        "Agreement Strongly disagree..Strongly agree (somewhat variant); higher = stronger agreement.",
        local_missing=["does not apply to my neighborhood"],
    ),
    "agree_lotr_1_5": _t(
        {
            "i disagree a lot": 1,
            "i disagree a little": 2,
            "i neither agree nor disagree": 3,
            "i agree a little": 4,
            "i agree a lot": 5,
        },
        "LOT-R optimism agreement; higher = stronger agreement as worded.",
    ),
    "describes_me_0_4": _t(
        {
            "does not describe me at all": 0,
            "does not describe me": 1,
            "neutral": 2,
            "describes me": 3,
            "describes me very well": 4,
        },
        "Self-description Does not describe..Describes very well; higher = more characteristic.",
    ),
    # ---- intensity / amount / distress --------------------------------------
    "intensity_0_4": _t(
        {
            "not at all": 0,
            "a little bit": 1,
            "moderately": 2,
            "quite a bit": 3,
            "extremely": 4,
        },
        "Intensity Not at all..Extremely; higher = more intense.",
    ),
    "amount_0_4": _t(
        {
            "not at all": 0,
            "a little": 1,
            "a moderate amount": 2,
            "very much": 3,
            "an extreme amount": 4,
        },
        "Amount Not at all..An extreme amount; higher = more.",
    ),
    "distress_0_4": _t(
        {
            "not distressing at all. it was a positive experience": 0,
            "not distressing, a neutral experience": 1,
            "a bit distressing": 2,
            "quite distressing": 3,
            "very distressing": 4,
        },
        "Distress positive..Very distressing; higher = more distress.",
    ),
    # ---- health / quality / ability -----------------------------------------
    "health_1_5": _t(
        {"excellent": 5, "very good": 4, "good": 3, "fair": 2, "poor": 1},
        "Self-rated health/quality Poor..Excellent; higher = better.",
    ),
    "ability_extent_0_4": _t(
        {"completely": 4, "mostly": 3, "moderately": 2, "a little": 1, "not at all": 0},
        "Ability/extent Not at all..Completely; higher = greater ability.",
    ),
    "severity_0_4": _t(
        {"none": 0, "mild": 1, "moderate": 2, "severe": 3, "very severe": 4},
        "Symptom severity None..Very severe; higher = worse.",
    ),
    "difficulty_0_4": _t(
        {
            "unable to do": 0,
            "with much difficulty": 1,
            "with some difficulty": 2,
            "with a little difficulty": 3,
            "without any difficulty": 4,
        },
        "PROMIS everyday-activity difficulty; higher = less difficulty / more able.",
    ),
    # ---- alcohol / substance frequency --------------------------------------
    "audit_freq_0_4": _t(
        {
            "never": 0,
            "monthly or less": 1,
            "2-4 times a month": 2,
            "2 to 4 per month": 2,
            "two to four times a month": 2,
            "2-3 times a week": 3,
            "2 to 3 per week": 3,
            "two to three times a week": 3,
            "4 or more times a week": 4,
            "4 or more per week": 4,
            "four or more times a week": 4,
        },
        "AUDIT-C drinking frequency; higher = more frequent.",
    ),
    "binge_freq_0_4": _t(
        {
            "never": 0,
            "never in the last year": 0,
            "less than monthly": 1,
            "monthly": 2,
            "weekly": 3,
            "daily or almost daily": 4,
        },
        "Binge/6+ drink frequency; higher = more frequent.",
    ),
    "subuse_lifestyle_0_4": _t(
        {
            "never": 0,
            "once or twice": 1,
            "monthly": 2,
            "weekly": 3,
            "daily or almost daily": 4,
        },
        "Substance-use frequency Never..Daily or almost daily; higher = more frequent.",
    ),
    "smoke_freq_0_2": _t(
        {"not at all": 0, "some days": 1, "every day": 2},
        "Current smoking/use frequency Not at all..Every day; higher = more frequent.",
    ),
    "current_use_0_3": _t(
        {
            "no, never": 0,
            "not currently, but in the past": 1,
            "yes, some days": 2,
            "yes, every day": 3,
        },
        "Use status No never..Yes every day; higher = more current/frequent use.",
    ),
    # ---- change (signed) -----------------------------------------------------
    "change_lessmore_signed": _t(
        {
            "less often than usual": -1,
            "the same as usual": 0,
            "more often than usual": 1,
        },
        "Signed change vs usual; negative = less, 0 = same, positive = more.",
    ),
    "change_muchmore_signed": _t(
        {
            "much less": -2,
            "a little less": -1,
            "about the same": 0,
            "a little more": 1,
            "much more": 2,
        },
        "Signed change Much less..Much more; negative = less, positive = more.",
    ),
    # ---- additional recurring scales ----------------------------------------
    "freq_never_always5_0_4": _t(
        {"never": 0, "rarely": 1, "sometimes": 2, "often": 3, "always": 4},
        "Frequency Never..Always (5 levels); higher = more frequent.",
    ),
    "recency_lifetime_0_2": _t(
        {
            "never": 0,
            "yes, but not in the last 12 months": 1,
            "yes, within the last 12 months": 2,
        },
        "Lifetime event recency; higher = happened and more recent.",
    ),
    "yes_maybe_no_0_2": _t(
        {"no": 0, "maybe": 1, "yes": 2},
        "Affirmation likelihood No/Maybe/Yes; higher = more likely.",
        confidence="medium",
    ),
    "ace_freq_0_2": _t(
        {"never": 0, "once": 1, "more than once": 2},
        "ACE-style frequency Never/Once/More than once; higher = more repeated.",
    ),
    "happiness_0_5": _t(
        {
            "extremely unhappy": 0,
            "very unhappy": 1,
            "moderately unhappy": 2,
            "moderately happy": 3,
            "very happy": 4,
            "extremely happy": 5,
        },
        "Happiness Extremely unhappy..Extremely happy; higher = happier.",
    ),
    "impact_0_3": _t(
        {"not at all": 0, "a little": 1, "somewhat": 2, "a lot": 3},
        "Impact/impairment Not at all..A lot; higher = more impact.",
    ),
    "food_insecurity_0_2": _t(
        {"never true": 0, "sometimes true": 1, "often true": 2},
        "Food-insecurity item Never/Sometimes/Often true; higher = more insecurity.",
    ),
    "likelihood_0_4": _t(
        {
            "very unlikely": 0,
            "unlikely": 1,
            "i do not know yet": 2,
            "likely": 3,
            "very likely": 4,
        },
        "Likelihood Very unlikely..Very likely; 'do not know yet' intermediate.",
    ),
    "spiritual_frequency_0_5": _t(
        {
            "never or almost never": 0,
            "once in a while": 1,
            "some days": 2,
            "most days": 3,
            "every day": 4,
            "many times a day": 5,
            "i do not believe in god (or a higher power)": 0,
            "i am not religious": 0,
        },
        "BMMRS spiritual/religious frequency; no belief coded lowest (0).",
        confidence="medium",
    ),
    "drink_count_band_midpoint": _t(
        {"1 or 2": 1.5, "3 or 4": 3.5, "5 or 6": 5.5, "7 to 9": 8.0, "10 or more": 10.0},
        "Daily drink-count bands (midpoints); open top bin -> 10.",
        confidence="medium",
    ),
    "visit_count_band_midpoint": _t(
        {
            "1": 1.0,
            "2-3": 2.5,
            "4-5": 4.5,
            "6-7": 6.5,
            "8-9": 8.5,
            "10-12": 11.0,
            "13-15": 14.0,
            "16 or more": 16.0,
        },
        "Number-of-visits bands (midpoints); top bin 16+ -> 16.",
        confidence="medium",
        local_missing=["don't know"],
    ),
}


# --------------------------------------------------------------------------- #
# Explicit per-item overrides (keyed by REDCap item_concept)
# --------------------------------------------------------------------------- #
# Each entry: {"rule": name, "confidence": .., "map": {norm_label: value},
#              "desc": ...}.  Overrides win over templates.

def _o(rule, map_, desc, confidence="high", local_missing=None):
    return {
        "rule": rule,
        "confidence": confidence,
        "desc": desc,
        "map": {norm(k): v for k, v in map_.items()},
        "local_missing": {norm(x) for x in (local_missing or [])},
    }


ITEM_OVERRIDES: dict[str, dict] = {
    # ---- education & income (ea_proxy.md) -----------------------------------
    "educationlevel_highestgrade": _o(
        "education_years_ea_proxy",
        {
            "never attended school or only attended kindergarten": 9,
            "grades 1 through 4 (primary school)": 9,
            "grades 5 through 8 (middle school)": 9,
            "grades 9 through 11 (some high school)": 10,
            "grade 12 or ged (high school graduate)": 13,
            "1 to 3 years after high school (some college, associate's degree, or technical school)": 15,
            "college 4 years or more (college graduate)": 18,
            "advanced degree (master's, doctorate, etc.)": 20,
        },
        "Approx education years (ea_proxy.md): low bins clamped to 9, GED=13, college grad=18, advanced=20.",
    ),
    # income annual band midpoints ($k). Codebook labels vary; matched on tail.
    "income_annualincome": _o(
        "income_midpoint_k",
        {
            "less than $10,000": 5.0,
            "$10,000- $24,999": 17.5,
            "$10,000 to $24,999": 17.5,
            "$25,000- $34,999": 30.0,
            "$25,000 to $34,999": 30.0,
            "$35,000- $49,999": 42.5,
            "$35,000 to $49,999": 42.5,
            "$50,000- $74,999": 62.5,
            "$50,000 to $74,999": 62.5,
            "$75,000-$99,999": 87.5,
            "$75,000 to $99,999": 87.5,
            "$100,000- $149,999": 125.0,
            "$100,000 to $149,999": 125.0,
            "$150,000- $199,999": 175.0,
            "$150,000 to $199,999": 175.0,
            "$200,000 or more": 250.0,
        },
        "Annual household income band midpoints in $k (ea_proxy.md); top code 250.",
    ),
    "alcohol_averagedailydrinkcount": _o(
        "drink_count_band_midpoint",
        {"1 or 2": 1.5, "3 or 4": 3.5, "5 or 6": 5.5, "7 to 9": 8.0, "10 or more": 10.0},
        "Typical drinks per drinking day (band midpoints); open top bin -> 10.",
        confidence="medium",
    ),
    "ukmh_j1": _o(
        "happiness_0_5",
        {
            "extremely unhappy": 0,
            "very unhappy": 1,
            "moderately unhappy": 2,
            "moderately happy": 3,
            "very happy": 4,
            "extremely happy": 5,
        },
        "General happiness in the past month; higher = happier.",
    ),
    # ---- The Basics ----------------------------------------------------------
    "livingsituation_howmanylivingyears": _o(
        "living_years_midpoint",
        {
            "less than 1 year": 0.5,
            "1-2 years": 1.5,
            "3-5 years": 4.0,
            "6-10 years": 8.0,
            "11-20 years": 15.0,
            "more than 20 years": 25.0,
        },
        "Years lived at current address (band midpoints); top code 25.",
        confidence="medium",
    ),
    # ---- Lifestyle -----------------------------------------------------------
    "alcohol_drinkfrequencypastyear": _o(
        "audit_freq_0_4",
        {
            "never": 0,
            "monthly or less": 1,
            "two to four times a month": 2,
            "two to three times a week": 3,
            "four or more times a week": 4,
        },
        "Past-year drinking frequency (AUDIT-C style); higher = more frequent.",
    ),
    "alcohol_6ormoredrinksoccurence": _o(
        "binge_freq_0_4",
        {
            "never in the last year": 0,
            "less than monthly": 1,
            "monthly": 2,
            "weekly": 3,
            "daily or almost daily": 4,
        },
        "Past-year 6+ drinks occurrence; higher = more frequent.",
    ),
    # ---- Overall Health ------------------------------------------------------
    "overallhealth_medicalformconfidence": _o(
        "confidence_1_5",
        {
            "extremely": 5,
            "quite a bit": 4,
            "somewhat": 3,
            "a little bit": 2,
            "not at all": 1,
        },
        "Confidence filling out medical forms; higher = more confident.",
    ),
    "overallhealth_everydayactivities": _o(
        "ability_extent_0_4",
        {"completely": 4, "mostly": 3, "moderately": 2, "a little": 1, "not at all": 0},
        "Extent able to carry out everyday physical activities; higher = more able.",
    ),
    "overallhealth_averagefatigue7days": _o(
        "severity_0_4",
        {"none": 0, "mild": 1, "moderate": 2, "severe": 3, "very severe": 4},
        "Average fatigue past 7 days; higher = worse.",
    ),
    "overallhealth_ovaryremovalhistory": _o(
        "ovary_removal_count_0_2",
        {
            "no": 0,
            "yes, but only one ovary or part of one ovary": 1,
            "yes sectioned": 1,
            "yes, both ovaries": 2,
            "yes both": 2,
        },
        "Ovary removal extent; higher = more removed. 'don't know one or both' -> missing.",
        confidence="medium",
        local_missing=["yes, but don't know whether one or both ovaries"],
    ),
    # ---- Healthcare Access & Utilization ------------------------------------
    "insurance_healthcarecoverage": _o(
        "better_worse_signed",
        {"better": 1, "about the same": 0, "worse": -1},
        "Coverage vs a year ago; signed better/same/worse.",
        confidence="medium",
    ),
    "healthadvice_spokentoprofessional": _o(
        "healthcare_recency_0_5",
        {
            "never": 0,
            "more than 5 years ago": 1,
            "more than 2 years, but not more than 5 years ago": 2,
            "more than 1 year, but not more than 2 years ago": 3,
            "more than 6 months, but not more than 1 year ago": 4,
            "6 months or less": 5,
        },
        "Recency of last professional contact; higher = more recent; never = 0.",
        confidence="medium",
    ),
    "cantaffordcare_worriedaboutpaying": _o(
        "worry_0_2",
        {"not at all worried": 0, "somewhat worried": 1, "very worried": 2},
        "Worry about affording care; higher = more worried.",
    ),
    "healthproviderracereligion_howimportant": _o(
        "importance_0_3",
        {
            "not important at all": 0,
            "slightly important": 1,
            "somewhat important": 2,
            "very important": 3,
        },
        "Importance of provider race/religion concordance; higher = more important.",
    ),
    # (HCAU '# of visits' items are handled by the visit_count_band_midpoint
    #  template, matched on the shared 1 / 2-3 / ... / 16 or more answer set.)
    # ---- Social Determinants of Health --------------------------------------
    "ips_1": _o(
        "housing_density_1_5",
        {
            "detached single-family housing": 1,
            "townhouses, row house, apartments, or condos of 2-3 stories": 2,
            "mix of single-family residences and townhouses, row houses, apartments, or condos": 3,
            "apartments or condos of 4-12 stories": 4,
            "apartments or condos of more than 12 stories": 5,
        },
        "Neighbourhood housing density; higher = denser housing form.",
        confidence="medium",
    ),
    # bmmrs_1/2/3 spiritual-frequency items are covered by the shared
    # spiritual_frequency_0_5 template (no per-item override needed).
    "nhs_xx": _o(
        "religious_service_frequency_0_5",
        {
            "i am not religious": 0,
            "never or almost never": 0,
            "less than once per month": 1,
            "1 to 3 times per month": 2,
            "once a week": 3,
            "more than once a week": 4,
        },
        "Religious-service attendance; 'I am not religious' -> lowest (0).",
        confidence="medium",
    ),
    "chis_1_xx": _o(
        "english_proficiency_0_3",
        {"not at all": 0, "not well": 1, "well": 2, "very well": 3},
        "Self-rated English speaking proficiency; higher = better.",
    ),
    # ---- COPE ----------------------------------------------------------------
    "copect_40_xx15_a": _o(
        "symptom_onset_month_ordinal",
        {
            "january or february 2020": 1,
            "march or april 2020": 2,
            "may or june 2020": 3,
            "july or august 2020": 4,
            "september or october 2020": 5,
            "november or december 2020": 6,
            "january or february 2021": 7,
        },
        "COVID symptom-onset 2-month bin; higher = later in the pandemic.",
        confidence="medium",
    ),
    "cdc_covid_xx_a": _o(
        "dose_count_1_2",
        {"1": 1, "2": 2},
        "Number of vaccine doses received; ordinal count.",
    ),
    "audit_c_1": _o(
        "audit_freq_0_4",
        {
            "never": 0,
            "monthly or less": 1,
            "2-4 times a month": 2,
            "2-3 times a week": 3,
            "4 or more times a week": 4,
        },
        "Past-month drinking frequency (AUDIT-C item 1); higher = more frequent.",
    ),
    "audit_c_3": _o(
        "binge_freq_0_4",
        {
            "never": 0,
            "less than monthly": 1,
            "monthly": 2,
            "weekly": 3,
            "daily or almost daily": 4,
        },
        "Past-month 6+ drinks frequency (AUDIT-C item 3); higher = more frequent.",
    ),
    "lot_r_1": _o(
        "agree_lotr_1_5",
        {
            "i disagree a lot": 1,
            "i disagree a little": 2,
            "i neither agree nor disagree": 3,
            "i agree a little": 4,
            "i agree a lot": 5,
        },
        "LOT-R optimism item; higher = stronger agreement.",
    ),
    "dmfs_xx_1": _o(
        "likelihood_0_4",
        {
            "very unlikely": 0,
            "unlikely": 1,
            "i do not know yet": 2,
            "likely": 3,
            "very likely": 4,
        },
        "Vaccination likelihood; higher = more likely; 'do not know yet' intermediate.",
    ),
    "dmfs_xx_1_additionaldose": _o(
        "likelihood_maybe_0_5",
        {
            "very unlikely": 0,
            "unlikely": 1,
            "not likely now, but maybe later": 2,
            "i do not know yet": 3,
            "likely": 4,
            "very likely": 5,
        },
        "Additional-dose likelihood with an explicit 'maybe later' step; higher = more likely.",
        confidence="medium",
    ),
    "msds_15": _o(
        "days_last5_midpoint",
        {
            "none of the days (0 days)": 0.0,
            "a few days (1-2 days)": 1.5,
            "most days (3-4 days)": 3.5,
            "every day": 5.0,
        },
        "Days in last 5 (midpoints); 'I don't know' -> missing.",
        local_missing=["i don't know"],
    ),
    "msds_16": _o(
        "social_change_signed",
        {
            "a lot less than normal": -2,
            "somewhat less than normal": -1,
            "about the same as normal": 0,
            "more than normal": 1,
            "a lot more than normal": 2,
        },
        "Signed change in social interaction vs normal.",
    ),
    "msds_18": _o(
        "hygiene_frequency_1_4",
        {"rarely": 1, "sometimes": 2, "most of the time": 3, "all of the time": 4},
        "Adherence frequency Rarely..All of the time; higher = more frequent.",
    ),
    "cdc_covid_19_9_xx25": _o(
        "testing_access_0_2",
        {
            "no": 0,
            "no, i tried and was unable to be tested": 1,
            "yes": 2,
        },
        "COVID testing access: no attempt/unable/tested; 'I don't know' -> missing.",
        confidence="medium",
        local_missing=["i don't know"],
    ),
    "rand_alp_csq_xx_1": _o(
        "remote_childcare_0_3",
        {
            "no, at care, school, or college full time": 0,
            "yes, but not at home": 1,
            "yes, at home part of the time": 2,
            "yes, at home full time": 3,
        },
        "At-home childcare/schooling disruption; higher = more at-home.",
        confidence="medium",
    ),
    # ---- Emotional Health History & Well-Being ------------------------------
    "mhqukb_8": _o(
        "duration_fraction_day_1_4",
        {
            "less than half the day": 1,
            "about half of the day": 2,
            "most of the day": 3,
            "all day long": 4,
        },
        "Fraction of day feelings lasted; higher = larger fraction.",
    ),
    "mhqukb_9": _o(
        "frequency_episode_1_3",
        {"less often": 1, "almost every day": 2, "every day": 3},
        "Episode frequency; higher = more frequent.",
    ),
    "mhqukb_21": _o(
        "duration_month_midpoint",
        {
            "less than a month": 0.5,
            "between one and three months": 2.0,
            "over three months, but less than six months": 4.5,
            "over six months, but less than 12 months": 9.0,
            "one to two years": 18.0,
            "over two years": 30.0,
        },
        "Episode duration in months (midpoints); >2 years top-coded to 30.",
        confidence="medium",
    ),
    "mhqukb_23": _o(
        "experience_problem_severity_0_2",
        {
            "does not sound like me": 0,
            "yes, but has not caused problems in relationships": 1,
            "yes, and this has caused problems in work or social relationships": 2,
        },
        "Symptom characteristic + problem severity; higher = more/worse.",
        confidence="medium",
    ),
    "mhqukb_24": _o(
        "episode_count_1_2",
        {"one": 1, "several": 2},
        "Number of qualifying episodes; higher = more.",
    ),
    "mhqukb_14": _o(
        "appetite_change_signed",
        {
            "decreased appetite": -1,
            "no changes in appetite": 0,
            "increased appetite": 1,
        },
        "Signed appetite change; 0 = none, +1 = increase, -1 = decrease.",
        confidence="medium",
    ),
    "mhqukb_15": _o(
        "weight_change_signed",
        {
            "lost weight": -1,
            "stayed about the same or was on a diet": 0,
            "both gained and lost some weight during the episode": 0,
            "gained weight": 1,
        },
        "Signed weight change; +1 = gain, -1 = loss, mixed/same = 0.",
        confidence="medium",
    ),
    # ---- Behavioral Health & Personality ------------------------------------
    "mhqukb_46": _o(
        "episode_duration_days_midpoint",
        {
            "less than 24 hours": 0.5,
            "at least a day, but less than four days": 2.5,
            "at least four days in a row but less than a week": 5.5,
            "a week or more": 10.0,
        },
        "Longest high/irritable period duration in days (midpoints); a week+ -> 10.",
        confidence="medium",
    ),
    "mhqukb_47": _o(
        "impairment_0_1",
        {
            "no problems": 0,
            "needed treatment or caused problems with work, relationships, finances, the law or other aspects of life": 1,
            "needed treatment or caused problems with work": 1,
        },
        "Problem/impairment from high periods; higher = more impairment.",
        confidence="medium",
    ),
    "mhqukb_54": _o(
        "distress_0_4",
        {
            "not distressing at all. it was a positive experience": 0,
            "not distressing, a neutral experience": 1,
            "a bit distressing": 2,
            "quite distressing": 3,
            "very distressing": 4,
        },
        "Distress from unusual experiences; higher = more distress.",
    ),
    "cidi5_19": _o(
        "count_band_midpoint",
        {
            "1 or 2": 1.5,
            "more than 2 but less than 10": 6.0,
            "10 or more times": 10.0,
        },
        "Lifetime panic-attack count band midpoints; open top bin -> 10.",
        confidence="medium",
    ),
}

# --------------------------------------------------------------------------- #
# Deliberately nominal single-select items (binary one-vs-rest only)
# --------------------------------------------------------------------------- #
NOMINAL_ITEMS = {
    # geography / identity / non-ordered category
    "thebasics_birthplace",
    "birthplace_birthplacestate",
    "birthplace_countryborn",
    "biologicalsexatbirth_sexatbirth",
    "homeown_currenthomeown",
    "livingsituation_currentliving",
    "maritalstatus_currentmaritalstatus",
    "gender_genderidentity",
    "genderidentity_sexualitycloserdescription",
    "sexualorientation_sexualorientation",
    "healthinsurance_healthinsurancetype",
    "insurance_insurancetypeupdate",
    "employment_employmentstatus",
    "yesnone_menstrualstoppedreason",
    "healthadvice_placeforhealthadvice",
    "healthadvice_whatkindofplace",
    # covid nominal
    "cdc_covid_xx_b",
    "cdc_covid_xx_b_firstdose",
    "cdc_covid_xx_b_seconddose",
    "cdc_covid_xx_b_additionaldose",
    "c19corset_59",  # blood type
    "cdc_covid_19_9_xx24",  # test positive yes/no/waiting
    "mhqukb_11",  # mood worse morning/evening/no variation
    "section_participation",  # survey-readiness admin
    "cdc_covid_xx",  # vaccine yes/no/trial
    "cdc_covid_xx_firstdose",  # vaccine yes/no/trial (Minute Survey)
    "cdc_covid_xx_type_dose3",  # booster / additional full dose / other (type)
    "dmfs_28",  # flu-plan yes/no/not sure/received-ago (mixed)
    "mhqukb_27",  # postpartum yes/no/NA
    "ace_5",  # parents separated/divorced/not-married (nominal 3rd cat)
    "overallhealth_menstrualstopped",  # menopause status categories
    "copect_50_xx19",  # duration unit selector Weeks/Months/Years (admin)
}

# Nominal answer-signatures (normalised, sorted tuple of core labels) that are
# never ordinal regardless of which item they appear on (e.g. US state lists).
NOMINAL_SIGNATURE_HINTS = [
    {"alabama", "alaska", "arizona", "california", "texas"},  # US states
    {"child", "friend", "parent or guardian", "relative", "spouse or partner"},
    {"married", "divorced", "widowed", "separated", "never married"},
    {"pfizer (pfizer-biontech)", "moderna", "other"},
    {"a", "b", "ab", "o"},
    {"weeks", "months", "years"},  # duration unit selector (admin)
    {"booster dose", "additional full dose", "other"},  # vaccine dose type
]


# --------------------------------------------------------------------------- #
# Sensitive-topic detection (flag, do not exclude)
# --------------------------------------------------------------------------- #
SENSITIVE_PATTERNS = [
    (r"sexual orientation|gay|lesbian|bisexual|straight|heterosexual|queer|asexual", "sexual_orientation"),
    (r"gender identity|transgender|two.?spirit|nonbinary", "gender_identity"),
    (r"suicid|kill(ing)? yourself|end your life|self.?harm", "suicidality_self_harm"),
    (r"sexual (intercourse|activity|partner|assault|abuse)|forced.*sex|rape", "sexual_behavior_or_trauma"),
    (r"abuse|violence|trauma|assault|hit you|hurt you|threatened", "trauma_violence"),
    (r"drug use|cannabis|marijuana|cocaine|heroin|opioid|hallucinogen|inhalant|stimulant|illicit", "substance_use"),
    (r"alcohol|drink(ing)? .*alcohol|binge", "alcohol"),
    (r"menstrua|menopaus|pregnan|ovar|hysterectomy|birth control|abortion|reproductive", "reproductive_menstrual"),
    (r"depress|anxiet|bipolar|mania|psychosis|schizophren|panic|ptsd|mental health", "mental_health"),
    (r"covid|vaccin|coronavirus|pandemic", "covid_related"),
    (r"income|poverty|food (run out|didn't last)|afford|financial", "financial_hardship"),
    (r"incarcerat|prison|jail|arrest|police|legal (trouble|problem)", "justice_involvement"),
    (r"immigra|citizenship|country.*born|birthplace", "immigration_origin"),
    (r"religio|god|spiritual|higher power|faith", "religion"),
    (r"disab|deaf|blind|difficulty (walking|dressing|concentrating)", "disability"),
]


def sensitive_topics(question_label: str, survey: str) -> list[str]:
    text = f"{survey} {question_label}".lower()
    hits = []
    for pat, tag in SENSITIVE_PATTERNS:
        if re.search(pat, text):
            hits.append(tag)
    return hits
