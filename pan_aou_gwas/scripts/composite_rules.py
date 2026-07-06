#!/usr/bin/env python3
"""Config for building validated composite scores from the codebook Scoring sheets.

The scoring sheets give per-item answer->value scoring but do NOT reverse-key the
reverse-worded items, so reverse handling is specified here explicitly. Each
composite is a prorated sum: mean(available item scores) x n_items, requiring at
least MIN_ITEM_FRACTION of items answered. Items are merged across survey
administrations by question text (so one GAD-7, one PSS, etc.).

BFI-2-XS is special: it yields the five Big Five domain scores (3 items each),
not a single total, using the standard BFI-2-XS domain/reverse structure.
"""

# Canonical instrument name (as it appears in composite_items_manifest.tsv) ->
# phenotype slug. Instruments not listed here are available in the manifest but
# not auto-built (see SPECSHEET; a few, e.g. Social Cohesion / Neighborhood
# Disorder, are left out because the sheet mixes item valence ambiguously).
SUM_INSTRUMENTS = {
    "GAD": "gad7_anxiety",
    "PHQ9": "phq9_depression",
    "ASRS": "asrs_adhd",
    "ACE": "ace_adversity",
    "AUDIT-C": "auditc_alcohol",
    "Alcohol Use Disorders Identification Test-Concise (AUDIT-C)": "auditc_alcohol",
    "IES": "ies_event_impact",
    "Everyday Discrimination Scale": "everyday_discrimination",
    "Daily Spritual Experience Scale Short Form": "daily_spiritual_experience",
    "Brief Health Literacy Scale": "health_literacy",
    "BRCS": "brcs_resilience",
    "MOS Social Support - Tangible Support": "social_support_tangible",
    "RAND Moss Social Support Survey Tangible Support Subscale": "social_support_tangible",
    "RAND Moss Social Support Survey": "social_support",
    "Promis Physical Health": "promis_physical_health",
    "PROMIS Mental Health": "promis_mental_health",
    "Perceived Stress Scale": "pss_perceived_stress",
    "PSS": "pss_perceived_stress",
    "UCLA Loneliness Scale": "ucla_loneliness",
    "UCLA LONELINESS SCALE": "ucla_loneliness",
}

# Reverse-worded items, matched by a substring of the (lowercased) question text.
# Applied on the item's own observed min/max: reversed = (min + max) - value.
REVERSE_TEXT_FRAGMENTS = [
    # PSS positively-worded items (4,5,7,8)
    "confident about your ability",
    "things were going your way",
    "able to control irritations",
    "on top of things",
    # UCLA ULS-8 positively-worded items (reverse for loneliness)
    "outgoing person",
    "find companionship when",
]

# BFI-2-XS Big Five domains (Behavioral Health, single administration).
# (item_code, is_reverse). Reverse on the 1-5 scale: 6 - value.
BFI2_DOMAINS = {
    "bfi2_extraversion": [("bfi2xs_1", True), ("bfi2xs_6", False), ("bfi2xs_11", False)],
    "bfi2_agreeableness": [("bfi2xs_2", False), ("bfi2xs_7", True), ("bfi2xs_12", False)],
    "bfi2_conscientiousness": [("bfi2xs_3", True), ("bfi2xs_8", True), ("bfi2xs_13", False)],
    "bfi2_neuroticism": [("bfi2xs_4", False), ("bfi2xs_9", False), ("bfi2xs_14", True)],
    "bfi2_openness": [("bfi2xs_5", False), ("bfi2xs_10", True), ("bfi2xs_15", False)],
}

MIN_ITEM_FRACTION = 0.8

# Composites the scoring sheets group with mixed item valence (or fail to score
# cleanly). Built directly from the survey items, reusing their ordinal scores,
# with opposite-valence items reverse-keyed (flipped on their own min/max).
# Format: slug -> (description, [(item_concept, reverse), ...]).
EXPLICIT_COMPOSITES = {
    "social_cohesion": (
        "Neighborhood social cohesion; higher = more cohesion.",
        [("scns_1", False), ("scns_2", False), ("scns_3", False), ("scns_4", False)],
    ),
    "neighborhood_disorder": (
        "Perceived neighborhood disorder (order items reversed); higher = more disorder.",
        [("nds_1", False), ("nds_2", False), ("nds_3", False), ("nds_4", False),
         ("nds_7", False), ("nds_8", False), ("nds_9", False), ("nds_10", False), ("nds_11", False),
         ("nds_5", True), ("nds_6", True), ("nds_12", True), ("nds_13", True)],
    ),
    "neighborhood_physical_disorder": (
        "Physical disorder subscale (order items reversed); higher = more disorder.",
        [("nds_1", False), ("nds_2", False), ("nds_3", False), ("nds_4", False),
         ("nds_5", True), ("nds_6", True)],
    ),
    "neighborhood_social_disorder": (
        "Social disorder subscale (order items reversed); higher = more disorder.",
        [("nds_7", False), ("nds_8", False), ("nds_9", False), ("nds_10", False), ("nds_11", False),
         ("nds_12", True), ("nds_13", True)],
    ),
    "neighborhood_walkability": (
        "PANES neighborhood walkability (crime-safety items reversed); higher = more walkable.",
        [("ips_2", False), ("ips_3", False), ("ips_4", False), ("ips_5", False), ("ips_6", False),
         ("ips_7", True), ("ips_16", True)],
    ),
    "hunger_vital_sign": (
        "Hunger Vital Sign food-insecurity screener; higher = more food insecurity.",
        [("hvs_1", False), ("hvs_2", False)],
    ),
}
