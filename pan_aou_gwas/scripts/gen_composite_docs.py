#!/usr/bin/env python3
"""Generate the validated-composite-score reference section from the manifest.

Reads metadata/composite_items_manifest.tsv and composite_rules.py and writes a
markdown section documenting each scale: item count, item questions, per-item
scoring, reverse-keyed items, and how items combine. Output is spliced into
SPECSHEET.md (section 11c).
"""

from __future__ import annotations

import csv
import re
from collections import OrderedDict
from pathlib import Path

import composite_rules as CR

ROOT = Path(__file__).resolve().parent.parent


def norm_q(t: str) -> str:
    t = re.sub(r"<[^>]+>", " ", str(t or ""))
    t = t.translate(str.maketrans({"‘": "'", "’": "'", "“": '"', "”": '"'}))
    return re.sub(r"\s+", " ", t).strip().lower().strip(" ?.\"'")


def norm_a(t: str) -> str:
    return re.sub(r"\s+", " ", str(t or "")).strip().lower()


# display title -> (manifest instrument names, built-phenotype-slug or None)
GROUPS = OrderedDict([
    ("GAD-7 — Generalized Anxiety Disorder scale (anxiety)", (["GAD"], "gad7_anxiety")),
    ("PHQ-9 — Patient Health Questionnaire (depression)", (["PHQ9"], "phq9_depression")),
    ("PSS — Perceived Stress Scale", (["PSS", "Perceived Stress Scale"], "pss_perceived_stress")),
    ("ACE — Adverse Childhood Experiences", (["ACE"], "ace_adversity")),
    ("IES — Impact of Event Scale (event-related distress)", (["IES"], "ies_event_impact")),
    ("ASRS — Adult ADHD Self-Report Scale (Part A screener)", (["ASRS"], "asrs_adhd")),
    ("UCLA / ULS-8 — Loneliness", (["UCLA LONELINESS SCALE", "UCLA Loneliness Scale"], "ucla_loneliness")),
    ("Everyday Discrimination Scale", (["Everyday Discrimination Scale"], "everyday_discrimination")),
    ("MOS Social Support (RAND) + Tangible subscale",
     (["RAND Moss Social Support Survey", "MOS Social Support - Tangible Support",
       "RAND Moss Social Support Survey Tangible Support Subscale"], "social_support")),
])

REVERSE_FRAG = CR.REVERSE_TEXT_FRAGMENTS


def main() -> None:
    rows = list(csv.DictReader(open(ROOT / "metadata/composite_items_manifest.tsv"), delimiter="\t"))

    out = []
    out.append("## 11c. Validated composite score definitions\n")
    out.append(
        "Each composite is a **prorated sum**: mean(available item scores) × n_items, "
        "requiring valid answers for more than half of items. Reverse-worded items (flagged per scale) are "
        "flipped on their own min/max before summing. Items are matched to survey responses "
        "by question text and merged across survey administrations. PHQ-9 and GAD-7 pool "
        "EHHWB and COPE administrations with EHHWB priority and a `from_cope` covariate; "
        "PSS-10 pools SDOH and COPE administrations with SDOH priority and the same source "
        "covariate. The score is then "
        "inverse-normal-transformed and residualized like any quantitative trait (§4.1). "
        "Phenotype ids are prefixed `comp_`.\n"
    )

    for title, (insts, slug) in GROUPS.items():
        items = OrderedDict()  # norm_q -> (question, {norm_answer: (label, value)})
        for r in rows:
            if r["instrument"] not in insts:
                continue
            k = norm_q(r["question"])
            items.setdefault(k, (r["question"], OrderedDict()))
            items[k][1][norm_a(r["answer_label"])] = (r["answer_label"], r["value"])
        if not items:
            continue
        n = len(items)
        scales = {tuple((v[0], v[1]) for v in amap.values()) for _q, amap in items.values()}
        revs = [q for k, (q, _a) in items.items() if any(f in k for f in REVERSE_FRAG)]

        out.append(f"### {title}\n")
        out.append(f"- **Items:** {n}")
        if len(scales) == 1:
            one = list(items.values())[0][1]
            def fmt(v):
                fv = float(v)
                return int(fv) if fv == int(fv) else fv
            out.append("- **Per-item scoring:** " + ", ".join(f"{lbl} = {fmt(val)}" for lbl, val in one.values()))
        else:
            out.append(f"- **Per-item scoring:** {len(scales)} answer scales across items (shown per item below)")
        combo = f"prorated sum of {n} items"
        combo += f"; {len(revs)} reverse-keyed" if revs else "; no reverse-keyed items"
        out.append(f"- **Total score:** {combo}")
        out.append(f"- **Auto-built:** {'yes (comp_' + slug + ')' if slug else 'no — documented only; items are also GWASed individually as ordinal phenotypes'}")
        if slug == "pss_perceived_stress":
            out.append("- **Pooling:** SDOH is primary; COPE fills COPE-only responses. The GWAS residualization includes `from_cope`.")
        out.append("- **Questions:**")
        for k, (q, amap) in items.items():
            rev = " *(reverse-keyed)*" if any(f in k for f in REVERSE_FRAG) else ""
            if len(scales) == 1:
                out.append(f"    - {q}{rev}")
            else:
                sc = ", ".join(f"{lbl}={val}" for lbl, val in amap.values())
                out.append(f"    - {q}{rev}  — [{sc}]")
        out.append("")

    # explicit mixed-valence composites, item labels from the ordinal manifest
    code_qtext, code_vals = {}, {}
    for r in csv.DictReader(open(ROOT / "metadata/ordinal_mapping_manifest.tsv"), delimiter="\t"):
        code_qtext[r["item_concept"]] = r["field_label"]
        code_vals.setdefault(r["item_concept"], set()).add(r["ordinal_value"])
    out.append("### Neighborhood, walkability & food-insecurity composites\n")
    out.append(
        "Built directly from the survey items (reusing their ordinal scores), because the scoring "
        "sheet groups these with mixed item valence. Opposite-valence items are reverse-keyed.\n"
    )
    for slug, (desc, items) in CR.EXPLICIT_COMPOSITES.items():
        present = [(c, rev) for c, rev in items if c in code_qtext]
        if len(present) < 2:
            continue
        n_rev = sum(1 for _c, rev in present if rev)
        out.append(f"#### comp_{slug}\n")
        out.append(f"- {desc}")
        out.append(f"- **Items:** {len(present)}; **reverse-keyed:** {n_rev}; prorated sum")
        out.append("- **Questions:**")
        for c, rev in present:
            out.append(f"    - {code_qtext[c]}" + (" *(reverse-keyed)*" if rev else ""))
        out.append("")

    (ROOT / "metadata/COMPOSITE_SCORES.md").write_text("\n".join(out) + "\n")
    print(f"Wrote metadata/COMPOSITE_SCORES.md ({len(out)} lines)")


if __name__ == "__main__":
    main()
