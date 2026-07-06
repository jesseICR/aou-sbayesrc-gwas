# Specsheet: Pan-UKB-style All of Us survey + physical-measurement GWAS (HapMap3, residualize-first PLINK2)

**Project:** `aou_panallofus_hapmap3_residualized_gwas`
**Date:** 2026-07-06
**Dataset:** All of Us Controlled Tier, cdr v9 (`C2025Q4R6` by default)
**Ancestry stratum:** classified European-ancestry WGS participants, unrelated (KING < 0.0441941)
**SNP panel:** HapMap3 HQ bfile already built on-platform (see §2)
**Association engine:** PLINK2 `--glm` **linear** GWAS on **pre-residualized** phenotypes, **no covariates in PLINK2**
**Downstream use:** LDSC-ready per-phenotype summary statistics for later h²/rg scans
**Deliverable of *this* run:** the GWAS summary-statistic files only. Downstream LDSC/rg (old §21) is deferred.

This file supersedes `../../draft_specsheet.md` for the current run. It bakes in the decisions made
after the covariate experiment (`~/projects/ukgwas/covariate_experiment`) and after the HapMap3 HQ
bfile was built in the AoU environment.

---

## 0. What changed from the draft specsheet

1. **HapMap3 HQ bfile is already built** in the workspace bucket; this pipeline consumes it directly
   (§2). No variant QC step is run here.
2. **Residualize-first is the method for *every* phenotype**, binary and continuous alike. Covariates
   are regressed out of the phenotype *before* PLINK2; PLINK2 then runs a covariate-free linear GWAS.
   The covariate experiment showed this is genetically identical to logistic/linear-with-covariates
   and ~30× faster (§4, §5).
3. **Covariates are `age_c`, `sex_c`, `age_c:sex_c`, and `PC1–PC10`** (dropped `age_c^2`, its
   interaction, and PC11–PC20 from the draft).
4. **No variant-level filters in PLINK2.** No `--mac`, no `--geno`, no `--hwe`. The bfile's variants
   already passed a MAF floor, a missingness ceiling, and a frequency-concordance check when it was
   built (§2). The `--threads` flag is dropped (the AoU box caps at 4 cores; PLINK2 autodetects).
5. **Minimum sample sizes:** binary phenotypes require **≥ 200 cases AND ≥ 200 controls**; ordinal /
   numeric phenotypes require **N ≥ 500** and ≥ 3 observed levels.
6. **Comprehensive ordinal coverage.** Every survey question in the v9 codebooks that can *arguably*
   be put on an ordinal scale is mapped. The mappings are machine-readable data files (§6, §7).
   Sensitive and stretched/uncertain questions are flagged, never silently dropped (§9).
7. **Downstream LDSC/rg (draft §21) is out of scope for this run.** We stop at summary statistics.
8. **PFHH conditions also get a genetic-relatedness-weighted family-burden sumscore** (self = 1,
   first-degree = 0.5, grandparent = 0.25, summed) as a quantitative liability-proxy phenotype
   alongside the binary `self_has_<condition>` (§11.2).
9. **EA/income use the exact repo conversions** (`setup_ea_gwas.py`/`setup_income_gwas.py`) and the
   same **age ≥ 26** restriction (§7.3).
10. **Cognitive + EA-proxy scores are GWASed** as external quantitative traits — the 4 ETM task
   scores plus teacher_z, SES-EA proxy, fine-tuned proxies, direct-XGBoost proxy, and the final
   g-EA proxy (§11b), residualized on sex + PCs.

---

## 1. Interpretation

This is a Pan-UKB-*style* phenotype-wide GWAS factory, not a SAIGE replication. It follows Pan-UKB
design principles (one GWAS per phenotype; all eligible survey-derived phenotypes; PHESANT-like
binary/ordinal/continuous handling; inverse-rank-normal transform for quantitative traits;
per-phenotype summary statistics) but the engine is PLINK2 fixed-effect linear regression in an
unrelated European sample. Because PLINK2 does not model relatedness, the unrelated-sample
restriction is mandatory.

---

## 2. Genotypes: the pre-built HapMap3 HQ bfile

Already built in the AoU workspace bucket by `make_hapmap3_bfile_hq.sh`:

```text
gs://<workspace-bucket>/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.{bed,bim,fam}
```

Verified contents:

```text
variants:  1,140,557   (autosomes 1–22, biallelic HapMap3 SNPs)
samples:     303,903   (all classified-European WGS participants)
requested HapMap3 rsids:      1,154,522
present in WGS metrics:       1,154,381
final HQ SNPs:                1,140,557
```

Variant filters **already applied** when the bfile was built (so **no** PLINK2 variant filters here):

```text
MAF                 >= 0.007   computed in the unrelated-EUR PCA-fit sample (fit_pca_iids)
variant missingness <= 0.01    computed in the classified-European sample
freq concordance    |fit-PCA ALT freq  -  SBayesRC snp.info ALT freq| <= 0.03
```

No MAC threshold, no separate `--geno`, no HWE filter is needed or applied.

---

## 3. Sample set (GWAS keep-list)

The GWAS analysis sample is the **unrelated European** subset, reusing the main pipeline's outputs:

```text
classified European ancestry            sbayesrc_genotypes/europeans/classified_european_iids.txt
AND confident genetic sex               sbayesrc_genotypes/genetic_sex/sex_covar.txt   (sex_01 in {0,1})
AND not in identical-component excl.    sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
AND unrelated at KING < 0.0441941       sbayesrc_genotypes/pca_eur/fit_pca_iids.txt    (PCA-fit unrelated set)
```

`fit_pca_iids.txt` is the third-degree-unrelated European set used to fit PCA (≈ 252,774 samples). It
is the primary keep-list. PLINK2 `--glm` is a fixed-effect model and cannot absorb kinship, so
close relatives must be removed up front; PCs and covariates handle residual structure.

Per-phenotype, the analysis sample is further intersected with participants who have a codeable
response / valid measurement and complete covariates (§5).

---

## 4. Method: residualize-first, covariate-free PLINK2

The covariate experiment (`~/projects/ukgwas/covariate_experiment`) established that, for downstream
h²/rg, running PLINK2 with covariates and running a covariate-free PLINK2 on a phenotype that was
pre-residualized on the same covariates are **genetically interchangeable**:

- continuous traits: `rg(with-covar, resid-first) = 1.0000` across 7 traits; all Δh² ≈ 0 (p ≈ 0.98–1.0).
- binary traits: `rg(logistic, LPM) = 1.00` for common traits; only the rare large-effect trait
  (red hair, 4.5% cases) drifts to 0.983 — the regime where you would switch to Firth/SAIGE anyway.
- speed: residualize-first, covariate-free PLINK2 was **~30× faster** for continuous traits and up to
  **~200×** faster for binary traits (median ~158×), because covariate-free linear regression is cheap.

So every phenotype here is reduced to a single pre-residualized quantitative vector and run through
one covariate-free PLINK2 linear pass.

### 4.1 Continuous / ordinal / numeric phenotypes

```text
1. raw       = mapped ordinal value (§6) OR parsed numeric value (§8) OR measurement (§10)
2. rint_raw  = inverse-rank-normal-transform(raw)          # rankdata(avg) -> norm.ppf((r-0.5)/n)
3. resid     = residual( rint_raw ~ 1 + age_c + sex_c + age_c:sex_c + PC1..PC10 )   # OLS
4. PLINK2 linear GWAS of `resid`, no covariates.
```

### 4.2 Binary phenotypes (linear probability model)

```text
1. y01    = 1 for the "case" answer, 0 for the eligible "control" answer(s), else missing (§6.1)
2. resid  = residual( y01 ~ 1 + age_c + sex_c + age_c:sex_c + PC1..PC10 )           # OLS, LPM
3. PLINK2 linear GWAS of `resid`, no covariates.                                    # NOT logistic
```

Binary phenotypes are **not** inverse-normal-transformed (they are LPM residuals); continuous/ordinal
phenotypes **are** INT'd before residualizing. This matches the covariate-experiment arms exactly
(`binary_resid_nocovar`, `rint_raw_resid_nocovar`).

### 4.3 IRNT and residualization definitions

```text
IRNT(x):  ranks = scipy.stats.rankdata(x, method="average")
          z     = scipy.stats.norm.ppf((ranks - 0.5) / len(x))
residual: beta  = lstsq([1, covars], y);  resid = y - [1, covars] @ beta
```

Both the raw and the residualized phenotype vectors are written to disk for auditability.

---

## 5. Covariates

```text
age_c              = age_at_event - mean(age_at_event within this phenotype's analysis sample)
sex_c              = sex_01 - 0.5                         # confirmed genetic sex, 0/1
age_c_sex_c_inter  = age_c * sex_c
PC1_AVG ... PC10_AVG                                       # from pca_eur/aou_projected.sscore
```

`age_at_event` is age at the selected survey response date (survey phenotypes) or the measurement
date (physical measurements). Covariates are centered within each phenotype's own complete-case
analysis sample. A phenotype's analysis sample = GWAS keep-list ∩ codeable response ∩ non-missing
`age`, `sex_01`, `PC1..PC10`.

---

## 6. Survey phenotype construction

Every codeable closed-ended survey item yields **binary one-vs-rest** phenotypes; every arguably
ordinal single-select item *additionally* yields an **ordinal linear** phenotype.

### 6.1 Binary response GWAS

```text
single-select question:  one binary phenotype per valid answer option
    case    = participant chose this option
    control = participant answered the question and chose another valid option
checkbox question:       one binary phenotype per option
    case    = option selected
    control = question answered (checkbox shown) and this option not selected
missing (never control): Skip, Prefer Not To Answer, Don't Know, Not sure, branch-not-asked,
                         invalid concept 2000000010, free-text-only, suppressed, conflicting dup
```

Exact complements (Yes vs No) are both kept for a complete atlas and flagged `exact_complement` in
the phenotype manifest.

### 6.2 Ordinal linear GWAS

For every single-select question with a defensible ordered scale (§7): map answer → numeric value,
set non-responses to missing, then run the §4.1 continuous pipeline. Both `<pheno>_raw` and the
residualized vector are written. Ordinal phenotypes with fewer than 3 observed levels in the analysis
sample fall back to binary-only.

### 6.3 Numeric survey GWAS

Free-numeric survey entries (counts, ages, durations, minutes): parse `value_as_number`, range-check
against the codebook min/max, drop impossible values, then run the §4.1 continuous pipeline (§8).

---

## 7. Ordinal mappings (machine-readable)

All ordinal mappings live in three generated data files under `metadata/` so the pipeline needs no
hand-editing on the AoU box. They are **answer-text driven**: an OMOP answer concept's name tail
matches the codebook answer label, so mapping on the normalized label works without a REDCap↔OMOP
crosswalk (this is how `setup_ses_ea_proxy_gwas.py::ordinal_value_from_answer` already works).

```text
metadata/ordinal_answer_templates.tsv   the shared answer-phrase rule library (rule -> label -> value)
metadata/ordinal_mapping_manifest.tsv   per (question, answer) -> ordinal value, one row each
metadata/survey_question_manifest.tsv   every question: its disposition, rule, sensitivity, counts
```

These are regenerated from the v9 codebook by `scripts/parse_codebooks.py` + `scripts/build_manifests.py`
(the rule knowledge base is `scripts/ordinal_rules.py`).

### 7.1 Mapping philosophy

```text
- Use the codebook scoring sheet / conventional scale direction.
- Higher value = more symptoms (GAD/PHQ/PSS/ACE), stronger agreement, more frequent behavior,
  better self-rated health/quality, or more of the named construct.
- Signed change scales: negative = less/decrease/loss, 0 = same/no change, positive = more/increase.
- Count/duration bands -> midpoints; open top bins -> lower bound or a documented top code.
- Genuinely nominal categories are NOT forced onto a scale -> one-vs-rest binary only.
- IRNT is applied after mapping and after sample/covariate filtering (§4.1).
- Education years and household income use the ea_proxy.md anchors.
```

### 7.2 Ordinal rule library (from `ea_proxy.md`, made explicit)

These are the shared answer-phrase templates. Each row family is a normalized `answer label → value`
table in `metadata/ordinal_answer_templates.tsv`; the count is how many v9 questions use it.

| Rule (template) | Scale | Direction | Example item |
| --- | --- | --- | --- |
| `phq_gad_0_3` | Not at all → Nearly every day | higher = more days | PHQ/GAD items |
| `freq_never_veryoften_0_4` | Never → Very often | higher = more frequent | ASRS/BFI frequency |
| `freq_pss_0_4` | Never → Very often (PSS) | higher = more frequent | PSS items |
| `freq_never_always5_0_4` | Never → Always | higher = more frequent | emotional-problem 7 days |
| `freq_never_always_0_4` | Never → Always (w/ Most of the time) | higher = more frequent | SDOH frequency |
| `freq_never_often_0_3` | Never → Often | higher = more frequent | loneliness items |
| `freq_event_0_5` | Never → Almost every day | higher = more frequent | discrimination |
| `freq_event_0_3` | Never → Almost every day (4-lvl) | higher = more frequent | discrimination (short) |
| `freq_covid_contact_0_3` | Only a few times → Daily | higher = more frequent | COPE contact |
| `time_none_all_0_4` | None → All of the time | higher = more of the time | MOS social support |
| `time_all_none_0_4` | All → None of the time | higher = more of the time | SF-style items |
| `days_last5_midpoint` | 0 / 1-2 / 3-4 / every day | days (midpoints) | MSDS |
| `agree_bfi_1_5` | Disagree strongly → Agree strongly | higher = more agreement | BFI-2-XS |
| `agree_1_4` | Strongly disagree → Strongly agree (no mid) | higher = more agreement | SDOH agreement |
| `agree_1_5` / `agree_neutral_1_5` | Strongly disagree → Strongly agree | higher = more agreement | SCNS/SDOH |
| `agree_somewhat_1_4` | Strongly disagree → Strongly agree (somewhat) | higher = more agreement | SDOH neighborhood |
| `agree_lotr_1_5` | LOT-R optimism | higher = more agreement | LOT-R |
| `intensity_0_4` | Not at all → Extremely | higher = more intense | IES-R / COVID impact |
| `amount_0_4` | Not at all → An extreme amount | higher = more | COPE impact |
| `distress_0_4` | Positive → Very distressing | higher = more distress | psychosis-experience distress |
| `health_1_5` | Poor → Excellent | higher = better | PROMIS global health |
| `confidence_1_5` | Not at all → Extremely confident | higher = more confident | medical-form confidence |
| `ability_extent_0_4` | Not at all → Completely | higher = more able | everyday activities |
| `severity_0_4` | None → Very severe | higher = worse | fatigue |
| `difficulty_0_4` | Unable to do → Without any difficulty | higher = more able | PROMIS |
| `audit_freq_0_4` | AUDIT-C drinking frequency | higher = more frequent | alcohol |
| `binge_freq_0_4` | 6+ drinks frequency | higher = more frequent | alcohol |
| `subuse_lifestyle_0_4` | Never → Daily or almost daily | higher = more frequent | substance use |
| `smoke_freq_0_2` | Not at all / Some days / Every day | higher = more frequent | smoking |
| `current_use_0_3` | No never → Yes every day | higher = more current use | smoking status |
| `change_lessmore_signed` | Less / Same / More | signed | COPE behavior change |
| `change_muchmore_signed` | Much less → Much more | signed | COPE change |
| `recency_lifetime_0_2` | Never / not-recent / recent | higher = more recent | partner-violence recency |
| `spiritual_frequency_0_5` | BMMRS Never → Many times a day | higher = more frequent | BMMRS |
| `happiness_0_5` | Extremely unhappy → Extremely happy | higher = happier | well-being |
| `count_band_midpoint`, `drink_count_band_midpoint`, `visit_count_band_midpoint` | count bands | midpoints | panic attacks, drinks, visits |

Plus per-question overrides for scales no template captures — education years, income $k, English
proficiency, religious-service frequency, symptom-onset month, likelihood, worry, importance,
housing density, ovary-removal count, and the mhqukb duration/appetite/weight/experience items — all
enumerated with explicit per-answer values in `metadata/ordinal_mapping_manifest.tsv`.

### 7.3 Education and income anchors (ea_proxy.md)

These reproduce the exact `EA_MAPPING` / `INCOME_MAPPING` in the repo's
`setup_ea_gwas.py` and `setup_income_gwas.py` (values verified identical).

```text
Education (question "highest grade", -> years):
  never attended / grades 1-4 / grades 5-8 = 9;  grades 9-11 = 10;  grade 12 or GED = 13;
  1-3 yrs after HS = 15;  college graduate = 18;  advanced degree = 20
Income (annual household, -> $k midpoints, top-coded 250):
  <10k=5, 10-25k=17.5, 25-35k=30, 35-50k=42.5, 50-75k=62.5, 75-100k=87.5,
  100-150k=125, 150-200k=175, >200k=250
```

Also matching those scripts: **every education and income phenotype (ordinal and
each binary one-vs-rest) restricts to respondents aged ≥ 26 at the survey response**
(`--min-age-at-survey 26`), so people who may not have finished education / are
early-career are excluded.

---

## 8. Numeric survey phenotypes

53 free-numeric survey entries across the included surveys (household size, cigarettes/day, years
smoked, ages at first/last episode, episode/attempt counts, IPAQ minutes/hours/days, COVID test
counts, etc.) are parsed as numbers, range-checked against the codebook min/max, impossible values
dropped, then INT'd and residualized (§4.1). Companion "too many to count / one episode ran into the
next" categorical top-codes are treated as missing for the numeric phenotype (they remain valid
one-vs-rest binary answers). SSN, phone, address, and other PII numeric fields are excluded.

---

## 9. Flagged questions (for your consideration — nothing dropped silently)

`metadata/flagged_questions.tsv` lists every question that needs human review before release. Reasons:

```text
sensitive                    272   sensitive topic (still analyzed; flag = release-review tier)
medium_confidence_ordinal     39   ordinal mapping is defensible but a stretch (bands, signed,
                                    religiosity, recency, housing density, symptom-onset month, ...)
uncertain_ordinal              1   >=3 options, no confident scale -> binary only + review
                                    (`cu_covid`, "what type of household", partially ordered by size)
```

Sensitive topics tagged (analyzed, `sensitive_release=true`): sexual orientation, gender identity,
suicidality/self-harm, sexual behavior/trauma, trauma/violence, substance use, alcohol,
reproductive/menstrual, mental health, COVID-related, financial hardship, justice involvement,
immigration/origin, religion, disability.

### 9.1 Deliberate exclusions to reconsider

The draft spec (and this run) exclude these blocks by design. They are **available** if you want them
— call them out explicitly rather than assume:

```text
Personal Medical History  (453 single-select, self medical-history follow-ups: still-seeing-doctor,
                           age-at-diagnosis band, currently-on-meds per condition) -> excluded as
                           draft §11.1 do-not-run items. The "age when first told" life-stage bands
                           (Child 0-11 / Adolescent 12-17 / Adult 18-64 / Older adult 65-74) ARE
                           ordinal if you decide to include them.
Family Health History     (relative-specific "who in your family has X") -> excluded (not a clean
                           participant phenotype).
Personal & Family Health  broad family-history items -> excluded EXCEPT the 33 self-history
  History                 phenotypes in metadata/pfhh_self_allowlist.tsv (neuro/mental/substance +
                           fibromyalgia + recent fracture), run as binary self_has_<condition> (§11).
race / ethnicity / PII    excluded (12 items): not GWAS phenotypes in an EUR-stratified analysis.
free-text / date-only     excluded (237 free-text, 33 date-only): no derived phenotype here.
```

If you want any of these included, they are already parsed in `metadata/survey_item_inventory.tsv`;
flip the disposition rule in `scripts/build_manifests.py`.

---

## 10. Physical measurements

Extracted from OMOP `measurement`, unit-converted, impossible values dropped, same-date multiples
averaged after outlier removal, closest-to-baseline record used for the primary phenotype. All run
through the §4.1 continuous pipeline.

```text
height_cm                          100–250
bmi_kg_m2                          12–80    exclude records where pregnant at measurement date
systolic_bp_mmhg                   70–260
diastolic_bp_mmhg                  40–160   require SBP & DBP from same visit; drop SBP<=DBP
pulse_pressure_mmhg  = SBP - DBP   15–150
mean_arterial_pressure = DBP + (SBP-DBP)/3   50–180
heart_rate_or_pulse_bpm            30–220
```

BMI pregnancy exclusion is per-record (exclude only the affected BMI measurement, not the
participant), using a visit-level pregnancy flag or same-day pregnancy-status observation.

---

## 11. PFHH self-history allowlist

`metadata/pfhh_self_allowlist.tsv` — 33 allowlisted conditions (neuro/nervous-system,
mental-health/substance-use, fibromyalgia, recent fracture). Each condition yields **two**
phenotypes.

### 11.1 `pfhh_self_has_<condition>` — binary (self-only)

```text
case    = participant selected "Self" on "Including yourself, who in your family has had <condition>?"
control = completed the relevant PFHH category screen and did not self-report that condition
missing = did not complete the section / refused / denominator unrecoverable
```

### 11.2 `pfhh_burden_<condition>` — quantitative genetic-relatedness-weighted family burden

An aggregate genetic-liability proxy: the sum, over the relations the participant selected for that
condition, of each relation's coefficient of relationship to the participant.

```text
weight(Self)                                   = 1.00
weight(first-degree: parent/sibling/son/daughter) = 0.50   each
weight(second-degree: grandparent, half-sibling)  = 0.25   each
score  = sum of selected relations' weights          e.g. Self + Mother + Grandparent = 1.75
       = 0 for participants who completed the category screen but reported no one (or were not
         shown the condition question because no family member has it)
missing = answered only PMI (Skip/Prefer not/Don't know), or did not complete the screen
```

The score is then INT'd and residualized like any quantitative trait (§4.1). It is a within-family
liability aggregate for the proband, **not** a phenotype about any individual relative — no
mother/father/sibling GWAS is run on its own. Applied to all 33 allowlisted conditions (the same set
as §11.1); trivially extends to the full PFHH condition list if the diagnosis restriction is lifted.
The `pfhh_burden_*` phenotypes are heavily zero-inflated; IRNT handles the ties, and the ≥ 3 observed
levels rule (§12) drops any condition too rare to score.

---

## 11b. External cognitive / EA-proxy score GWAS

The ea_proxy workflow produces pre-computed continuous scores that are GWASed here
as external quantitative traits. Because they are already age/sex-normalized
upstream, they are INT'd and residualized on **sex_c + PC1..PC10 only** (no age
term) — matching the repo's final g-EA proxy GWAS covariates. The registry is
`metadata/external_scores.tsv` (score file paths + column names; missing files are
skipped with a warning). Phenotype ids are prefixed `cog_`.

ETM cognitive task scores (recommended per-task summary score, `*_z_age_sex`):

```text
dd_patience_z_age_sex     Delay Discounting patience = -lnk (negative log mean discount rate k);
                          higher = more patient / less delay discounting.
gradcpt_perf_z_age_sex    GradCPT sustained attention = PC1 of d-prime, -log(RT CV), -log(median RT);
                          higher = better sustained attention.
flanker_efficiency_z_age_sex  Flanker inhibitory control = official AoU Flanker rate-correct score (0-100);
                          higher = better attentional control.
emorecog_perf_z_age_sex   Emotional Recognition = PC1 of accuracy, -log(RT CV), -log(median RT);
                          higher = better emotion recognition.
```

EA / SES / cognitive proxy scores:

```text
teacher_z                  EA years residualized on yob, sex, yob:sex, z-scored (the EA teacher label).
ses_ea_proxy_z             cross-fit XGBoost survey/area-SES -> EA-years proxy (no genotypes).
gradcpt_flanker_finetuned_ea_proxy_z   SES-EA boosters fine-tuned toward the GradCPT+Flanker mean.
gradcpt_flanker_direct_xgb_proxy_z     scratch XGBoost survey -> GradCPT+Flanker proxy (not fine-tuned).
g4_finetuned_ea_proxy_z    SES-EA boosters fine-tuned toward the 4-domain ETM-g factor.
gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z   the final selected g-EA proxy (headline).
```

## 12. Phenotype-eligibility thresholds

```text
Binary:            cases >= 200  AND  controls >= 200
Ordinal/numeric:   N >= 500  AND  >= 3 observed levels  AND  no single level > 98%  AND finite variance
```

(Recommended-for-LDSC tiers, e.g. h² z-scores, are a downstream concern and not gated here.)

---

## 13. PLINK2 command (the only association call)

One covariate-free linear pass per phenotype over the whole autosomal bfile:

```bash
plink2 \
  --bfile ${HM3_BFILE} \
  --keep sample_qc/unrelated_eur.keep \
  --pheno phenotypes/${PHENO}.resid.pheno.tsv \
  --pheno-name ${PHENO}_resid \
  --glm allow-no-covars \
  --no-input-missing-phenotype \
  --out gwas/${PHENO}/${PHENO}
```

No `--mac`, `--geno`, `--hwe`, or `--threads` (see §0, §2). Output: `${PHENO}.${PHENO}_resid.glm.linear`.
The keep-list is per-phenotype (analysis complete-cases); the phenotype file already holds only the
residualized values for those samples.

---

## 14. Candidate GWAS counts (from the v9 codebook, this classification)

Pre-QC candidate counts; actual runnable counts are lower after the §12 N/case filters.

| Category | Questions | Binary | Ordinal | Numeric | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| The Basics | 29 | 268 | 2 | 2 | 272 |
| Life Functioning (+Basics disability) | 6 | 12 | 0 | 0 | 12 |
| Lifestyle | 36 | 105 | 18 | 5 | 128 |
| Overall Health | 24 | 95 | 13 | 3 | 111 |
| Healthcare Access & Utilization | 57 | 202 | 20 | 0 | 222 |
| Social Determinants of Health | 81 | 395 | 77 | 1 | 473 |
| COVID-19 Participant Experience (COPE) | 222 | 914 | 129 | 25 | 1068 |
| Minute Survey on COVID-19 Vaccines | 72 | 499 | 3 | 0 | 502 |
| Emotional Health History & Well-Being | 103 | 311 | 58 | 8 | 377 |
| Behavioral Health & Personality | 60 | 176 | 25 | 9 | 210 |
| PFHH self-history allowlist (binary self_has) | 33 | 33 | 0 | 0 | 33 |
| PFHH relatedness-burden sumscore (quant) | 33 | 0 | 33 | 0 | 33 |
| Physical measurements | 9 | 0 | 0 | 9 | 9 |
| Cognitive / EA-proxy external scores | 10 | 0 | 0 | 10 | 10 |
| **TOTAL** | **782** | **3010** | **378** | **72** | **3460** |

(The 33 `pfhh_burden_*` sumscores and the 10 `cog_*` external scores are counted as
quantitative traits. Physical measurements now include pulse pressure and MAP.)

---

## 15. Output files

```text
metadata/survey_item_inventory.tsv        every codebook item, field type, options (audit source)
metadata/survey_question_manifest.tsv     every question: disposition, ordinal rule, sensitivity
metadata/ordinal_answer_templates.tsv     shared ordinal rule library (rule -> label -> value)
metadata/ordinal_mapping_manifest.tsv     per (ordinal question, answer) -> value
metadata/flagged_questions.tsv            sensitive + uncertain/stretched items for review
metadata/pfhh_self_allowlist.tsv          33 self_has_<condition> phenotypes
phenotypes/<pheno>.raw.pheno.tsv          raw + residualized phenotype vectors (audit)
phenotypes/<pheno>.resid.pheno.tsv        the vector PLINK2 reads
gwas/<pheno>/<pheno>.<pheno>_resid.glm.linear      per-phenotype summary statistics
gwas/<pheno>/<pheno>.sumstats.tsv.gz      lightweight columns (rsid,a1,a1freq,beta,se,p,n,...)
metadata/phenotype_manifest.tsv           one row per run phenotype: N, cases/controls, rule, paths
metadata/run_manifest.json               inputs, covariates, sample counts, timing
```

Lightweight sumstats columns: `chrom pos rsid a1 a2 a1freq beta se t_or_z p n phenotype trait_type
case_count control_count`.

---

## 16. How to run on the AoU Researcher Workbench

```bash
# From inside the AoU Verily Jupyter terminal, after the main pipeline has produced
# the European keep-list, genetic sex, PCs, and the HapMap3 HQ bfile.
cd aou-sbayesrc-gwas/pan_aou_gwas

# (once, off-platform or on) regenerate the phenotype manifests from the codebook:
python3 scripts/parse_codebooks.py
python3 scripts/build_manifests.py

# extract survey/measurement data, build residualized phenotypes, run PLINK2:
bash run_pan_aou_gwas.sh --setup-only     # extract + build phenotypes + manifests, no GWAS
bash run_pan_aou_gwas.sh --smoke          # a few phenotypes end-to-end as a smoke test
bash run_pan_aou_gwas.sh                   # full run
```
