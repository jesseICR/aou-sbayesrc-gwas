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
11. **Validated composite scores are built** from the codebook Scoring sheets — GAD-7, PHQ-9, PSS,
   ACE, IES, ASRS, AUDIT-C, Everyday Discrimination, loneliness, social support, PROMIS, resilience,
   spiritual experience, health literacy, the BFI-2 Big Five domains, and the neighborhood social
   cohesion / disorder / walkability and Hunger Vital Sign composites (§11c), as prorated sums with
   per-scale reverse-keying (opposite-valence items flipped).
12. **Wearable (Fitbit) phenotypes** — mean daily steps, sedentary/active minutes, sleep duration,
   efficiency, and a chronotype (sleep-onset) proxy (§10b), on the Fitbit subcohort.
13. **Derived psychiatric phenotypes** (§11d) from the UKB-MHQ/CIDI/PCL items — probable MDD (incl.
   recurrent), probable bipolar/mania, psychotic experiences, PCL-PTSD, and lifetime suicidal
   ideation / attempt (all sensitive).
14. **Acculturation index** (§11e) — US-born + English-at-home + English proficiency.
15. **Geographic / political state-cluster membership** (§11f) — 12 binary Census-region / subregion /
   political-grouping phenotypes (primarily capture residual geographic structure).

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
AND pan-AoU binary sex covariate         data/pan_aou_gwas_work/sample_qc/pan_aou_sex_covar.txt
AND not in identical-component excl.    sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
AND unrelated at KING < 0.0441941       sbayesrc_genotypes/pca_eur/fit_pca_iids.txt    (PCA-fit unrelated set)
```

`pan_aou_sex_covar.txt` is built by `scripts/build_pan_aou_sex_covar.py`. It starts from the main
pipeline's strict `genetic_sex/sex_covar.txt` and then adds the following pre-specified rows from
`genetic_sex/sex_ploidy_qc.tsv`: assigned sex at birth Male with DRAGEN `X0`/`XO` is coded male,
and skipped/prefer-not-to-answer sex-at-birth rows with DRAGEN `XX`/`XY` are coded from the DRAGEN
binary ploidy. All other nonbinary, missing, discordant, or noncanonical sex/ploidy combinations
remain excluded from the pan-AoU GWAS sample.

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
sex_c              = sex_01 - 0.5                         # pan-AoU binary sex covariate, 0/1
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

### 9.2 SES-EA XGBoost feature-source supplement

The primary survey metadata is codebook-derived, then linked back to live AoU
question IDs at runtime. The SES-EA and direct GradCPT/Flanker XGBoost models
were trained from live v9 question IDs, and some of those IDs do not round-trip
through the codebook text matcher. `metadata/ea_proxy_feature_sources.tsv`
therefore acts as a supplemental include-list for source questions used by those
models. It only fills question IDs missing from the normal manifest; codebook
metadata remains authoritative wherever both sources exist.

The supplemental source covers the same survey questions used by the XGBoost
feature contract. Ordinal supplemental rows use the same answer-text parser as
the SES-EA setup code; one-hot source questions are emitted as one-vs-rest binary
phenotypes; numeric source questions are emitted as continuous phenotypes when
the response text parses as a number. Technical XGBoost columns such as
survey-taken flags, age-at-survey columns, genetic sex, and curated nonresponse
indicators are not treated as participant traits unless explicitly added as such.

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

## 10b. Wearable (Fitbit) phenotypes

Per-person averages over the participant's Fitbit wear, requiring ≥ 10 valid days, run through the
§4.1 continuous pipeline (residualized on the full covariate set; age at wear is known):

```text
fitbit_mean_daily_steps     mean daily step count (days with steps > 0)
fitbit_sedentary_minutes    mean daily sedentary minutes
fitbit_active_minutes       mean daily fairly+very active minutes
fitbit_sleep_minutes        mean nightly minutes asleep (main sleep)
fitbit_sleep_efficiency     mean minute_asleep / minute_in_bed
fitbit_chronotype_sleep_onset  mean main-sleep onset clock hour (chronotype proxy;
                            onset before noon wrapped to [24,36), higher = later/evening)
```

Chronotype uses the `sleep_level` start times (earliest main-sleep segment per night); confirm that
table's schema on-platform.

Sourced from the AoU `activity_summary` and `sleep_daily_summary` Fitbit tables (a smaller wearable
subcohort). The orchestrator extracts them and skips gracefully if the tables are absent
(`PAN_AOU_SKIP_FITBIT=1` to disable). Confirm the table/column names on-platform.

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

## 11b.1 ZIP3 socioeconomic context GWAS

AoU `ds_zip_code_socioeconomic` provides the latest ZIP3-level socioeconomic
context row per participant. The pipeline treats the seven numeric fields as
quantitative contextual phenotypes, then applies the standard continuous-trait
pipeline (§4.1): inverse-rank-normal transform followed by residualization on
age at observation, sex, age×sex, and PC1..PC10. Raw ZIP3 and ACS vintage are
kept in the local extract for auditability but are not GWASed.

Phenotype ids:

```text
zip3_deprivation_index
zip3_median_income
zip3_fraction_poverty
zip3_fraction_assisted_income
zip3_fraction_no_health_ins
zip3_fraction_vacant_housing
zip3_fraction_high_school_edu
```

## 11c. Validated composite score definitions

Each composite is a **prorated sum**: mean(available item scores) × n_items, requiring ≥ 80% of items answered. Reverse-worded items (flagged per scale) are flipped on their own min/max before summing. Items are matched to survey responses by question text and merged across survey administrations. The score is then inverse-normal-transformed and residualized like any quantitative trait (§4.1). Phenotype ids are prefixed `comp_`.

### GAD-7 — Generalized Anxiety Disorder scale (anxiety)

- **Items:** 7
- **Per-item scoring:** Not at all = 0, Several days = 1, Over half the days = 2, Nearly all days = 3
- **Total score:** prorated sum of 7 items; no reverse-keyed items
- **Auto-built:** yes (comp_gad7_anxiety)
- **Questions:**
    - Feeling nervous, anxious, or on edge
    - Not being able to stop or control worrying
    - Worrying too much about different things
    - Trouble relaxing
    - Being so restless that it's hard to sit still
    - Becoming easily annoyed or irritable
    - Feeling afraid as if something awful might happen

### PHQ-9 — Patient Health Questionnaire (depression)

- **Items:** 9
- **Per-item scoring:** Not at all = 0, Several days = 1, Over half the days = 2, Nearly all days = 3
- **Total score:** prorated sum of 9 items; no reverse-keyed items
- **Auto-built:** yes (comp_phq9_depression)
- **Questions:**
    - Little interest or pleasure in doing things
    - Feeling down, depressed, or hopeless
    - Trouble falling or staying asleep, or sleeping too much
    - Feeling tired or having little energy
    - Poor appetite or overeating
    - Feeling bad about yourself or that you are a failure or have let yourself or your family down
    - Trouble concentrating on things, such as reading the newspaper or watching television
    - Moving or speaking so slowly that other people could have noticed? Or the opposite - being so fidgety or restless that you have been moving around a lot more than usual
    - Thoughts that you would be better off dead or of hurting yourself in some way

### PSS — Perceived Stress Scale

- **Items:** 10
- **Per-item scoring:** 2 answer scales across items (shown per item below)
- **Total score:** prorated sum of 10 items; 4 reverse-keyed
- **Auto-built:** yes (comp_pss_perceived_stress)
- **Questions:**
    - In the last month, how often have you been upset because of something that happened unexpectedly?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt that you were unable to control the important things in your life?  — [Never=0.0, Almost never=1.0, Sometime=2.0, Fairly often=3.0, Often=4.0]
    - In the last month, how often have you felt nervous and "stressed?"  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt confident about your ability to handle your personal problems? *(reverse-keyed)*  — [Never=0.0, Almost never=1.0, Sometime=2.0, Fairly often=3.0, Often=4.0]
    - In the last month, how often have you felt that things were going your way? *(reverse-keyed)*  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you found that you could not cope with all the things that you had to do?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you been able to control irritations in your life? *(reverse-keyed)*  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt that you were on top of things? *(reverse-keyed)*  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you been angered because of things that were outside of your control?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt difficulties were piling up so high that you could not overcome them?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]

### ACE — Adverse Childhood Experiences

- **Items:** 11
- **Per-item scoring:** 3 answer scales across items (shown per item below)
- **Total score:** prorated sum of 11 items; no reverse-keyed items
- **Auto-built:** yes (comp_ace_adversity)
- **Questions:**
    - During your first 18 years of life, did you live with anyone who was depressed, mentally ill, or suicidal? (ACE category: Mentally ill household member)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, did you live with anyone who was a problem drinker or alcoholic? (ACE category: Substance abuse in household)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, did you live with anyone who used illegal street drugs or who abused prescription medications? (ACE category: Substance abuse in household)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, did you live with anyone who served time or was sentenced to serve time in a prison, jail, or other correctional facility? (ACE category: Incarcerated household member)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, were your parents separated or divorced? (ACE category: Parental separation/divorce)  — [Yes=1.0, No=0.0, Parents not married=0.0]
    - During your first 18 years of life, how often did your parents or adults in your home ever slap, hit, kick, punch or beat each other up? (ACE category: Violence between adults in household)  — [Never=0.0, Once=1.0, More than once=1.0]
    - Before age 18, how often did a parent or adult in your home ever hit, beat, kick, or physically hurt you in any way? Do not include spanking. (ACE category: Physical abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did a parent or adult in your home ever swear at you, insult you, or put you down? (ACE category: Emotional abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did anyone at least 5 years older than you or an adult, ever touch you sexually? (ACE category: Sexual abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did anyone at least 5 years older than you or an adult, try to make you touch them sexually? (ACE category: Sexual abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did anyone at least 5 years older than you or an adult, force you to have sex? (ACE category: Sexual abuse)  — [Never=0.0, Once=1.0, More than once=1.0]

### IES — Impact of Event Scale (event-related distress)

- **Items:** 6
- **Per-item scoring:** Not at all = 0, A little bit = 1, Moderately = 2, Quite a bit = 3, Extremely = 4
- **Total score:** prorated sum of 6 items; no reverse-keyed items
- **Auto-built:** yes (comp_ies_event_impact)
- **Questions:**
    - In the past 7 days, I thought about COVID-19 when I didn't mean to.
    - In the past 7 days, I felt watchful or on-guard.
    - In the past 7 days, other things kept making me think about COVID-19.
    - In the past 7 days, I was aware that I still had a lot of feelings about COVID-19, but I didn't deal with them.
    - In the past 7 days, I tried not to think about COVID-19.
    - In the past 7 days, I had trouble concentrating.

### ASRS — Adult ADHD Self-Report Scale (Part A screener)

- **Items:** 6
- **Per-item scoring:** 2 answer scales across items (shown per item below)
- **Total score:** prorated sum of 6 items; no reverse-keyed items
- **Auto-built:** yes (comp_asrs_adhd)
- **Questions:**
    - How often do you have trouble wrapping up the final details of a project, once the challenging parts have been done?  — [Never=0.0, Rarely=0.0, Sometimes=1.0, Often=1.0, Very often=1.0]
    - How often do you have difficulty getting things in order when you have to do a task that requires organization?  — [Never=0.0, Rarely=0.0, Sometimes=1.0, Often=1.0, Very often=1.0]
    - How often do you have problems remembering appointments or obligations?  — [Never=0.0, Rarely=0.0, Sometimes=1.0, Often=1.0, Very often=1.0]
    - When you have a task that requires a lot of thought, how often do you avoid or delay getting started?  — [Never=0.0, Rarely=0.0, Sometimes=0.0, Often=1.0, Very often=1.0]
    - How often do you fidget or squirm with your hands or feet when you have to sit down for a long time?  — [Never=0.0, Rarely=0.0, Sometimes=0.0, Often=1.0, Very often=1.0]
    - How often do you feel overly active and compelled to do things, like you were driven by a motor?  — [Never=0.0, Rarely=0.0, Sometimes=0.0, Often=1.0, Very often=1.0]

### UCLA / ULS-8 — Loneliness

- **Items:** 8
- **Per-item scoring:** 2 answer scales across items (shown per item below)
- **Total score:** prorated sum of 8 items; 2 reverse-keyed
- **Auto-built:** yes (comp_ucla_loneliness)
- **Questions:**
    - I lack companionship  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - There is no one I can turn to  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - I am an outgoing person *(reverse-keyed)*  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0]
    - I feel left out  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - I feel isolated from others  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - I can find companionship when I want it *(reverse-keyed)*  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0]
    - I am unhappy being so withdrawn  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - People are around me but not with me  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]

### Everyday Discrimination Scale

- **Items:** 9
- **Per-item scoring:** Almost everyday = 6, At least once a week = 5, A few times a month = 4, A few times a year = 3, Less than once a year = 2, Never = 1
- **Total score:** prorated sum of 9 items; no reverse-keyed items
- **Auto-built:** yes (comp_everyday_discrimination)
- **Questions:**
    - You are treated with less courtesy than other people are.
    - You are treated with less respect than other people are.
    - You receive poorer service than other people at restaurants or stores.
    - People act as if they are afraid of you.
    - People act as if they're better than you are.
    - You are called names or insulted.
    - You are threatened or harassed.
    - People act as if they think you are not smart.
    - People act as if they think you are dishonest.

### MOS Social Support (RAND) + Tangible subscale

- **Items:** 9
- **Per-item scoring:** None of the time = 1, A little of the time = 2, Some of the time = 3, Most of the time = 4, All of the time = 5
- **Total score:** prorated sum of 9 items; no reverse-keyed items
- **Auto-built:** yes (comp_social_support)
- **Questions:**
    - Someone to help you if you were confined to bed
    - Someone to take you to the doctor if you needed it
    - Someone to prepare your meals if you were unable to do it yourself
    - Someone to help with daily chores if you were sick
    - Someone to take you to the doctor if you need it
    - Someone to have a good time with
    - Someone to turn to for suggestions about how to deal with a personal problem
    - Someone who understands your problems
    - Someone to love and make you feel wanted

### Neighborhood, walkability & food-insecurity composites

Built directly from the survey items (reusing their ordinal scores), because the scoring sheet groups these with mixed item valence. Opposite-valence items are reverse-keyed.

#### comp_social_cohesion

- Neighborhood social cohesion; higher = more cohesion.
- **Items:** 4; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - People around here are willing to help their neighbors.
    - People in my neighborhood generally get along with each other.
    - People in my neighborhood can be trusted.
    - People in my neighborhood share the same values.

#### comp_neighborhood_disorder

- Perceived neighborhood disorder (order items reversed); higher = more disorder.
- **Items:** 13; **reverse-keyed:** 4; prorated sum
- **Questions:**
    - There is a lot of graffiti in my neighborhood.
    - My neighborhood is noisy.
    - Vandalism is common in my neighborhood.
    - There are lot of abandoned buildings in my neighborhood.
    - There are too many people hanging around on the streets near my home.
    - There is a lot of crime in my neighborhood.
    - There is too much drug use in my neighborhood.
    - There is too much alcohol use in my neighborhood.
    - I'm always having trouble with my neighbors.
    - My neighborhood is clean. *(reverse-keyed)*
    - People in my neighborhood take good care of their houses and apartments. *(reverse-keyed)*
    - In my neighborhood, people watch out for each other. *(reverse-keyed)*
    - My neighborhood is safe. *(reverse-keyed)*

#### comp_neighborhood_physical_disorder

- Physical disorder subscale (order items reversed); higher = more disorder.
- **Items:** 6; **reverse-keyed:** 2; prorated sum
- **Questions:**
    - There is a lot of graffiti in my neighborhood.
    - My neighborhood is noisy.
    - Vandalism is common in my neighborhood.
    - There are lot of abandoned buildings in my neighborhood.
    - My neighborhood is clean. *(reverse-keyed)*
    - People in my neighborhood take good care of their houses and apartments. *(reverse-keyed)*

#### comp_neighborhood_social_disorder

- Social disorder subscale (order items reversed); higher = more disorder.
- **Items:** 7; **reverse-keyed:** 2; prorated sum
- **Questions:**
    - There are too many people hanging around on the streets near my home.
    - There is a lot of crime in my neighborhood.
    - There is too much drug use in my neighborhood.
    - There is too much alcohol use in my neighborhood.
    - I'm always having trouble with my neighbors.
    - In my neighborhood, people watch out for each other. *(reverse-keyed)*
    - My neighborhood is safe. *(reverse-keyed)*

#### comp_neighborhood_walkability

- PANES neighborhood walkability (crime-safety items reversed); higher = more walkable.
- **Items:** 7; **reverse-keyed:** 2; prorated sum
- **Questions:**
    - Many shops, stores, markets or other places to buy things I need are within easy walking distance of my home. Would you say that you...
    - It is within a 10-15 minute walk to a transit stop (such as bus, train, trolley, or tram) from my home. Would you say that you...
    - There are sidewalks on most of the streets in my neighborhood. Would you say that you...
    - There are facilities to bicycle in or near my neighborhood, such as special lanes, separate paths or trails, or shared use paths for cycles and pedestrians. Would you say that you...
    - My neighborhood has several free or low-cost recreation facilities, such as parks, walking trails, bike paths, recreation centers, playgrounds, public swimming pools, etc. Would you say that you...
    - The crime rate in my neighborhood makes it unsafe to go on walks at night. Would you say that you... *(reverse-keyed)*
    - The crime rate in my neighborhood makes it unsafe to go on walks during the day. Would you say that you... *(reverse-keyed)*

#### comp_hunger_vital_sign

- Hunger Vital Sign food-insecurity screener; higher = more food insecurity.
- **Items:** 2; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - Within the past 12 months, we worried whether our food would run out before we got money to buy more.
    - Within the past 12 months, the food we bought just didn't last and we didn't have money to get more.

#### comp_ptsd_pcl

- PTSD symptoms (abbreviated PCL, 5 items, 0-4 each); higher = more symptoms.
- **Items:** 5; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - In the past month, have you had repeated, disturbing memories, thoughts, or images of a stressful experience from the past?
    - In the past month, have you felt very upset when something reminded you of a stressful experience from the past?
    - In the past month, have you avoided activities or situations because they reminded you of a stressful experience from the past?
    - In the past month, have you felt distant or cut off from other people?
    - In the past month, have you felt irritable or had angry outbursts?

#### comp_subjective_wellbeing

- Subjective well-being (happiness + life meaning, UKB-style); higher = greater well-being.
- **Items:** 2; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - In general, how happy are you?
    - To what extent do you feel your life to be meaningful?

## 11d. Derived psychiatric phenotypes (UKB-MHQ / CIDI-SF / PCL)

AoU imported the UKB Mental Health Questionnaire, so the published algorithmic phenotypes are
derivable from the raw items. These are screening-level derivations (documented assumptions below),
**all sensitive** (mental health / suicidality → sensitive release tier). Binary unless noted.

```text
psych_psychotic_experiences_any     any Yes to voices / thought-insertion / paranoia (cidi5_21/22/23)
psych_self_harm_ideation_lifetime   ss_1  (ever thoughts of purposely hurting yourself)
psych_suicidal_ideation_lifetime    ss_2  (ever thoughts of killing yourself)
psych_suicide_attempt_lifetime      ss_3  (ever a suicide attempt)
psych_mania_episode_screen          (ever high/hyper OR irritable) AND >=3 manic symptoms (mhqukb_43/44/45)
psych_probable_bipolar              mania screen AND >=4-day duration (mhqukb_46) AND impairment (mhqukb_47)
psych_lifetime_depressed_episode    ever a >=2-week low-mood / anhedonia period (mhqukb_5/6)
psych_probable_recurrent_depression lifetime depressed episode AND several episodes (mhqukb_24)
ptsd_pcl (composite, §11c)          abbreviated PCL 5-item sum (pcl_1..5), continuous
```

The depression and bipolar derivations follow the UKB Smith et al. 2013 logic at the item level; they
are not the full CIDI symptom-count diagnoses, so they read as "probable"/"screen". Controls are
participants who completed the relevant module and do not meet criteria.

## 11e. Acculturation index

A cultural-assimilation score (`accult_index`, quantitative; higher = more acculturated):

```text
US-born (birthplace = USA)               + 1 / 0
English spoken at home (chis_1 = No)     + 1 / 0
English proficiency (chis_1_xx)          + 0..1  (imputed 1 for English-at-home speakers)
```

Summed to a 0–3 index, then INT'd and residualized (§4.1). Mainly informative across the immigrant /
language-minority gradient; US-born English-at-home participants sit at the ceiling.

## 11f. Geographic / political state-cluster membership

One binary GWAS per state cluster (member vs non-member), from `metadata/state_clusters.tsv` — 12
clusters: the 4 Census regions (Northeast, Midwest, South, West), 4 subregions (New England, Great
Lakes, Rocky Mountain, Sunbelt), and 4 political groupings (Swing, Solid Blue, Lean Blue, Solid Red).
Phenotype ids `geo_<cluster>`.

State source: the participant's **work-address state** (`employmentworkaddress_state`) by default —
the only participant-linked state in the survey, since home address is privacy-suppressed. For true
**residence state**, supply a `person_id,state,age` CSV (derived from the AoU ZIP3 geography) via
`PAN_AOU_STATE_CSV`; the builder uses it in preference. "Southern California" (Sunbelt) is
approximated by all of California; Alaska/Hawaii/DC are unassigned (controls).

**Interpretation caveat:** within EUR-unrelated samples, residualized on 10 PCs, these largely capture
**residual fine-scale geographic genetic structure** (the PCs already absorb the major geographic
gradients) — read them as geography/structure/migration signals, not trait biology. Political
clusters via work-state are a coarse proxy for where people actually live and vote.

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
| ZIP3 socioeconomic context (§11b.1) | 7 | 0 | 0 | 7 | 7 |
| Cognitive / EA-proxy external scores | 10 | 0 | 0 | 10 | 10 |
| Validated composite scores (§11c): scales + BFI-2 Big Five + neighborhood/walkability/hunger + PCL + well-being | 29 | 0 | 29 | 0 | 29 |
| Derived psychiatric phenotypes (§11d) | 8 | 8 | 0 | 0 | 8 |
| Acculturation index (§11e) | 1 | 0 | 1 | 0 | 1 |
| Geographic / political state clusters (§11f) | 12 | 12 | 0 | 0 | 12 |
| Wearable (Fitbit) phenotypes incl. chronotype (§10b) | 6 | 0 | 0 | 6 | 6 |
| **TOTAL** | **~852** | **3036** | **408** | **86** | **3531** |

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
metadata/ea_proxy_feature_sources.tsv     supplemental live v9 source questions from SES-EA/direct-XGB
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
