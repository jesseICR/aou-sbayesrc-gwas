# Cognitive Task Factor Score Spec

This document specifies the repo-tracked command for creating task-specific
Exploring the Mind cognitive scores among participants with a valid SES-EA proxy
phenotype.

Command:

```bash
bash run_etm_cog_task_factors.sh --stage-aggregate
```

The command does not run GWAS. It scores the already-defined SES-EA proxy
phenotype cohort from:

```text
sbayesrc_genotypes/regenie_input/ses_ea_proxy/all_scores.tsv
```

Both `role=oof` and `role=applied` samples are included. The scorer preserves
`role`, `fold_id`, `teacher_z`, `ses_ea_proxy_z`, and `ea_years` for diagnostics.

## Outputs

Primary/recommended task scores:

- `dd_patience_z_age_sex`, sourced from official `-lnk`
- `gradcpt_perf_z_age_sex`, sourced from PC1 of `dprime`,
  `-log(cv_rtc)`, and `-log(median_rtc)`
- `flanker_efficiency_z_age_sex`, sourced from the Flanker efficiency score
  selected by the simple-score rule

Individual-level outputs stay in the local on-platform scratch directory:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/
```

Aggregate diagnostics can be staged to:

```text
/home/jupyter/workspace/workspace-bucket/sbayesrc_genotypes/regenie_input/ses_ea_proxy/scrap/etm_cog_task_factors/
```

The user-facing wide table is:

```text
etm_cog_task_factors_recommended_wide.tsv
```

It has one row per SES-EA proxy-cohort sample, with missing values for ETM tasks
the participant did not take.

## Final Score Decisions

The scorer still computes factor and sensitivity diagnostics, but the
programmatic recommended outputs are fixed to the three-score set above.

Delay Discounting:

- Use official `-lnk`.
- `lnk` is the natural log of the mean discounting factor across the four
  substantive delay conditions:

  ```text
  lnk = log(mean(exp(two_weeks_lnk),
                 exp(one_month_lnk),
                 exp(one_year_lnk),
                 exp(ten_years_lnk)))
  ```

- The four-delay factor is cleaner as a latent-condition score, but the official
  `-lnk` aggregate was stronger against both the SES-EA proxy and teacher-label
  diagnostics and avoids overbuilding a score that already exists in the task
  output.

GradCPT:

- Use PC1 of `dprime`, `-log(cv_rtc)`, and `-log(median_rtc)`.
- One-factor FA is attempted, but PCA is the primary fallback when FA does not
  converge cleanly.
- The component-accuracy sensitivity score is not recommended when no-go
  accuracy has a weak loading.

Flanker:

- Attempt the one-factor efficiency/interference candidate for diagnostics.
- If the blend fails loading rules, use the predeclared efficiency score.
- Interference remains diagnostic unless its split indicators form a stable
  coherent score.

## Cohort, Age, And Sex

The scoring cohort is the final SES-EA proxy cohort in `all_scores.tsv`. The
script does not independently rebuild European ancestry, EA-label, sex, age, or
sample-QC filters.

Age at test is computed for each task sitting as:

```text
DATE_DIFF(test_start_date_time, person.birth_datetime, DAY) / 365.25
```

`birth_datetime` comes from the main CDR `person` table, and
`test_start_date_time` comes from the ETM task metadata.

Sex uses the confirmed genetic-sex covariate already produced by the genetics
pipeline. In this repo, `base_covar.txt` contains `sex_c`, where:

```text
sex_c = confirmed_genetic_sex_01 - 0.5
```

The scorer expects `sex_c` to contain only `-0.5` and `0.5`; any other value is
reported as an error.

## Valid Sitting Rules

Use the first valid completed sitting per `person_id` and task. Sort valid
sittings by `test_start_date_time`, then `sitting_id` as a deterministic
tie-breaker.

Delay Discounting valid sitting:

```text
flag_median_rt == 0
flag_catch_trials == 0
test_restarted == false when available
```

GradCPT valid sitting:

```text
flag_trial_flags == 0
flag_non_response == 0
flag_omission_error_rate == 0
test_restarted == false when available
```

Flanker valid sitting:

```text
flag_accuracy == 0
flag_trial_flags == 0
test_restarted == false when available
```

`any_timeouts` is not a primary exclusion. It is reported and can be used for
sensitivity analyses.

Delay Discounting catch trials are attention/validity checks. They are used for
QC through `flag_catch_trials`; `catch_score` is not a factor indicator.

## Score Inputs

### Delay Discounting

Recommended score:

```text
-lnk
```

Diagnostic four-delay factor inputs:

```text
-two_weeks_lnk
-one_month_lnk
-one_year_lnk
-ten_years_lnk
```

Validation-only:

```text
score
catch_score
mean_rt
median_rt
sd_rt
```

Do not use `log(lnk)`. The `lnk` fields are already natural-log discounting
parameters and can be negative.

### GradCPT

Primary PC inputs:

```text
dprime
-log(cv_rtc)
-log(median_rtc)
```

If `-log(median_rtc)` loads opposite `dprime` and `-log(cv_rtc)`, drop
`median_rtc` from the primary PC and report it as a speed diagnostic.

Sensitivity component factor:

```text
logit(go_accuracy)
logit(nogo_accuracy)
-log(cv_rtc)
-log(median_rtc)
```

Validation-only:

```text
go_accuracy
nogo_accuracy
score
crit
```

The primary PC does not include both `dprime` and its component accuracy
rates. `score` is not included with `nogo_accuracy` because it duplicates no-go
accuracy.

### Flanker

One-factor candidate inputs:

```text
log(rcs_congruent + eps)
log(rcs_incongruent + eps)
-accuracy_interference
-median_rt_interference
```

Validation-only:

```text
-rcs_interference
score
accuracy
median_rtc
```

Do not include `rcs_interference` in the same primary factor as
`rcs_congruent` and `rcs_incongruent`, because `rcs_interference` is derived
from the condition RCS values.

If the one-factor Flanker score is unstable or incoherent, use predeclared split
scores:

```text
flanker_efficiency:
  log(rcs_congruent + eps)
  log(rcs_incongruent + eps)
  score = mean of z-scored aligned indicators

flanker_interference:
  -accuracy_interference
  -median_rt_interference
  score = mean of z-scored aligned indicators if both align
  otherwise report the indicators separately and mark the composite unstable
```

Validate `flanker_efficiency` against `score`, and validate
`flanker_interference` against `-rcs_interference`.

## Transformations

Use fixed transformations, not normality-test-selected transformations.

Rules:

- Delay `lnk` fields: reverse sign only.
- RT and RT variability: `log(x)`, then reverse sign when lower is better.
- Accuracy/proportion variables: clipped logit after clipping to `[0.001, 0.999]`.
- RCS: `log(x + eps)`.
- Interference differences: reverse sign only.

For RCS:

```text
eps = 0.5 * min positive RCS in the task reference set after QC
```

The script saves `eps` in diagnostics.

After transformation, winsorize each indicator at the 0.5th and 99.5th
percentiles, then z-score it. Z-scoring is required because factor analysis and
PCA are scale-sensitive.

## Missingness

Primary rule:

```text
Require all primary indicators for a task to be non-missing after QC and transform.
```

Sensitivity rule:

```text
Mean imputation is allowed only after transform/winsorization/z-scoring,
only for non-core indicators, and must write imputed_indicator_count.
```

The current implementation uses the strict complete-case primary rule.

## Factor Model

Primary model:

```text
sklearn.decomposition.FactorAnalysis(n_components=1)
```

PCA is used only as fallback or diagnostic.

Factor construction is unsupervised with respect to SES-EA and teacher labels.
Do not select a cognitive score because it has a higher SES-EA or teacher
correlation.

The factor score is oriented so higher means better performance, or less
temporal discounting for Delay Discounting.

## Age/Sex Norming

Do not residualize primary indicators before factor fitting. Compute the task
factor first, then residualize/norm the final factor:

```text
factor_raw ~ sex_c + age_at_test + age_at_test^2
```

The final score is the z-scored residual from that model.

Diagnostic sensitivity:

```text
Fit the same factor model on age/sex-residualized indicators.
Compare that score with the primary final-only-adjusted score.
```

If the correlation is high, final-only age/sex norming is adequate. If not,
age/sex structure materially changes the factor definition and must be reported.

## Redundancy And Acceptance Rules

Before fitting a primary factor, compute the direction-aligned indicator
correlation matrix. Drop near-duplicates with:

```text
abs(r) > 0.95
```

Use fixed construct priority, not SES-EA/teacher correlations, to decide which
indicator to keep.

Accept a one-factor score only if:

- expected-good indicators load in the expected direction
- major loadings are not weak, default `abs(loading) >= 0.20`
- no obvious split creates opposite-signed construct dimensions
- the factor correlates sensibly with the official/simple score

If a factor correlates `>= 0.95` with the official/simple score, interpret the
simple score as primary regardless of SES-EA/teacher correlations. Keep the
factor as a sensitivity score.

Official/simple comparisons:

- Delay Discounting: `-lnk` and `score`
- GradCPT: `dprime` and `score`
- Flanker efficiency/performance: `score`
- Flanker interference: `-rcs_interference`

## Diagnostics

The command writes:

- factor scores
- selected and dropped indicators
- indicator missingness/exclusions
- repeat-sitting counts
- RCS `eps` parameters
- direction-aligned indicator correlation matrices
- factor-analysis loadings
- PCA diagnostic loadings
- factor-vs-simple-score correlations
- SES-EA proxy/teacher correlations, combined and by `role`
- age/sex model coefficients
- age/sex-residualized-indicator sensitivity
- device/language/version sensitivity

Device/language/version diagnostics use:

```text
score_raw ~ age + age^2 + sex_c + response_device + touch + test_language + test_version
```

The reported diagnostic is the correlation between the age/sex-only score and
the score additionally adjusted for administration metadata. Small categories
are collapsed to avoid tiny-cell diagnostic tables.
