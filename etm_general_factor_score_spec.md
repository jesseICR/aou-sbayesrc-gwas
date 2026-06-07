# ETM General Cognitive/Performance Factor Score

This spec defines a downstream command for creating one ETM general
cognitive/performance factor score from the already-built task-specific ETM
scores. It is intended for SES-EA proxy-phenotyped participants who have at
least one valid ETM task score.

The retained primary score is a point estimate of a three-domain ETM general
factor from Delay Discounting, GradCPT, and Flanker. When Emotional Recognition
is available, the command also writes a four-domain ETM general-factor score for
comparison and possible downstream use. Neither score should be described as
definitive global IQ or full psychometric `g`.

## Command

Preferred command:

```bash
bash run_etm_g_from_task_scores.sh --stage-aggregate
```

Optional follow-on from task scoring:

```bash
bash run_etm_cog_task_factors.sh --stage-aggregate --make-etm-g
```

This step does not query ETM tables, run GWAS, commit or push, use SES-EA or
teacher labels to choose loadings, residualize the final score for age/sex,
impute missing task scores, or output posterior uncertainty columns.

## Inputs

Primary task-score input:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_cog_task_factors_recommended_wide.tsv
```

Authoritative proxy cohort input:

```text
sbayesrc_genotypes/regenie_input/ses_ea_proxy/all_scores.tsv
```

Confirmed genetic sex and YOB-like covariates are read from:

```text
sbayesrc_genotypes/regenie_input/ses_ea_proxy/base_covar.txt
```

The retained three-domain task inputs, in fixed order, are:

```text
DD       = dd_patience_z_age_sex
GradCPT  = gradcpt_perf_z_age_sex
Flanker  = flanker_efficiency_z_age_sex
```

The optional four-domain model adds:

```text
EmoRecog = emorecog_perf_z_age_sex
```

All inputs are already oriented so higher means better performance, or less
discounting/more patience for Delay Discounting. The ETM-g step applies no
additional age/sex residualization because these task inputs are already
task-level age/sex-normalized.

If a future task-score run emits a coherent `flanker_perf_z_age_sex`, the
general-factor command can use it with:

```bash
--flanker-input flanker_perf_z_age_sex
```

The default remains `flanker_efficiency_z_age_sex`. Flanker interference is not
included in the primary ETM-g model because it would overweight Flanker relative
to DD and GradCPT.

## Measurement Model

The model is estimated only in all-three-task complete cases:

```text
complete_case_ref_3 = nonmissing DD and nonmissing GradCPT and nonmissing Flanker
complete_case_ref_4 = nonmissing DD and nonmissing GradCPT and nonmissing Flanker and nonmissing EmoRecog
```

The complete-case reference sample is used to:

- re-center and re-scale the selected task inputs
- fit one-factor Gaussian factor analysis
- estimate loadings and uniquenesses
- compute observed-task-pattern scoring weights
- z-scale final point estimates

For each task input `j`, compute complete-case mean and sample SD, then
standardize every participant:

```text
x_ij = (task_ij - mu_j_cc) / sd_j_cc
```

Fit:

```python
FactorAnalysis(n_components=1, random_state=0)
```

on the complete-case standardized matrix. Extract:

```text
lambda_vec_3 = [lambda_DD, lambda_GradCPT, lambda_Flanker]
lambda_vec_4 = [lambda_DD, lambda_GradCPT, lambda_Flanker, lambda_EmoRecog]
```

Orient the factor so higher means better ETM performance:

```text
if sum(lambda_vec) < 0: lambda_vec = -lambda_vec
```

Primary `etm_g_z` is accepted only if all three-domain loadings are positive,
at least two loadings are `>= 0.20`, and the complete-case task correlations do
not show a severe negative contradiction to a positive manifold. If the
three-domain model fails these rules, the command withholds default `etm_g_z`;
`--force-three-domain-g` writes `etm_g_z_forced` for debugging instead.

The four-domain `etm_g4_z` is accepted only if all four loadings are positive,
at least three loadings are `>= 0.20`, and the all-four complete-case task
correlations do not show a severe negative contradiction to a positive manifold.

## Scoring Missing-Task Patterns

All SES-EA proxy-cohort rows are preserved. Participants with one, two, or
three observed three-domain task scores are scored for `etm_g_z`. For
`etm_g4_z`, participants with one to four observed task scores among DD,
GradCPT, Flanker, and Emotional Recognition are scored. Participants with no
observed task score retain missing scores.

For participant `i`, let `O_i` be the observed task subset. For that subset:

```text
Sigma_O = lambda_O lambda_O' + diag(psi_O)
b_O     = Sigma_O^{-1} lambda_O
g_hat_i = x_iO' b_O
```

This is the regression factor-score point estimate. One-task participants are
intentionally shrunk toward zero because a single domain provides less
information about the general factor. Missing task scores are not imputed, and
scores are not z-scaled separately by missingness pattern.

Final scaling uses the complete-case `g_hat` distribution:

```text
etm_g_z_i = (g_hat_i - mean(g_hat_cc)) / sample_sd(g_hat_cc)
```

## Outputs

Local individual-level outputs:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_general_factor/
  etm_general_factor_scores_wide.tsv
  etm_general_factor_scores_scored_only.tsv
```

Required score columns include:

```text
IID
person_id
role
fold_id
etm_g_z
etm_g_hat
accepted_three_domain_g
n_tasks_observed
task_pattern
has_dd
has_gradcpt
has_flanker
etm_g4_z
etm_g4_hat
accepted_four_domain_g
n_tasks_observed_four_domain
task_pattern_four_domain
has_four_emorecog
```

The command also writes the selected task scores, `flanker_input_source`,
SES-EA proxy diagnostics columns, `sex_c`, and available YOB-like columns. It
does not write posterior SD, posterior variance, or standard-error columns.

Local aggregate diagnostics:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_general_factor/diagnostics/
```

With `--stage-aggregate`, diagnostics only are copied to:

```text
/home/jupyter/workspace/workspace-bucket/sbayesrc_genotypes/regenie_input/ses_ea_proxy/scrap/etm_cog_task_factors/etm_general_factor/
```

Individual-level score files are not staged by default.

## Diagnostics

The scorer writes:

```text
etm_g_reference_standardization.tsv
etm_g_task_missingness_counts.tsv
etm_g_task_pattern_counts.tsv
etm_g_complete_case_task_correlations.tsv
etm_g_pair_available_task_correlations.tsv
etm_g_fa_loadings.tsv
etm_g_uniquenesses.tsv
etm_g_model_implied_correlations.tsv
etm_g_residual_correlations.tsv
etm_g_scoring_weights_by_pattern.tsv
etm_g_score_distribution_by_pattern.tsv
etm_g_age_sex_validation.tsv
etm_g_external_validation_correlations.tsv
etm_g_comparison_scores_summary.tsv
etm_attention_exec_diagnostic_summary.tsv
etm_g_factor_summary.tsv
etm_g_three_vs_four_comparison.tsv
etm_g4_reference_standardization.tsv
etm_g4_complete_case_task_correlations.tsv
etm_g4_pair_available_task_correlations.tsv
etm_g4_fa_loadings.tsv
etm_g4_uniquenesses.tsv
etm_g4_scoring_weights_by_pattern.tsv
etm_g4_external_validation_correlations.tsv
etm_g_reproducibility_params.tsv
```

Age/sex and SES-EA/teacher/EA-years diagnostics are descriptive only. They do
not alter the phenotype, loadings, scoring weights, task inclusion, or fallback
rules.

The predeclared GradCPT/Flanker diagnostic score is `etm_attention_exec_z`. It
uses the same complete-case loading/scoring logic on GradCPT and Flanker only.
It is not the primary phenotype unless the three-domain model is not accepted.

## Reproducibility

The command saves enough metadata to reproduce the score:

- input paths, sizes, mtimes, and hashes when cheap to compute
- task column order and Flanker source
- complete-case reference N
- complete-case task means and SDs
- FA random state, loadings, and uniquenesses
- uniqueness floor use
- orientation sign and acceptance flag
- pattern-specific scoring weights
- complete-case `g_hat` mean and SD
- software versions

The command is idempotent. If outputs and reproducibility parameters already
match the current inputs, it skips. If existing outputs do not match, it fails
unless `--force` is passed.
