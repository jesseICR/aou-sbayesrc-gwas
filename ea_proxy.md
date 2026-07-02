# EA Proxy, Income, ETM Cognitive Scores, and Final GradCPT/Flanker Proxy

This document is the methods/runbook record for downstream phenotype workflows
that build on the core AoU SBayesRC/REGENIE genotype setup in `README.md`. The
workflows stay on-platform and regenerate AoU-derived individual-level files in
the workspace bucket or local scratch. Generated score tables, model JSONs,
phenotype files, covariates, and REGENIE outputs are intentionally not tracked
in git.

## Downstream phenotype command overview

These commands are deliberately separate from `get_genotypes.sh`. They assume
the main pipeline has completed through Step 13, but they do not require the
optional height GWAS to have run.

Cheap setup-only checks for simple survey phenotypes:

```bash
bash run_ea_gwas.sh --setup-only
bash run_income_gwas.sh --setup-only
```

Full GWAS submissions for those simple survey phenotypes:

```bash
bash run_ea_gwas.sh
bash run_income_gwas.sh
```

Final selected cdrv9 g-EA proxy GWAS:

```bash
bash run_g_ea_proxy_v9_pipeline.sh --preflight-only --force-final
bash run_g_ea_proxy_v9_pipeline.sh --smoke --force-final
bash run_g_ea_proxy_v9_pipeline.sh --skip-setup --chroms 1-22 --force-final
```

That workflow uses the final fold-safe no-teacher phenotype:

```text
regenie_input/g_ea_proxy_sbayesrc7m/phen.txt
gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z
```

The final GWAS covariates are `sex_c` and `PC1_AVG` through `PC10_AVG`.
`yob_c` and `yob_c_sex_c_inter` are not included. After construction, the final
phenotype had only small age/year-of-birth association, so the GWAS covariate
set was kept to sex plus ancestry PCs.

The final GWAS writes lightweight summary files here:

```text
regenie_output/g_ea_proxy_sbayesrc7m_gwas/lightweight/
```

The simpler EA and income GWAS commands use:

```text
Samples:
  our classified European IIDs
  AND confident genetic sex in genetic_sex/sex_covar.txt
  AND not in the Step 12 identical-component size >=3 exclusion list
  AND has a codeable trait response
  AND age at the relevant The Basics survey response >= 26 years

Genotypes:
  Step 1: gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1
  Step 2: gwas_genotypes/step2_wgs_pfiles/chr{1..22}

REGENIE:
  quantitative trait
  rank-inverse normal transform enabled by default
  lightweight per-chromosome outputs created by default
```

The shared REGENIE runner names Step 1 and Step 2 output files with the GWAS
output name as the prefix. For example, EA outputs are
`step1/ea_gwas_step1_pred.list` and
`step2/chr1/chr1_ea_gwas.regenie.gz`; income outputs are
`step1/income_gwas_step1_pred.list` and
`step2/chr1/chr1_income_gwas.regenie.gz`.

Future runs also write compact outputs:

```text
regenie_output/<trait>/lightweight/chr1.<trait>.regenie_lite.tsv.gz
regenie_output/<trait>/lightweight/regenie_lite.summary.tsv
```

The lightweight columns are `rsid`, `allele1`, `a1freq`, `n`, `beta`, `se`,
and `log10p`.

The educational-attainment phenotype is AoU The Basics question `1585940`,
"Education Level: Highest Grade." `PMI: Skip` and `PMI: Prefer Not To Answer`
are treated as missing. Samples younger than `GWAS_MIN_AGE_AT_SURVEY=26` at
the selected survey response are excluded so participants who may not have
completed education yet are not included.

| AoU answer concept | AoU answer | EA years |
|---:|---|---:|
| `1585941` | Highest Grade: Never Attended | 9.0 |
| `1585942` | Highest Grade: One Through Four | 9.0 |
| `1585943` | Highest Grade: Five Through Eight | 9.0 |
| `1585944` | Highest Grade: Nine Through Eleven | 10.0 |
| `1585945` | Highest Grade: Twelve Or GED | 13.0 |
| `1585946` | Highest Grade: College One to Three | 15.0 |
| `1585947` | Highest Grade: College Graduate | 18.0 |
| `1585948` | Highest Grade: Advanced Degree | 20.0 |

EA covariates:

```text
yob_c                  = fractional_year_of_birth - mean(fractional_year_of_birth)
sex_c                  = sex_01 - 0.5
yob_c_sex_c_inter      = yob_c * sex_c
PC1_AVG ... PC10_AVG
```

Fractional year of birth is computed from `person.birth_datetime`; for example
July 1, 1950 is approximately `1950.5`.

The household-income phenotype is AoU The Basics question `1585375`, "Income:
Annual Income." `PMI: Skip` and `PMI: Prefer Not To Answer` are treated as
missing. Samples younger than `GWAS_MIN_AGE_AT_SURVEY=26` at the selected
survey response are excluded so early-career participants are not included.
Values are annual household income in thousands of dollars.

| AoU answer concept | AoU answer | income_k |
|---:|---|---:|
| `1585376` | Annual Income: less 10k | 5.0 |
| `1585377` | Annual Income: 10k 25k | 17.5 |
| `1585378` | Annual Income: 25k 35k | 30.0 |
| `1585379` | Annual Income: 35k 50k | 42.5 |
| `1585380` | Annual Income: 50k 75k | 62.5 |
| `1585381` | Annual Income: 75k 100k | 87.5 |
| `1585382` | Annual Income: 100k 150k | 125.0 |
| `1585383` | Annual Income: 150k 200k | 175.0 |
| `1585384` | Annual Income: more 200k | 250.0 |

Income covariates:

```text
age_c                  = age_at_income_survey - mean(age_at_income_survey)
age_c_sq               = age_c^2
sex_c                  = sex_01 - 0.5
age_c_sex_c_inter      = age_c * sex_c
PC1_AVG ... PC10_AVG
```

The setup scripts write workspace-specific sample counts and answer counts to
`{ea,income}_gwas.summary.tsv` and `{ea,income}_answer_counts.tsv`.

### ses_ea_proxy primary setup and scoring

`run_ses_ea_proxy_gwas.sh` builds the primary SES-EA proxy scores and the
matching REGENIE input files. Setup/scoring is the default behavior; it does
not submit REGENIE unless `--run-gwas` is passed explicitly.

Recommended run order for the final selected phenotype:

```bash
# 1. Build genotype/sample-QC/PCA/sex/genotype inputs.
nohup bash get_genotypes.sh > logs/run.log 2>&1 &

# 2. Build the revised SES-EA proxy scores and REGENIE input files.
# The sixth/applied model excludes fit_pca relatives of the applied cohort.
SES_EA_PROXY_GWAS_INPUT_NAME=ses_ea_proxy_v2_kinholdout \
SES_EA_PROXY_GWAS_OUTPUT_NAME=ses_ea_proxy_v2_kinholdout \
bash run_ses_ea_proxy_gwas.sh --setup-only

# 3. Build ETM cognitive task scores from the revised proxy-score outputs.
SES_EA_PROXY_GWAS_INPUT_NAME=ses_ea_proxy_v2_kinholdout \
bash run_etm_cog_task_factors.sh --stage-aggregate

# 4. Fine-tune the saved SES-EA proxy boosters toward GradCPT+Flanker.
SES_EA_PROXY_GWAS_INPUT_NAME=ses_ea_proxy_v2_kinholdout \
bash run_gradcpt_flanker_finetuned_ea_proxy.sh --stage-aggregate

# 5. Train the direct scratch XGBoost GradCPT/Flanker survey proxy.
SES_EA_PROXY_GWAS_INPUT_NAME=ses_ea_proxy_v2_kinholdout \
bash run_gradcpt_flanker_direct_xgb_proxy.sh --stage-aggregate

# 6. Build and run the final g-EA proxy GWAS.
bash run_g_ea_proxy_v9_pipeline.sh --skip-setup --chroms 1-22 --force-final
```

The optional ETM general-factor command can still be run as a diagnostic, but
it is not required for the final selected phenotype:

```bash
SES_EA_PROXY_GWAS_INPUT_NAME=ses_ea_proxy_v2_kinholdout \
bash run_etm_g_from_task_scores.sh --stage-aggregate
```

During development we backfilled kinship holdout into an existing
`ses_ea_proxy_v2` output by refitting only the sixth/final applied model and
copying the five OOF fold models unchanged. That repair path was a one-off
workspace acceleration and is not part of the reproducible public pipeline. A
fresh run uses `run_ses_ea_proxy_gwas.sh --setup-only`, which now applies the
same kinship holdout inside the normal setup path.

```bash
bash run_ses_ea_proxy_gwas.sh --setup-only
```

This command is downstream of the genotype/sample-prep pipeline. Run
`get_genotypes.sh` first and let it complete through the European ancestry,
PCA, confirmed genetic sex, identical-component sample-QC, and final REGENIE
genotype-input steps. The proxy setup consumes these outputs:

```text
sbayesrc_genotypes/europeans/classified_european_iids.txt
sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
sbayesrc_genotypes/pca_eur/aou_projected.sscore
sbayesrc_genotypes/genetic_sex/sex_covar.txt
sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.fam
```

The goal of this step is to produce a non-genetic proxy for the
education-attainment teacher label and then review out-of-sample model
performance plus covariate correlations before deciding whether to run GWAS.
For that reason, setup stops after score generation by default; REGENIE is an
explicit opt-in follow-up.

Samples are restricted to participants who are classified European, have
confirmed genetic sex in `genetic_sex/sex_covar.txt`, are not in the
identical-component size `>=3` sample-QC exclusion list, have a codeable
education response from The Basics question `1585940`, and were at least 26
years old at that response. Confirmed genetic sex is also included as an
XGBoost model feature.

The proxy score uses a cross-fit design:

```text
1. Split eligible fit_pca_iids into 5 seeded folds.
2. For each fold, fit the EA residualization OLS and XGBoost model on the
   other four folds only, then predict the held-out fold.
3. Fit a sixth model for the applied cohort. This model starts from eligible
   `fit_pca_iids`, but excludes any fit-PCA sample with a direct KING edge to
   an applied-cohort sample at `KINSHIP >= 0.0441941`.
4. Apply the sixth model to eligible classified-European samples that were not
   in fit_pca_iids. Those applied samples are never used to train or tune the
   sixth model.
```

The EA teacher label is mapped to years from answers `1585941` through
`1585948`. The finalized mapping deliberately clamps the sparse lowest
education bins to 9 years and uses a slightly less compressed college scale:

| AoU answer concept | AoU answer | EA years |
|---:|---|---:|
| `1585941` | Highest Grade: Never Attended | 9.0 |
| `1585942` | Highest Grade: One Through Four | 9.0 |
| `1585943` | Highest Grade: Five Through Eight | 9.0 |
| `1585944` | Highest Grade: Nine Through Eleven | 10.0 |
| `1585945` | Highest Grade: Twelve Or GED | 13.0 |
| `1585946` | Highest Grade: College One to Three | 15.0 |
| `1585947` | Highest Grade: College Graduate | 18.0 |
| `1585948` | Highest Grade: Advanced Degree | 20.0 |

This mapping was chosen after inspecting categorical proxy and ETM
GradCPT+Flanker diagnostics anchored to Twelve/GED = 13 and Advanced Degree =
20. Those diagnostics suggested that the very low education bins are too small
and compressed to justify extreme year values in this cohort, while College
Graduate is closer to 18 than 17. In each cross-fit fold, EA years are
residualized on `yob_c`, `sex_c`, and `yob_c * sex_c` using only the
four-fifths training pool, then z-scored using that training-pool residual mean
and SD. The sixth model uses the same residualization procedure, fit only on
the kinship-clean final-model training subset.

The current `ses_ea_proxy_v2_kinholdout` setup run produced:

```text
Final eligible proxy cohort: 280,101
OOF / fit_pca samples:       252,774
Applied samples:              27,327
Feature columns:                 711
Feature hash: b171d4724414e2933df1cc8cc7b8f2834dded1005f3b60a23a2784d3ccf9f9c8
```

The final applied-model kinship holdout used the same KING threshold as the
third-degree PCA pruning step:

```text
KING threshold:                              0.0441941
Applied seed samples:                          27,327
Candidate fit_pca samples:                    252,774
KING edges at or above threshold:              87,885
Excluded fit_pca relatives of applied samples: 14,341
Final model training samples:                 238,433
```

This kinship holdout affects only the sixth/final model used for the applied
cohort. The OOF fold models still train on the other four folds of
`fit_pca_iids`, so their out-of-fold validation remains directly comparable to
earlier runs.

OOF correlations with the revised `teacher_z`, by fold:

| Model | Test N | Train N | Best rounds | Pearson r vs teacher_z | Spearman r |
|---|---:|---:|---:|---:|---:|
| fold 0 | 39,882 | 158,384 | 1300 | 0.6817 | 0.6464 |
| fold 1 | 39,478 | 158,788 | 1024 | 0.6758 | 0.6392 |
| fold 2 | 39,835 | 158,431 | 1290 | 0.6795 | 0.6434 |
| fold 3 | 39,607 | 158,659 | 1166 | 0.6828 | 0.6442 |
| fold 4 | 39,464 | 158,802 | 1367 | 0.6740 | 0.6410 |

Overall:

```text
OOF Pearson r vs teacher_z:      0.6788
OOF Spearman r vs teacher_z:     0.6429
Applied Pearson r vs teacher_z:  0.6901
Applied Spearman r vs teacher_z: 0.6648
```

Survey features come from The Basics, Lifestyle, Overall Health, Healthcare
Access & Utilization, Personal and Family Health History, Social Determinants
of Health, and Behavioral Health & Personality. BHP is read from the off-cycle
Mental Health / Well-Being CDR dataset; override `WORKSPACE_MHWB_CDR` if the
dataset name differs. The Washington Group disability items are sourced from
The Basics. ZIP3-derived socioeconomic features come from
`ds_zip_code_socioeconomic`; raw ZIP codes are not used. The highest-grade
education item itself is excluded from the XGBoost feature matrix because it
defines the teacher label. Personal and Family Health History includes
allowlisted mental-health/substance-use family-history indicators: the
family-condition question keeps None, ADHD, alcohol use disorder, drug use
disorder, and autism spectrum disorder, and the alcohol/drug relative-specific
questions keep self, parent, sibling, grandparent, son, and daughter
indicators. PMI missing answers such as Skip, Prefer Not, and Don't Know are
treated as missing/nonresponse, not as negative family history.

Feature extraction is multi-select safe. For each person/question, the setup
selects the latest survey timestamp, keeps all answer rows from that timestamp,
and one-hot encodes every retained `answer_concept_id` for nominal and
multi-select fields. Ordered Likert-style fields are encoded as a single
ordinal numeric feature. Continuous survey values are parsed as numeric
features.

Missing-data handling is explicit:

```text
Did not take survey:
  Survey item features remain NaN, took_<survey>=0, age_at_<survey>=NaN.

Took survey:
  took_<survey>=1 and age_at_<survey> is included as a model feature.

PMI: Skip (903096):
  Value remains NaN for native XGBoost missing routing and is counted
  separately; no nonresponse indicator is set.

PMI: Prefer Not To Answer (903079) / PMI: Don't Know (903087):
  Value remains NaN. Curated high-value fields also get a *_nonresponse
  indicator where valid answers are 0 and Prefer Not/Don't Know are 1.
  Not asked or did not take the survey remains NaN.

Lifestyle branching:
  Never/no parent answers set not-applicable downstream smoking, alcohol, and
  substance-use count/frequency fields to 0 rather than NaN.
```

XGBoost receives `DMatrix(missing=np.nan)`, so survey nonparticipation and
unanswered fields route through missing-native tree splits. The setup writes
audit tables for these choices:

```text
feature_manifest.resolved.tsv
missing_data_handling.tsv
pmi_missingness_counts.tsv
branch_recoding_summary.tsv
feature_missingness.tsv
feature_counts.tsv
```

#### Resolved SES-EA survey feature contract

The resolved `ses_ea_proxy_v2_kinholdout` feature manifest contains 213
included survey/area items before expansion into 711 XGBoost columns:

| Encoding | Included items | Meaning |
|---|---:|---|
| `one_hot` | 100 | Nominal or multi-select answers become one binary column per retained answer. |
| `ordinal` | 98 | Ordered answer concepts become one numeric ordered feature. |
| `numeric` | 12 | Free numeric responses, survey ages, sex, and area-SES values remain numeric. |
| `allowlisted_one_hot` | 3 | Curated PFHH multi-select items use a fixed answer allowlist. |

The table below is generated from `feature_manifest.resolved.tsv` and records
the included survey item/source and how it is encoded before XGBoost sees it.
The direct GradCPT/Flanker scratch-XGBoost benchmark uses this same contract,
then additionally adds education item `1585940` as revised numeric EA years and
one-hot education-response features.

| Question concept/code | Survey/source | Encoding | Item |
|---|---|---|---|
| `1703874` | Behavioral Health and Personality | `ordinal` | During the past 6 months How often do you have trouble wrapping up the final details of a project, once the challenging parts have been done? |
| `1703875` | Behavioral Health and Personality | `ordinal` | I am someone who has difficulty getting started on tasks. |
| `1703878` | Behavioral Health and Personality | `ordinal` | During the past 6 months How often do you have problems remembering appointments or obligations? |
| `1703881` | Behavioral Health and Personality | `ordinal` | I am someone who is original, comes up with new ideas. |
| `1703892` | Behavioral Health and Personality | `ordinal` | I am someone who tends to be disorganized. |
| `1703893` | Behavioral Health and Personality | `ordinal` | I am someone who worries a lot. |
| `1703896` | Behavioral Health and Personality | `ordinal` | I am someone who is emotionally stable, not easily upset. |
| `1703903` | Behavioral Health and Personality | `ordinal` | I am someone who is full of energy. |
| `1703904` | Behavioral Health and Personality | `ordinal` | I am someone who is compassionate, has a soft heart. |
| `1703911` | Behavioral Health and Personality | `ordinal` | I am someone who is sometimes rude to others. |
| `1703912` | Behavioral Health and Personality | `ordinal` | I am someone who is reliable, can always be counted on. |
| `1703913` | Behavioral Health and Personality | `ordinal` | During the past 6 months How often do you fidget or squirm with your hands or feet when you have to sit down for a long time? |
| `1703914` | Behavioral Health and Personality | `ordinal` | During the past 6 months How often do you have difficulty getting things in order when you have to do a task that requires organization? |
| `1703916` | Behavioral Health and Personality | `ordinal` | During the past 6 months When you have a task that requires a lot of thought, how often do you avoid or delay getting started? |
| `1703918` | Behavioral Health and Personality | `ordinal` | I am someone who is dominant, acts as a leader. |
| `1703919` | Behavioral Health and Personality | `ordinal` | I am someone who has little interest in abstract ideas. |
| `1703925` | Behavioral Health and Personality | `ordinal` | During the past 6 months How often do you feel overly active and compelled to do things, like you were driven by a motor? |
| `1703926` | Behavioral Health and Personality | `ordinal` | I am someone who tends to be quiet. |
| `1703929` | Behavioral Health and Personality | `ordinal` | I am someone who tends to feel depressed, blue. |
| `1703930` | Behavioral Health and Personality | `ordinal` | I am someone who is fascinated by art, music, or literature. |
| `1703932` | Behavioral Health and Personality | `ordinal` | I am someone who assumes the best about people. |
| `43528660` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To General Doctor |
| `43528661` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Medical Specialist |
| `43528662` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Dental Care |
| `43528663` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Emergency Care |
| `43528664` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Healthcare Provider |
| `43528665` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Bought Rx From Other Country |
| `43528666` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Alternative Therapies |
| `43529903` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Child Care |
| `43529904` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Elderly Care |
| `43529905` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Time Off Work |
| `43529906` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Transportation |
| `43529973` | Healthcare Access & Utilization | `one_hot` | Health Advice: Nurse Practitioner, Physician Assistant, or Midwife Visits |
| `43529974` | Healthcare Access & Utilization | `one_hot` | Health Advice: Dentist or Orthodontist Visits |
| `43529975` | Healthcare Access & Utilization | `one_hot` | Health Advice: OB/GYN Visits |
| `43529976` | Healthcare Access & Utilization | `one_hot` | Health Advice: Medical Specialist Visits |
| `43529977` | Healthcare Access & Utilization | `one_hot` | Health Advice: Mental Health Professional Visits |
| `43529978` | Healthcare Access & Utilization | `one_hot` | Health Advice: Traditional Healer Visits |
| `43530268` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Rural Area |
| `43530399` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Chiropractor |
| `43530400` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Dentist or Orthodontist |
| `43530401` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To OB/GYN |
| `43530402` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Mental Health Professional |
| `43530403` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Eye Doctor |
| `43530404` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Nurse Practitioner |
| `43530405` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Physical Therapist, Speech Therapist, Respiratory Therapist, Audiologist, or Occupational Therapist |
| `43530406` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Podiatrist |
| `43530407` | Healthcare Access & Utilization | `ordinal` | Health Advice: Spoken To Traditional Healer |
| `43530408` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Eyeglasses |
| `43530409` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Follow-up Care |
| `43530410` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Mental Health Counseling |
| `43530411` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Prescription Medicines |
| `43530412` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Specialist |
| `43530413` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Lower Cost Rx To Save Money |
| `43530415` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Delayed Filling Rx To Save Money |
| `43530416` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Skipped Med To Save Money |
| `43530417` | Healthcare Access & Utilization | `ordinal` | Can't Afford Care: Took Less Med To Save Money |
| `43530418` | Healthcare Access & Utilization | `ordinal` | Insurance: Insurance Accepted |
| `43530437` | Healthcare Access & Utilization | `one_hot` | Health Advice: Asked For Opinion |
| `43530438` | Healthcare Access & Utilization | `one_hot` | Health Advice: Ease of Understanding |
| `43530439` | Healthcare Access & Utilization | `one_hot` | Health Advice: Respected By Provider |
| `43530557` | Healthcare Access & Utilization | `one_hot` | Can't Afford Care: Worried About Paying |
| `43530559` | Healthcare Access & Utilization | `one_hot` | Insurance: Healthcare Coverage |
| `43530562` | Healthcare Access & Utilization | `one_hot` | Health Advice: Place for Health Advice |
| `43530583` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Can't Afford Co-pay |
| `43530584` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Had To Pay Out Of Pocket |
| `43530585` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Deductible Too High |
| `43530588` | Healthcare Access & Utilization | `one_hot` | Health Advice: General Doctor Visits |
| `43530589` | Healthcare Access & Utilization | `one_hot` | Health Advice: Chiropractor Visits |
| `43530590` | Healthcare Access & Utilization | `one_hot` | Health Advice: Podiatrist Visits |
| `43530591` | Healthcare Access & Utilization | `one_hot` | Health Advice: Eye Doctor Visits |
| `43530592` | Healthcare Access & Utilization | `one_hot` | Health Advice: Physical Therapist, Speech Therapist, Respiratory Therapist, Audiologist, or Occupational Therapist Visits |
| `43530593` | Healthcare Access & Utilization | `one_hot` | Health Advice: What Kind Of Place |
| `43530594` | Healthcare Access & Utilization | `ordinal` | Delayed Medical Care: Nervous |
| `43530595` | Healthcare Access & Utilization | `one_hot` | Health Advice: Spoken To Professional |
| `1585636` | Lifestyle | `one_hot` | Recreational Drug Use: Which Drugs Used |
| `1585650` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Marijuana 3 Month Use |
| `1585656` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Cocaine 3 Month Use |
| `1585668` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Other Stimulant 3 Month Use |
| `1585674` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Inhalant 3 Month Use |
| `1585680` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Sedative 3 Month Use |
| `1585686` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Hallucinogen 3 Month Use |
| `1585692` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Street Opioid 3 Month Use |
| `1585698` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Prescription Opioid 3 Month Use |
| `1585704` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Other 3 Month Use |
| `1585857` | Lifestyle | `ordinal` | Smoking: 100 Cigs Lifetime |
| `1585860` | Lifestyle | `ordinal` | Smoking: Smoke Frequency |
| `1585864` | Lifestyle | `numeric` | Smoking: Daily Smoke Starting Age |
| `1585867` | Lifestyle | `one_hot` | Smoking: Serious Quit Attempt |
| `1585870` | Lifestyle | `numeric` | Attempt Quit Smoking: Completely Quit Age |
| `1585873` | Lifestyle | `numeric` | Smoking: Number Of Years |
| `1586159` | Lifestyle | `numeric` | Smoking: Current Daily Cigarette Number |
| `1586162` | Lifestyle | `numeric` | Smoking: Average Daily Cigarette Number |
| `1586166` | Lifestyle | `ordinal` | Electronic Smoking: Electric Smoke Participant |
| `1586169` | Lifestyle | `ordinal` | Electronic Smoking: Electric Smoke Frequency |
| `1586174` | Lifestyle | `ordinal` | Cigar Smoking: Cigar Smoke Participant |
| `1586177` | Lifestyle | `ordinal` | Cigar Smoking: Current Cigar Frequency |
| `1586182` | Lifestyle | `ordinal` | Hookah Smoking: Hookah Smoke Participant |
| `1586185` | Lifestyle | `ordinal` | Hookah Smoking: Current Hookah Frequency |
| `1586190` | Lifestyle | `ordinal` | Smokeless Tobacco: Smokeless Tobacco Participant |
| `1586193` | Lifestyle | `ordinal` | Smokeless Tobacco: Smokeless Tobacco Frequency |
| `1586198` | Lifestyle | `ordinal` | Alcohol: Alcohol Participant |
| `1586201` | Lifestyle | `ordinal` | Alcohol: Drink Frequency Past Year |
| `1586207` | Lifestyle | `one_hot` | Alcohol: Average Daily Drink Count |
| `1586213` | Lifestyle | `ordinal` | Alcohol: 6 or More Drinks Occurrence |
| `903058` | Lifestyle | `one_hot` | Past 3 Month Use Frequency: Prescription Stimulant 3 Month Use |
| `1585711` | Overall Health | `ordinal` | Overall Health: General Health |
| `1585717` | Overall Health | `ordinal` | Overall Health: General Quality |
| `1585723` | Overall Health | `ordinal` | Overall Health: General Physical Health |
| `1585729` | Overall Health | `one_hot` | Overall Health: General Mental Health |
| `1585735` | Overall Health | `ordinal` | Overall Health: Social Satisfaction |
| `1585741` | Overall Health | `one_hot` | Overall Health: Everyday Activities |
| `1585748` | Overall Health | `one_hot` | Overall Health: Average Fatigue 7 Days |
| `1585754` | Overall Health | `ordinal` | Overall Health: General Social |
| `1585760` | Overall Health | `one_hot` | Overall Health: Emotional Problem 7 Days |
| `1585766` | Overall Health | `one_hot` | Overall Health: Medical Form Confidence |
| `1585772` | Overall Health | `one_hot` | Overall Health: Health Material Assistance |
| `1585778` | Overall Health | `one_hot` | Overall Health: Difficult Understand Info |
| `1585815` | Overall Health | `one_hot` | Overall Health: Outside Travel 6 Month |
| `1740660` | Personal and Family Health History | `one_hot` | Including yourself, who in your family has had attention-deficit/hyperactivity disorder (ADHD)? Select all that apply. |
| `43529217` | Personal and Family Health History | `allowlisted_one_hot` | Have you or anyone in your family ever been diagnosed with the following mental health or substance use conditions? Think only of the people you are related to by blood. Select all that apply. |
| `836827` | Personal and Family Health History | `allowlisted_one_hot` | Including yourself, who in your family has had alcohol use disorder? Select all that apply. |
| `836851` | Personal and Family Health History | `allowlisted_one_hot` | Including yourself, who in your family has had a drug use disorder? Select all that apply. |
| `40192381` | Social Determinants of Health | `one_hot` | In the last month, how often have you felt that you were unable to control the important things in your life? |
| `40192384` | Social Determinants of Health | `ordinal` | How much you agree or disagree that your neighborhood is safe? |
| `40192386` | Social Determinants of Health | `ordinal` | How much you agree or disagree that people in your neighborhood take good care of their houses and apartments? |
| `40192388` | Social Determinants of Health | `one_hot` | How often do you have someone to prepare your meals if you were unable to do it yourself? |
| `40192390` | Social Determinants of Health | `ordinal` | How often do you feel that you are unhappy being so withdrawn? |
| `40192396` | Social Determinants of Health | `one_hot` | In the last month, how often have you been angered because of things that were outside of your control? |
| `40192397` | Social Determinants of Health | `ordinal` | How often do you feel that there is no one you can turn to? |
| `40192398` | Social Determinants of Health | `ordinal` | How often do you feel left out? |
| `40192399` | Social Determinants of Health | `one_hot` | How often do you have someone who understands your problems? |
| `40192400` | Social Determinants of Health | `ordinal` | How much you agree or disagree that in your neighborhood people watch out for each other? |
| `40192401` | Social Determinants of Health | `one_hot` | How often do you feel deep inner peace or harmony? |
| `40192402` | Social Determinants of Health | `one_hot` | Think about the place you live. Do you have problems with any of the following? Select all that apply. |
| `40192404` | Social Determinants of Health | `ordinal` | How much you agree or disagree that you are always having trouble with your neighbors? |
| `40192410` | Social Determinants of Health | `one_hot` | My neighborhood has several free or low-cost recreation facilities, such as parks, walking trails, bike paths, recreation centers, playgrounds, public swimming pools, etc. Would you say that you... |
| `40192411` | Social Determinants of Health | `one_hot` | How much you agree or disagree that people in your neighborhood generally get along with each other? |
| `40192412` | Social Determinants of Health | `ordinal` | How much you agree or disagree that vandalism is common in your neighborhood? |
| `40192414` | Social Determinants of Health | `one_hot` | The crime rate in my neighborhood makes it unsafe to go on walks during the day. Would you say that you... |
| `40192415` | Social Determinants of Health | `one_hot` | How often do you feel that you are spiritually touched by the beauty of creation? |
| `40192417` | Social Determinants of Health | `one_hot` | How much you agree or disagree that people in your neighborhood share the same values? |
| `40192419` | Social Determinants of Health | `one_hot` | In the last month, how often have you felt confident about your ability to handle your personal problems? |
| `40192420` | Social Determinants of Health | `ordinal` | How much you agree or disagree that there is a lot of graffiti in your neighborhood? |
| `40192426` | Social Determinants of Health | `one_hot` | Within the past 12 months, were you worried whether the food you had bought just didn't last and you didn't have money to get more? |
| `40192431` | Social Determinants of Health | `one_hot` | There are facilities to bicycle in or near my neighborhood, such as special lanes, separate paths or trails, or shared use paths for cycles and pedestrians. Would you say that you... |
| `40192436` | Social Determinants of Health | `one_hot` | Many shops, stores, markets or other places to buy things I need are within easy walking distance of my home. Would you say that you... |
| `40192437` | Social Determinants of Health | `one_hot` | There are sidewalks on most of the streets in my neighborhood. Would you say that you... |
| `40192439` | Social Determinants of Health | `one_hot` | How often do you have someone to have a good time with? |
| `40192440` | Social Determinants of Health | `one_hot` | It is within a 10-15 minute walk to a transit stop (such as bus, train, trolley, or tram) from my home. Would you say that you... |
| `40192441` | Social Determinants of Health | `numeric` | In the last 12 months, how many times have you or your family moved from one home to another? Number of moves in past 12 months: |
| `40192442` | Social Determinants of Health | `one_hot` | How often do you have someone to help you if you were confined to bed? |
| `40192443` | Social Determinants of Health | `one_hot` | How often do you desire to be closer to or in union with God (or a higher power)? |
| `40192445` | Social Determinants of Health | `one_hot` | In the last month, how often have you felt that you were on top of things? |
| `40192446` | Social Determinants of Health | `one_hot` | How often do you have someone to love and make you feel wanted? |
| `40192449` | Social Determinants of Health | `one_hot` | In the last month, how often have you been able to control irritations in your life? |
| `40192452` | Social Determinants of Health | `one_hot` | In the last month, how often have you been upset because of something that happened unexpectedly? |
| `40192456` | Social Determinants of Health | `ordinal` | How much you agree or disagree that your neighborhood is clean? |
| `40192457` | Social Determinants of Health | `ordinal` | How much you agree or disagree that there is too much drug use in your neighborhood? |
| `40192458` | Social Determinants of Health | `one_hot` | What is the main type of housing in your neighborhood? |
| `40192462` | Social Determinants of Health | `one_hot` | In the last month, how often have you felt difficulties were piling up so high that you could not overcome them? |
| `40192463` | Social Determinants of Health | `one_hot` | How much you agree or disagree that people around here are willing to help their neighbor? |
| `40192469` | Social Determinants of Health | `ordinal` | How much you agree or disagree that there are lot of abandoned buildings in your neighborhood? |
| `40192470` | Social Determinants of Health | `one_hot` | How often do you go to religious meetings or services? |
| `40192471` | Social Determinants of Health | `one_hot` | How often do you feel God's (or a higher power's) love for you, directly or through others? |
| `40192475` | Social Determinants of Health | `one_hot` | How often do you find strength and comfort in your religion? |
| `40192476` | Social Determinants of Health | `ordinal` | How much you agree or disagree that there is too much alcohol use in your neighborhood? |
| `40192480` | Social Determinants of Health | `one_hot` | How often do you have someone to take you to the doctor if you need it? |
| `40192491` | Social Determinants of Health | `one_hot` | In the last month, how often have you felt nervous and "stressed"? |
| `40192492` | Social Determinants of Health | `one_hot` | The crime rate in my neighborhood makes it unsafe to go on walks at night. Would you say that you... |
| `40192493` | Social Determinants of Health | `ordinal` | How much you agree or disagree that there is a lot of crime in your neighborhood? |
| `40192494` | Social Determinants of Health | `ordinal` | How often do you feel that people are around you but not with you? |
| `40192498` | Social Determinants of Health | `one_hot` | How often do you feel God's (or a higher power's) presence? |
| `40192499` | Social Determinants of Health | `one_hot` | How much you agree or disagree that people in your neighborhood can be trusted? |
| `40192500` | Social Determinants of Health | `ordinal` | How much you agree or disagree that there are too many people hanging around on the streets near your home? |
| `40192501` | Social Determinants of Health | `ordinal` | How often do you feel isolated from others? |
| `40192504` | Social Determinants of Health | `ordinal` | How often do you feel that you are an outgoing person? |
| `40192506` | Social Determinants of Health | `one_hot` | In the last month, how often have you found that you could not cope with all the things that you had to do? |
| `40192507` | Social Determinants of Health | `ordinal` | How often do you feel lack companionship? |
| `40192511` | Social Determinants of Health | `one_hot` | How often do you have someone to help you with daily chores if you were sick? |
| `40192516` | Social Determinants of Health | `ordinal` | How often do you fell that you can find companionship when you want it? |
| `40192517` | Social Determinants of Health | `one_hot` | Within the past 12 months, were you worried whether your food would run out before you got money to buy more? |
| `40192522` | Social Determinants of Health | `ordinal` | How much you agree or disagree that your neighborhood is noisy? |
| `40192525` | Social Determinants of Health | `one_hot` | In the last month, how often have you felt that things were going your way? |
| `40192528` | Social Determinants of Health | `one_hot` | How often do you have someone to turn to for suggestions about how to deal with a personal problem? |
| `1585357` | The Basics | `one_hot` | Gender Identity: Sexuality Closer Description |
| `1585370` | The Basics | `one_hot` | Home Own: Current Home Own |
| `1585375` | The Basics | `one_hot` | Income: Annual Income |
| `1585386` | The Basics | `ordinal` | Insurance: Health Insurance |
| `1585389` | The Basics | `one_hot` | Health Insurance: Health Insurance Type |
| `1585402` | The Basics | `one_hot` | Living Situation: Current Living |
| `1585852` | The Basics | `ordinal` | Active Duty: Active Duty Serve Status |
| `1585879` | The Basics | `one_hot` | Living Situation: How Many Living Years |
| `1585886` | The Basics | `ordinal` | Living Situation: Stable House Concern |
| `1585889` | The Basics | `one_hot` | Living Situation: How Many People |
| `1585890` | The Basics | `one_hot` | Living Situation: People Under 18 |
| `1585892` | The Basics | `one_hot` | Marital Status: Current Marital Status |
| `1585899` | The Basics | `one_hot` | The Basics: Sexual Orientation |
| `1585952` | The Basics | `one_hot` | Employment: Employment Status |
| `1586135` | The Basics | `one_hot` | The Basics: Birthplace |
| `43528428` | The Basics | `one_hot` | Health Insurance: Insurance Type Update |
| `903573` | The Basics | `one_hot` | Disability: Deaf |
| `903574` | The Basics | `one_hot` | Disability: Blind |
| `903575` | The Basics | `one_hot` | Disability: Difficulty Concentrating |
| `903576` | The Basics | `one_hot` | Disability: Walking Climbing |
| `903577` | The Basics | `one_hot` | Disability: Dressing Bathing |
| `903578` | The Basics | `one_hot` | Disability: Errands Alone |
| `zip3_ses` | zip3_ses_map | `numeric` | deprivation_index |
| `zip3_ses` | zip3_ses_map | `numeric` | median_income |
| `zip3_ses` | zip3_ses_map | `numeric` | fraction_poverty |
| `zip3_ses` | zip3_ses_map | `numeric` | fraction_assisted_income |
| `zip3_ses` | zip3_ses_map | `numeric` | fraction_no_health_ins |
| `zip3_ses` | zip3_ses_map | `numeric` | fraction_vacant_housing |

The main scoring outputs are:

```text
regenie_input/ses_ea_proxy/oof_scores.tsv
regenie_input/ses_ea_proxy/applied_scores.tsv
regenie_input/ses_ea_proxy/all_scores.tsv
regenie_input/ses_ea_proxy/fold_metrics.tsv
regenie_input/ses_ea_proxy/applied_metrics.tsv
regenie_input/ses_ea_proxy/proxy_covariate_correlations.tsv
regenie_input/ses_ea_proxy/xgboost_model_manifest.tsv
regenie_input/ses_ea_proxy/xgboost_feature_columns.json
regenie_input/ses_ea_proxy/xgboost_models/fold_0.json
regenie_input/ses_ea_proxy/xgboost_models/fold_1.json
regenie_input/ses_ea_proxy/xgboost_models/fold_2.json
regenie_input/ses_ea_proxy/xgboost_models/fold_3.json
regenie_input/ses_ea_proxy/xgboost_models/fold_4.json
regenie_input/ses_ea_proxy/xgboost_models/final_model.json
regenie_input/ses_ea_proxy/phen.txt
regenie_input/ses_ea_proxy/covar.txt
regenie_input/ses_ea_proxy/ses_ea_proxy_gwas.summary.tsv
```

`fold_metrics.tsv` reports the correlation of each held-out fold's OOF score
with the fold-safe teacher label. `applied_metrics.tsv` reports the same
correlation for the sixth-model applied samples. `proxy_covariate_correlations.tsv`
reports correlations with `teacher_z`, EA years, `yob_c`, `sex_c`, and
PC1 through PC10 by fold, OOF overall, sixth-model applied samples, and the
combined scored set. The XGBoost model JSON files are saved with the feature
column hash recorded in `xgboost_model_manifest.tsv`; reload them with the
same feature order in `xgboost_feature_columns.json` if continuing training or
auditing predictions.

The generated result files should be read as follows:

```text
ses_ea_proxy_gwas.summary.tsv
  One-row-per-metric setup summary: sample filters, final scored sample
  counts, feature counts, model score scale, OOF performance, and applied-set
  performance.

fold_metrics.tsv
  One row for each OOF fold. These are the primary performance checks because
  every sample in a fold was predicted by a model that did not train on it.

applied_metrics.tsv
  Performance for the sixth model applied to classified-European samples that
  were outside fit_pca_iids.

proxy_covariate_correlations.tsv
  Pearson and Spearman correlations between the proxy score and teacher_z,
  EA years, yob_c, sex_c, and PC1 through PC10, split by fold, OOF overall,
  applied samples, and combined scored samples.

feature_importance.tsv
  Gain/cover importance from the sixth model fit on all eligible fit_pca_iids.

feature_manifest.resolved.tsv
  Per-question manifest with question_concept_id, survey, item name, resolved
  encoding, include/exclude status, analog notes, and implementation notes.

xgboost_model_manifest.tsv
  Saved booster inventory. It records each model file, role, fold id, feature
  column hash, boost rounds, prediction sample count, and the teacher
  residualization parameters used by that model.
```

The saved Booster files can be loaded later with XGBoost for auditing or
continued training:

```python
import json
import xgboost as xgb

with open("xgboost_feature_columns.json") as f:
    feature_columns = json.load(f)

booster = xgb.Booster()
booster.load_model("xgboost_models/final_model.json")
```

If continuing training, rebuild the feature matrix with the same column order
from `xgboost_feature_columns.json` and pass the loaded booster as
`xgb_model=` to `xgb.train`.

The cognitive-score command is downstream of this proxy setup: it reads
`regenie_input/ses_ea_proxy/all_scores.tsv` and
`regenie_input/ses_ea_proxy/base_covar.txt`.

### ETM cognitive task factor scoring

`run_etm_cog_task_factors.sh` builds task-specific Exploring the Mind cognitive
scores for the already-defined `ses_ea_proxy` phenotype cohort. It consumes
`regenie_input/ses_ea_proxy/all_scores.tsv`, so both out-of-fold proxy scores
and sixth-model applied proxy scores are included. This is a scoring and
diagnostic command only; it does not run GWAS.

```bash
bash run_etm_cog_task_factors.sh --stage-aggregate
```

The scorer reads Flanker, GradCPT, Delay Discounting, and Emotional Recognition
from `WORKSPACE_ETM_CDR`. In the default cdrv9 workflow, that is the main
`C2025Q4R6` CDR because the ETM task tables are included there. The older v8
workflow used the `C_V8_R2_offcycle_etm` dataset. Set `WORKSPACE_ETM_CDR` or
pass `--etm-dataset PROJECT.DATASET` if the dataset name differs. Age at each
task is computed from ETM `test_start_date_time` and the main CDR
`person.birth_datetime`; confirmed genetic sex comes from the `sex_c` column
already present in `regenie_input/ses_ea_proxy/base_covar.txt`.

For each task, the scorer uses a fixed order of operations:

```text
1. Apply task-specific invalid-performance QC flags.
2. Exclude restarted tests when test_restarted is available.
3. Keep the first valid sitting per person/task by test_start_date_time,
   then sitting_id.
4. Build transformed indicators, winsorize them, and z-score them.
5. Fit the task-specific FA/PCA/simple score on those transformed indicators.
6. Only after the raw task score is formed, residualize:
     raw_task_score ~ sex_c + age_at_test + age_at_test^2
   then z-score the residual.
```

The individual measurement indicators are not residualized for age or sex
before FA/PCA. Age and confirmed genetic sex are removed only after the raw task
score is formed. This keeps the factor/PCA loadings tied to the task
measurements themselves rather than to age/sex-residualized inputs.

The primary valid-sitting filters are:

| Task | Excluded when any primary condition is true |
|---|---|
| Delay Discounting | `flag_median_rt != 0`, `flag_catch_trials != 0`, or `test_restarted` is true |
| GradCPT | `flag_trial_flags != 0`, `flag_non_response != 0`, `flag_omission_error_rate != 0`, or `test_restarted` is true |
| Flanker | `flag_accuracy != 0`, `flag_trial_flags != 0`, or `test_restarted` is true |
| Emotional Recognition | `flag_median_rtc != 0`, `flag_same_response != 0`, `flag_trial_flags != 0`, or `test_restarted` is true |

`any_timeouts` is retained as a diagnostic rather than an exclusion by default.

The command writes both a long diagnostic score table and a one-row-per-sample
recommended-score table. The recommended scores are:

```text
dd_patience_z_age_sex          # sourced from official -lnk
gradcpt_perf_z_age_sex         # PC1 of dprime + RT consistency/speed
flanker_efficiency_z_age_sex   # Flanker efficiency source selected by diagnostics
emorecog_perf_z_age_sex        # Emotional Recognition PC1 of score + RT consistency/speed
```

The score-selection rationale is:

```text
Delay Discounting:
  Use official -lnk. It summarizes all four delay-specific log-k fields as
  log(mean(k)), where k = exp(delay_lnk). The four-delay factor is still
  computed for diagnostics, but -lnk was stronger against both the SES-EA proxy
  and teacher-label diagnostics and is the cleaner official aggregate.

GradCPT:
  Use the PC1 fallback from dprime, -log(cv_rtc), and -log(median_rtc). The
  one-factor FA optimizer did not converge cleanly because the three-indicator
  correlation pattern sits on an invalid common-factor boundary, but PC1 was
  coherent and stable. The component-accuracy sensitivity score is computed but
  not recommended because no-go accuracy had a weak loading.

Flanker:
  Use the official AoU Flanker score after age/sex/age^2 norming. The broader
  one-factor efficiency/interference blend failed the loading rule. The
  predeclared split-efficiency composite was coherent, but it correlated
  >=0.95 with the official AoU score, so the simple official score wins by the
  predeclared simple-score rule.

Emotional Recognition:
  Score the task with the GradCPT-analog PC1 of score, -log(cv_rtc), and
  -log(median_rtc), after fixed transforms, winsorization, and z-scoring. The
  combined PC1 is then residualized on sex_c + age_at_test + age_at_test^2 and
  z-scored. The older per-emotion rate-correct efficiency factor, per-emotion
  accuracy factor, and simple score are still written as diagnostics.
```

The exact GradCPT and Flanker construction used for the final phenotype is:

```text
GradCPT:
  valid sitting:
    flag_trial_flags == 0
    flag_non_response == 0
    flag_omission_error_rate == 0
    test_restarted is false when available

  indicators:
    dprime
    -log(cv_rtc)
    -log(median_rtc)

  preprocessing:
    require positive cv_rtc and median_rtc before log transforms
    winsorize each transformed indicator at 0.5th / 99.5th percentiles
    z-score transformed indicators

  score:
    attempt one-factor FA
    use PC1 because FA fell onto an unstable/common-factor boundary
    orient PC1 so higher means better sustained attention/performance

  final norming:
    residualize raw PC1 on sex_c + age_at_test + age_at_test^2
    z-score the residual to write gradcpt_perf_z_age_sex

Flanker:
  valid sitting:
    flag_accuracy == 0
    flag_trial_flags == 0
    test_restarted is false when available

  broader candidate:
    log(rcs_incongruent + eps)
    log(rcs_congruent + eps)
    -accuracy_interference
    -median_rt_interference

  predeclared split-efficiency candidate:
    log(rcs_incongruent + eps)
    log(rcs_congruent + eps)
    unit mean of the two transformed, winsorized, z-scored indicators

  selected source:
    official AoU outcomes.score, named flanker_score in the code
    the score is "Overall rate correct score transformed to a 0-100 scale"
    this official score correlated 0.979 with the split-efficiency composite
    so it won by the simple-score >=0.95 rule

  final norming:
    residualize official Flanker score on sex_c + age_at_test + age_at_test^2
    z-score the residual to write flanker_efficiency_z_age_sex
```

Current aggregate diagnostics from the scored cohort are below. These are run
outputs, not pipeline constants; they are regenerated into the diagnostic TSVs
whenever the command is run on a new cohort/CDR.

| Recommended score | Source used by the pipeline | Scored N | Pearson r with `ses_ea_proxy_z` | Spearman r |
|---|---|---:|---:|---:|
| `dd_patience_z_age_sex` | official `-lnk` | 15,148 | 0.281 | 0.266 |
| `gradcpt_perf_z_age_sex` | PC1 of `dprime`, `-log(cv_rtc)`, `-log(median_rtc)` | 15,195 | 0.260 | 0.247 |
| `flanker_efficiency_z_age_sex` | selected Flanker efficiency/simple score | 14,629 | 0.250 | 0.240 |
| `emorecog_perf_z_age_sex` | PC1 of `score`, `-log(cv_rtc)`, `-log(median_rtc)` | 18,486 | 0.116 | 0.108 |

The main PC1 loadings used in the current run were:

| Score | Indicator | PC1 loading |
|---|---|---:|
| GradCPT | `dprime` | 0.661 |
| GradCPT | `-log(cv_rtc)` | 0.614 |
| GradCPT | `-log(median_rtc)` | 0.432 |
| Emotional Recognition | `score` | 0.649 |
| Emotional Recognition | `-log(cv_rtc)` | -0.423 |
| Emotional Recognition | `-log(median_rtc)` | 0.633 |

The negative Emotional Recognition `cv_rtc` loading is kept intentionally in
the recommended PC1 because flipping a single indicator after PCA would define a
different construct. It means this empirical score/RT PC is mostly accuracy plus
speed, with RT variability entering in the opposite direction for this task.

In the cdrv9 final EUR teacher-z cohort (`N = 280,101`), valid ETM sittings
were counted after applying the same task-specific QC filters used by the
scorer. Counts from 1 to 20 are suppressed under the AoU reporting rule; exact
zero counts are shown.

| Valid sittings | Delay Discounting | GradCPT | Flanker | Emotional Recognition |
|---:|---:|---:|---:|---:|
| 0 | 232,290 (82.931%) | 231,330 (82.588%) | 234,624 (83.764%) | 223,100 (79.650%) |
| 1 | 37,413 (13.357%) | 38,150 (13.620%) | 35,745 (12.761%) | 42,553 (15.192%) |
| 2 | 8,832 (3.153%) | 9,104 (3.250%) | 8,265 (2.951%) | 11,937 (4.262%) |
| 3 | 1,208 (0.431%) | 1,200 (0.428%) | 1,121 (0.400%) | 1,925 (0.687%) |
| 4 | 226 (0.081%) | 199 (0.071%) | 217 (0.077%) | 381 (0.136%) |
| 5 | 71 (0.025%) | 58 (0.021%) | 66 (0.024%) | 123 (0.044%) |
| 6 | 25 (0.009%) | 23 (0.008%) | 33 (0.012%) | 38 (0.014%) |
| 7 | suppressed | suppressed | suppressed | suppressed |
| 8 | suppressed | suppressed | suppressed | suppressed |
| 9 | suppressed | suppressed | suppressed | suppressed |
| 10 | suppressed | suppressed | suppressed | suppressed |
| 11 | suppressed | suppressed | suppressed | suppressed |
| 12 | suppressed | suppressed | 0 (0.000%) | suppressed |
| 13 | suppressed | suppressed | suppressed | suppressed |
| 14 | suppressed | suppressed | suppressed | suppressed |

Cumulative valid-sitting coverage in the same cohort was:

| Threshold | Delay Discounting | GradCPT | Flanker | Emotional Recognition |
|---|---:|---:|---:|---:|
| At least 1 | 47,811 (17.069%) | 48,771 (17.412%) | 45,477 (16.236%) | 57,001 (20.350%) |
| At least 2 | 10,398 (3.712%) | 10,621 (3.792%) | 9,732 (3.474%) | 14,448 (5.158%) |
| At least 3 | 1,566 (0.559%) | 1,517 (0.542%) | 1,467 (0.524%) | 2,511 (0.896%) |

Short-interval repeat valid sittings are uncommon, but they provide a useful
scoring sanity check. This is mostly a same-day or short-gap repeatability
diagnostic, not the longer-interval reliability estimate used below for the
latent cognitive ability lower-bound calculation. The diagnostic fits the
production scoring recipe on first valid sittings, applies those same
transforms/loadings/age-sex residualization parameters unchanged to each
person's second valid sitting, and then correlates first-vs-second task scores.

| Task score | Repeat pairs | Pearson r | Spearman r | Mean second - first z |
|---|---:|---:|---:|---:|
| `dd_patience_z_age_sex` | 223 | 0.781 | 0.771 | 0.036 |
| `gradcpt_perf_z_age_sex` | 64 | 0.892 | 0.904 | 0.150 |
| `flanker_efficiency_z_age_sex` | 173 | 0.776 | 0.772 | 0.121 |
| `emorecog_perf_z_age_sex` | 209 | 0.802 | 0.758 | 0.070 |

The repeat gaps are mostly same-day retries rather than year-scale retests:
65.0% for DD, 84.4% for GradCPT, 57.2% for Flanker, and 69.4% for Emotional
Recognition. These correlations support short-interval score repeatability,
but they should not be used as long-term task reliability estimates.

As a longer-interval cdrv9 check, we separately computed test-retest
correlations in the full classified-European ancestry set rather than only the
SES-EA proxy cohort. The sample universe was the pipeline's European keep-list
(`303,903` IIDs). The query used `C2025Q4R6` ETM tables and the same valid-
sitting QC filters as the scorer:

```text
GradCPT valid sitting:
  flag_trial_flags == 0
  flag_non_response == 0
  flag_omission_error_rate == 0
  test_restarted == false
  dprime is finite

Flanker valid sitting:
  flag_accuracy == 0
  flag_trial_flags == 0
  test_restarted == false
  official score is finite
```

For the main retest definition, each participant's first score is their
earliest valid sitting and the retest score is their first later valid sitting
more than 30 days after the first. Scores are the official primary task
outcomes, not the age/sex-normalized production proxy scores.

| Task | Score | Valid sittings | People with a valid sitting | People with 2+ valid sittings | Retest pairs >30d | Pearson r | Spearman r | Mean retest - first | Median gap days | Min gap days | Max gap days |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GradCPT | `dprime` | 64,389 | 51,197 | 10,999 | 10,399 | 0.6925 | 0.6888 | 0.0012 | 305.3 | 30.006 | 440.0 |
| Flanker | official `score` | 60,064 | 47,784 | 10,099 | 9,541 | 0.7241 | 0.7251 | 0.4202 | 298.9 | 30.033 | 644.9 |

A stricter sensitivity restricted to participants with exactly two valid
sittings, again more than 30 days apart, was nearly identical:

| Task | Score | Retest pairs >30d | Pearson r | Spearman r | Mean retest - first | Median gap days | Min gap days | Max gap days |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| GradCPT | `dprime` | 8,843 | 0.6889 | 0.6853 | -0.0042 | 329.0 | 30.006 | 440.0 |
| Flanker | official `score` | 8,029 | 0.7230 | 0.7230 | 0.3397 | 328.9 | 30.042 | 644.9 |

These longer-interval retest results are lower than the same-day-heavy
production-score repeat diagnostic above, as expected, but they are based on
much larger samples and are a better estimate of year-scale task stability.

Individual-level cognitive score files stay in:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/
```

Aggregate diagnostics can be staged to:

```text
regenie_input/ses_ea_proxy/scrap/etm_cog_task_factors/
```

Diagnostics include factor loadings, indicator correlation matrices, dropped
redundant indicators, factor-vs-simple-score correlations, SES-EA/teacher
correlations, repeat-sitting counts, and age/sex plus administration-metadata
sensitivity checks.

The most useful output files are:

```text
etm_cog_task_factors_recommended_wide.tsv
  One row per SES-EA proxy-cohort sample, with the four recommended
  age/sex-normalized cognitive scores and task-specific sitting metadata.

etm_cog_task_factors_recommended_sources.tsv
  The exact source score used for each recommended score.

etm_cog_task_factors_correlations.tsv
  Proxy, teacher, EA-years, and proxy-plus-teacher correlations for all
  candidate and recommended scores, split by combined/oof/applied role.

etm_cog_task_factors_recommended_cross_task_correlations.tsv
  Pairwise correlations among the recommended task scores and the existing
  three-domain ETM-g score when that output is present.

etm_cog_task_factors_emorecog_*.tsv
  Emotional Recognition-specific indicator, loading, simple-score, QC,
  timeout, age/sex, administration-metadata, and cross-task diagnostics.
```

### ETM general cognitive/performance factor scoring

`run_etm_g_from_task_scores.sh` builds downstream ETM general
cognitive/performance factor scores from the recommended task scores. It reads
the task-score output and the SES-EA proxy cohort files only. It does not query
ETM tables again and does not run GWAS.

```bash
bash run_etm_g_from_task_scores.sh --stage-aggregate
```

The primary three-domain score is retained for continuity:

```text
dd_patience_z_age_sex          # official -lnk, age/sex-normalized upstream
gradcpt_perf_z_age_sex         # GradCPT PC1 score, age/sex-normalized upstream
flanker_efficiency_z_age_sex   # Flanker efficiency score, age/sex-normalized upstream
```

When `emorecog_perf_z_age_sex` is present, the same command also writes a
four-domain score:

```text
etm_g4_z                       # DD + GradCPT + Flanker + Emotional Recognition
```

The ETM-g command fits a one-factor Gaussian factor model only in participants
with all selected task scores observed. Before model fitting, the selected task
scores are re-centered and re-scaled using the all-task complete-case reference
sample:

```text
x_ij = (task_ij - mean_j_complete_case) / sd_j_complete_case

x_DD       = lambda_DD       * g + error_DD
x_GradCPT  = lambda_GradCPT  * g + error_GradCPT
x_Flanker  = lambda_Flanker  * g + error_Flanker
x_EmoRecog = lambda_EmoRecog * g + error_EmoRecog   # four-domain model only
```

The learned loadings and uniquenesses are then applied to everyone with at
least one observed selected task score. Missing tasks are not imputed. For each
observed task pattern, the score is the regression factor-score point estimate:

```text
g_hat_i = lambda_O' * Sigma_O^{-1} * x_iO
Sigma_O = lambda_O lambda_O' + diag(psi_O)
```

One-task scores are intentionally shrunk toward zero. The command does not
divide by the sum of available weights and does not z-score within missingness
pattern. The final ETM-g score is z-scored using the complete-case `g_hat`
distribution. It is not residualized again for age or sex because the task
inputs were already age/sex-normalized upstream. Age/sex correlations and
regressions are written as validation diagnostics only.

Current complete-case FA loadings and external-validity diagnostics:

| Model | Complete-case N | Scored N | DD loading | GradCPT loading | Flanker loading | EmoRecog loading | Pearson r with proxy | Pearson r with teacher |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Three-domain `etm_g_z` | 11,165 | 19,613 | 0.170 | 0.658 | 0.706 | NA | 0.290 | 0.202 |
| Four-domain `etm_g4_z` | 11,077 | 22,121 | 0.154 | 0.789 | 0.588 | 0.321 | 0.274 | 0.188 |

Among complete cases only, the proxy correlations were:

| Model | Complete-case definition | N | Pearson r with proxy | Spearman r with proxy |
|---|---|---:|---:|---:|
| `etm_g_z` | DD + GradCPT + Flanker observed | 11,165 | 0.306 | 0.291 |
| `etm_g4_z` | DD + GradCPT + Flanker + EmoRecog observed | 11,077 | 0.298 | 0.282 |

The four-domain score passes the loading checks and adds coverage for people
with only Emotional Recognition observed. In the current run, however, it is
slightly weaker against the SES-EA proxy than the three-domain score. The
three-domain score is therefore kept as the continuity/primary ETM-g phenotype,
with `etm_g4_z` available as the four-domain coverage/sensitivity score.

The current four-task positive-manifold check was:

| Pair | Pair-available N | Pearson r |
|---|---:|---:|
| DD vs GradCPT | 12,487 | 0.114 |
| DD vs Flanker | 12,271 | 0.128 |
| DD vs Emotional Recognition | 13,887 | 0.057 |
| GradCPT vs Flanker | 11,766 | 0.466 |
| GradCPT vs Emotional Recognition | 13,459 | 0.270 |
| Flanker vs Emotional Recognition | 13,148 | 0.182 |

The age/sex validation checks on the final ETM-g scores were small:

| Score | N | Pearson r with `sex_c` | Pearson r with `yob_c` | Linear age/sex R^2 | Quadratic age/sex R^2 |
|---|---:|---:|---:|---:|---:|
| `etm_g_z` | 19,613 | -0.002 | 0.002 | 0.000006 | 0.000011 |
| `etm_g4_z` | 22,121 | -0.001 | 0.003 | 0.000011 | 0.000011 |

Local individual-level outputs stay in:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_general_factor/
```

The main files are:

```text
etm_general_factor_scores_wide.tsv
  One row per SES-EA proxy-cohort sample. `etm_g_z` is the three-domain score.
  `etm_g4_z` is present when the four-domain model passes the loading/positive-
  manifold checks and the person has at least one observed task score among the
  four domains.

etm_general_factor_scores_scored_only.tsv
  Convenience subset with only participants who have at least one observed ETM
  task score.
```

Aggregate diagnostics are written locally under:

```text
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_general_factor/diagnostics/
```

With `--stage-aggregate`, only diagnostic tables are staged to:

```text
regenie_input/ses_ea_proxy/scrap/etm_cog_task_factors/etm_general_factor/
```

Diagnostics include complete-case and pair-available task correlations, FA
loadings and uniquenesses, model-implied and residual correlations,
pattern-specific scoring weights, score distributions by missingness pattern,
age/sex validation, SES-EA/teacher/EA-years correlations, comparison scores,
the three-vs-four ETM-g comparison, and the GradCPT/Flanker
`etm_attention_exec_z` diagnostic score. These diagnostics do not choose or
alter the phenotype.

### G4-finetuned SES-EA proxy scoring

`run_g4_finetuned_ea_proxy.sh` is an optional downstream command that starts
from the six saved SES-EA proxy XGBoost boosters and continues training them
toward the four-domain ETM general cognitive/performance target, `etm_g4_z`.
It does not run GWAS, re-query BigQuery, retrain the original SES-EA proxy from
scratch, or overwrite the base booster JSON files.

```bash
bash run_g4_finetuned_ea_proxy.sh --stage-aggregate
```

The output score is:

```text
g4_finetuned_ea_proxy_z
```

This is a G4-targeted survey proxy initialized from the SES-EA proxy booster.
It should be interpreted as a cognitive-function-enriched predictive score,
not as an observed cognitive score or a causal mediation estimate.

The command consumes:

```text
regenie_input/ses_ea_proxy/all_scores.tsv
regenie_input/ses_ea_proxy/xgboost_feature_columns.json
regenie_input/ses_ea_proxy/xgboost_model_manifest.tsv
regenie_input/ses_ea_proxy/xgboost_models/fold_0.json ... fold_4.json
regenie_input/ses_ea_proxy/xgboost_models/final_model.json
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_cog_task_factors_recommended_wide.tsv
data/regenie/ses_ea_proxy_scrap/etm_cog_task_factors/etm_general_factor/etm_general_factor_scores_wide.tsv
```

It rebuilds the SES-EA survey feature matrix by reusing the same feature helper
code as `setup_ses_ea_proxy_gwas.py`, then requires the ordered feature list to
match `xgboost_feature_columns.json` exactly. The SHA-256 hash must match both
`xgboost_model_manifest.tsv` and the saved booster attributes. ETM task scores
and ETM-g values are labels and diagnostics only; they are never used as model
features.

The primary fine-tuning labels are participants with finite `etm_g4_z` who also
completed at least one of the two strongest ETM domains:

```text
finite(etm_g4_z)
AND (
  finite(gradcpt_perf_z_age_sex)
  OR finite(flanker_efficiency_z_age_sex)
)
```

Samples with ETM-g4 but neither GradCPT nor Flanker are retained for prediction
and diagnostics, but not used as primary fine-tuning labels. A sensitivity mode
can use all finite ETM-g4 labels:

```bash
bash run_g4_finetuned_ea_proxy.sh --target all-etm-g4 --stage-aggregate
```

The fine-tuning keeps the original cross-fit design. For OOF fold `k`, the
command loads `xgboost_models/fold_k.json`, fine-tunes only on target-eligible
`role=oof` samples with `fold_id != k`, and predicts the original held-out
`fold_id == k` samples. The final applied model loads `final_model.json`,
fine-tunes on all target-eligible OOF samples, and predicts only
`role=applied` samples. Thus, target-labeled samples are predicted by models
that did not fine-tune on those samples.

Continued XGBoost training appends trees to boosters whose raw outputs are on
the original SES-EA proxy raw-score scale. To keep the appended trees on that
scale, the ETM-g4 label is aligned within each fold's allowed fine-tuning rows
before training:

```text
y_g4_aligned =
  mean(base_score_raw_train)
  + sd(base_score_raw_train)
    * ((etm_g4_z - mean(etm_g4_z_train)) / sd(etm_g4_z_train))
```

The model trains on `y_g4_aligned`, but validation diagnostics still correlate
against the original `etm_g4_z`. The final phenotype is then standardized over
the full SES-EA proxy cohort:

```text
g4_finetuned_ea_proxy_z =
  zscore(g4_finetuned_ea_proxy_raw over all scored cohort rows)
```

Fine-tuned outputs are written locally under:

```text
data/regenie/ses_ea_proxy_scrap/g4_finetuned_ea_proxy/
```

and copied to the workspace bucket under:

```text
regenie_input/g4_finetuned_ea_proxy/
```

The main files are:

```text
g4_finetuned_ea_proxy_scores_wide.tsv
g4_finetuned_ea_proxy_oof_scores.tsv
g4_finetuned_ea_proxy_applied_scores.tsv
phen.txt
base_covar.txt
covar.txt
training_iids.txt
g4_finetuned_model_manifest.tsv
g4_finetuned_params.tsv
g4_finetuned_runtime_manifest.json
xgboost_models/fold_0_g4_finetuned.json
xgboost_models/fold_1_g4_finetuned.json
xgboost_models/fold_2_g4_finetuned.json
xgboost_models/fold_3_g4_finetuned.json
xgboost_models/fold_4_g4_finetuned.json
xgboost_models/final_model_g4_finetuned.json
```

With `--stage-aggregate`, aggregate diagnostic tables are also staged to:

```text
regenie_input/ses_ea_proxy/scrap/g4_finetuned_ea_proxy/
```

Those diagnostics include target-overlap counts, label-alignment means and SDs,
selected appended-tree rounds, before/after correlations with `etm_g4_z`,
teacher and EA-years retention, covariate correlations with YOB/sex/PCs, score
distributions in ETM and non-ETM samples, feature-importance comparisons, and
task-pattern validation. Review these diagnostics before using
`g4_finetuned_ea_proxy_z` as a downstream phenotype.

### GradCPT/Flanker fine-tuned SES-EA proxy scoring

`run_gradcpt_flanker_finetuned_ea_proxy.sh` uses the same continued-training
machinery as the G4 fine-tune, but it targets only the two strongest ETM task
domains:

```bash
bash run_gradcpt_flanker_finetuned_ea_proxy.sh --stage-aggregate
```

The fine-tuning label is defined only for samples with both task scores:

```text
gradcpt_flanker_mean_z = zscore(mean(gradcpt_perf_z_age_sex,
                                     flanker_efficiency_z_age_sex))
```

For each fold, the command excludes the held-out fold from fine-tuning, aligns
the GradCPT and Flanker component scores to that fold booster's raw prediction
mean and SD, averages the aligned components, and z-matches that combined label
to the same booster raw scale before appending trees. Validation diagnostics
are reported against the unaligned `gradcpt_flanker_mean_z`, plus `teacher_z`
and EA years.

The first-stage fine-tuned booster score is:

```text
gradcpt_flanker_finetuned_ea_proxy_z
```

The output needed by the final selected phenotype is the first-stage
fine-tuned booster score:

```text
gradcpt_flanker_finetuned_ea_proxy_z
```

Earlier runs also wrote a teacher-including linear calibration as a diagnostic
branch. That branch is retained in some output tables for audit/history, but it
is not the final selected phenotype. The final selected phenotype is built later
by `run_gradcpt_flanker_factor18_no_teacher_calibrated_proxy_gwas.sh`, using
this fine-tuned booster score together with the original SES-EA proxy and the
direct scratch-XGBoost score.

Outputs are written locally under:

```text
data/regenie/ses_ea_proxy_scrap/gradcpt_flanker_finetuned_ea_proxy/
```

and copied to:

```text
regenie_input/gradcpt_flanker_finetuned_ea_proxy/
```

If `SES_EA_PROXY_GWAS_INPUT_NAME` is not the legacy `ses_ea_proxy`, the wrapper
automatically suffixes the output name to avoid overwriting earlier runs. For
example, rebuilding from `ses_ea_proxy_v2_kinholdout` writes to:

```text
data/regenie/ses_ea_proxy_scrap/gradcpt_flanker_finetuned_ea_proxy_ses_ea_proxy_v2_kinholdout/
regenie_input/gradcpt_flanker_finetuned_ea_proxy_ses_ea_proxy_v2_kinholdout/
regenie_input/ses_ea_proxy_v2_kinholdout/scrap/gradcpt_flanker_finetuned_ea_proxy_ses_ea_proxy_v2_kinholdout/
```

The selected fine-tuning defaults for this two-domain target are:

```text
eta = 0.03
max_depth = 4
min_child_weight = 20
lambda = 2
max_rounds = 1000
early_stopping_rounds = 50
valid_fraction = 0.20
```

These settings append new trees to the saved SES-EA boosters; the original
booster trees are reused as a fixed initialization and are not modified. The
wrapper sets these defaults, but the `G4_FINETUNE_*` environment variables can
still override them for sensitivity runs.

During development, this target mode also wrote a teacher-including
second-stage linear calibration as a diagnostic branch:

```text
gradcpt_flanker_mean_z
  ~ teacher_z
  + ses_ea_proxy_z
  + gradcpt_flanker_finetuned_ea_proxy_z
```

That branch was useful for checking how much direct teacher-label information
changed GradCPT/Flanker prediction, but it is not the selected phenotype. The
selected phenotype uses the fine-tuned booster score as one of three inputs in
the final no-teacher calibration.

### Direct GradCPT/Flanker scratch XGBoost proxy

`run_gradcpt_flanker_direct_xgb_proxy.sh` is a downstream benchmark/candidate
that trains survey-to-cognition boosters from scratch rather than appending
trees to the SES-EA proxy boosters:

```bash
bash run_gradcpt_flanker_direct_xgb_proxy.sh --stage-aggregate
```

The command uses the same outer fold structure as the SES-EA proxy. OOF fold
`k` is trained on fit-PCA samples from folds other than `k` and predicts fold
`k`. The final applied model is trained only on fit-PCA samples with
`final_model_train_allowed == 1`, which excludes fit-PCA samples related to the
applied cohort at KING kinship `>= 0.0441941`, then predicts all applied
samples.

The scratch XGBoost model uses the original SES-EA proxy training procedure and
hyperparameters:

```text
eta = 0.05
max_depth = 6
min_child_weight = 20
lambda = 1
subsample = 0.8
colsample_bytree = 0.8
num_boost_round = 2000
early_stopping_rounds = 50
internal validation = 4-fold xgb.cv within the allowed training pool
```

The feature matrix starts with the exact SES-EA proxy feature contract, then
adds The Basics education item `1585940` as both a revised numeric years feature
and one-hot response indicators. This is intentional for this direct cognitive
proxy benchmark: unlike the original EA proxy, this model is not trained to
predict the education item itself.

The training target uses all samples with at least one of the two strongest ETM
task scores, rather than only participants who completed both:

```text
gradcpt_perf_z_age_sex
flanker_efficiency_z_age_sex
```

Among complete cases, the current GradCPT-Flanker correlation is about `0.466`,
so the equal-loading two-indicator model uses loading `sqrt(0.466) ~= 0.683`.
People with both tests receive the two-test regression factor score. People
with only one test receive a shrunken one-test estimate, approximately
`0.683 * observed_task_z`, instead of treating a single task as a full-strength
two-domain cognitive score. The resulting `gradcpt_flanker_factor_z` is
z-scored using the both-test complete-case reference.

For each scratch XGBoost fit, the training labels are rescaled to match the
mean and SD of `ses_ea_proxy_z` in that model's allowed target-labeled training
samples. This keeps the direct target on the same broad scale as the proxy
scores while preserving fold safety.

The direct score needed by the final selected phenotype is:

```text
gradcpt_flanker_direct_xgb_proxy_z
```

Earlier scratch-XGBoost runs also evaluated teacher-including three- and
four-variable calibrations as diagnostics. Those comparisons are useful history,
but the final selected phenotype uses the no-teacher calibration described
below.

## Final Selected Phenotype: No-Teacher GradCPT/Flanker Calibration

The final selected phenotype, and the phenotype used for the cdrv9 GWAS
reported here, is:

```text
gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z
```

It is built from three fold-safe survey-derived predictors:

```text
ses_ea_proxy_z
gradcpt_flanker_finetuned_ea_proxy_z
gradcpt_flanker_direct_xgb_proxy_z
```

The calibration target is `gradcpt_flanker_factor_z`, a missing-pattern-aware
GradCPT/Flanker factor target for everyone with at least one of the two task
scores. Participants with both tests contribute the two-test score; one-test
participants contribute a shrunken one-test factor estimate based on the
GradCPT-Flanker complete-case correlation. The final linear calibration excludes
`teacher_z` by design, so the selected phenotype is not a direct linear blend
with the education teacher label.

The calibration is fold-safe. For OOF fold `k`, the linear model is fit on
target-labeled OOF samples from the other four folds and predicts all samples in
fold `k`. The applied model is fit on target-labeled OOF samples allowed by the
kinship holdout and predicts the applied cohort. In the cdrv9 run, the resulting
raw prediction is z-scored over all 280,101 SES-EA proxy cohort rows.

### Fold-safe prediction and kinship safety

The final phenotype is designed so that each participant's score is produced by
models that did not train on that participant's own row. The participant's own
survey responses, area-SES features, genetic sex feature, and education-response
feature where applicable are used as predictors when scoring that participant.
Those same feature values, and that participant's own teacher/cognitive outcome
labels, are not part of the training data for the model that generates that
participant's prediction.

This rule is applied at every predictive layer:

```text
Original SES-EA proxy:
  OOF fold k is predicted by an XGBoost model trained on the other four
  fit_pca folds. The applied cohort is predicted by a sixth model trained on
  fit_pca samples after excluding fit-PCA relatives of applied samples at
  KING kinship >= 0.0441941.

GradCPT/Flanker fine-tuned EA proxy:
  OOF fold k is predicted by the corresponding fine-tuned fold model. Its
  cognitive fine-tuning labels come only from the other four OOF folds. The
  applied model is fine-tuned only on final_model_train_allowed OOF samples.

Direct GradCPT/Flanker scratch XGBoost proxy:
  OOF fold k is predicted by a scratch XGBoost model trained on GradCPT/Flanker
  labels from the other four OOF folds. The applied model is trained only on
  final_model_train_allowed OOF samples.

Final 3-variable linear calibration:
  OOF fold k is predicted by a linear model fit on target-labeled OOF samples
  from the other four folds. The applied cohort is predicted by a linear model
  fit only on target-labeled final_model_train_allowed OOF samples.
```

The OOF rows are members of the PCA-fit sample set, which was selected using
the third-degree KING relatedness cutoff. Thus, for OOF predictions, the other
fit-PCA folds should not contain relatives of the held-out sample at
`KINSHIP >= 0.0441941`, the cutoff used here for first-cousin/third-degree-or-
closer relatedness. For the applied rows, the `final_model_train_allowed` flag
additionally removes fit-PCA samples related to any applied sample at the same
`0.0441941` threshold before fitting the sixth/final models. In practical terms,
each final score comes from survey-based predictions made by models trained on
held-out and threshold-based kinship-clean training samples.

The selected run produced a finite phenotype for every row in the SES-EA proxy
cohort and used the same sample count as the final GWAS input:

```text
Total rows:                               280,101
OOF / fit_pca rows:                       252,774
Applied rows:                              27,327
Calibration target labels, either task:    55,528
OOF target labels:                         49,492
Applied target labels:                      6,036
Final-model kinholdout target labels:      46,744
```

Primary validation correlations were:

| Group | Target | N | Pearson r | Spearman r |
|---|---|---:|---:|---:|
| Full cohort | `teacher_z` | 280,101 | 0.641148 | 0.623209 |
| OOF | `teacher_z` | 252,774 | 0.640220 | 0.621428 |
| Applied | `teacher_z` | 27,327 | 0.646964 | 0.630312 |
| Either GradCPT or Flanker | `gradcpt_flanker_factor_z` | 55,528 | 0.393148 | 0.375928 |
| Both GradCPT and Flanker | `gradcpt_flanker_factor_z` | 38,720 | 0.403959 | 0.384074 |
| Both GradCPT and Flanker | `gradcpt_flanker_mean_z` | 38,720 | 0.404018 | 0.384114 |

Predictor-level GradCPT/Flanker validation correlations were:

| Group | Predictor | Target | N | Pearson r | Spearman r |
|---|---|---|---:|---:|---:|
| Applied both-task target | `ses_ea_proxy_z` | `gradcpt_flanker_mean_z` | 4,237 | 0.332125 | 0.303597 |
| Applied both-task target | `gradcpt_flanker_finetuned_ea_proxy_z` | `gradcpt_flanker_mean_z` | 4,237 | 0.376375 | 0.351244 |
| Applied both-task target | `gradcpt_flanker_direct_xgb_proxy_z` | `gradcpt_flanker_mean_z` | 4,237 | 0.416083 | 0.381262 |
| Combined both-task target | `ses_ea_proxy_z` | `gradcpt_flanker_mean_z` | 38,720 | 0.313799 | 0.293417 |
| Combined both-task target | `gradcpt_flanker_finetuned_ea_proxy_z` | `gradcpt_flanker_mean_z` | 38,720 | 0.359032 | 0.340223 |
| Combined both-task target | `gradcpt_flanker_direct_xgb_proxy_z` | `gradcpt_flanker_mean_z` | 38,720 | 0.397867 | 0.378689 |
| Combined either-task target | `ses_ea_proxy_z` | `gradcpt_flanker_factor_z` | 55,528 | 0.304134 | 0.285639 |
| Combined either-task target | `gradcpt_flanker_finetuned_ea_proxy_z` | `gradcpt_flanker_factor_z` | 55,528 | 0.350399 | 0.334417 |
| Combined either-task target | `gradcpt_flanker_direct_xgb_proxy_z` | `gradcpt_flanker_factor_z` | 55,528 | 0.387313 | 0.370425 |
| Final no-teacher phenotype | `gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z` | `gradcpt_flanker_factor_z` | 55,528 | 0.393148 | 0.375928 |
| Final no-teacher phenotype | `gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z` | `gradcpt_flanker_mean_z` | 38,720 | 0.404018 | 0.384114 |

Because the complete-case GradCPT/Flanker cohort is not education-response
balanced relative to the full teacher/EA proxy cohort, we also checked the
same correlations after matching the complete-case cohort back to the full
cohort's education-response distribution. The reference distribution was all
280,101 rows with `teacher_z` and an EA response. The comparison cohort was the
set of rows with finite `gradcpt_flanker_mean_z`. Sparse EA-response categories
that did not pass the `>20` reporting rule in the complete-case cohort were
suppressed and excluded from the stratified comparison. The final resampling
target size after that rule was 38,718, with fixed seed `20260702`.

Education-response distributions before and after matching were:

| EA response | Full % | GradCPT/Flanker % | Matched % | Action |
|---|---:|---:|---:|---|
| Advanced Degree | 31.970 | 39.873 | 31.970 | undersample |
| College Graduate | 28.433 | 30.926 | 28.434 | undersample |
| College One to Three | 25.539 | 22.338 | 25.539 | oversample |
| Twelve Or GED | 11.824 | 6.408 | 11.824 | oversample |
| Nine Through Eleven | 1.796 | 0.375 | 1.795 | oversample |
| Five Through Eight | 0.439 | 0.080 | 0.439 | oversample |

Correlations with `gradcpt_flanker_mean_z` increased after matching:

| Predictor | Original r | EA-matched resample r |
|---|---:|---:|
| `teacher_z` baseline | 0.2181 | 0.2529 |
| SES-EA proxy | 0.3138 | 0.3494 |
| Fine-tuned GradCPT/Flanker proxy | 0.3590 | 0.3887 |
| Direct XGBoost GradCPT/Flanker proxy | 0.3979 | 0.4214 |
| Final no-teacher GWAS phenotype | 0.4040 | 0.4289 |

A deterministic inverse-stratum-weighted Pearson check gave similar results:

| Predictor | Weighted Pearson r |
|---|---:|
| `teacher_z` baseline | 0.2536 |
| SES-EA proxy | 0.3464 |
| Fine-tuned GradCPT/Flanker proxy | 0.3889 |
| Direct XGBoost GradCPT/Flanker proxy | 0.4232 |
| Final no-teacher GWAS phenotype | 0.4302 |

Interpretation: the complete-case GradCPT/Flanker cohort is enriched for
higher-education responses and underrepresents lower-education responses.
Matching it back to the full EA-response distribution restores more
between-education contrast, so all predictor correlations rise. The final
no-teacher phenotype remains the strongest predictor in both the original and
matched comparisons, though its gain over the direct XGBoost proxy is modest.

As a reliability-aware validation, we also estimated conservative lower bounds
on the final phenotype's correlation with latent general cognitive ability.
This is not a point estimate of `r(proxy, g)`. A task retest correlation treats
both general ability and task-specific stable variance as reliable signal, so
dividing the proxy-task correlation by `sqrt(task reliability)` disattenuates
toward each task's stable true score. If the task contains stable task-specific
variance, the result is a lower bound for `r(proxy, g)` under the assumption
that the proxy reaches the task mainly through general ability.

This check used the exact production task z-scores that feed
`gradcpt_flanker_mean_z`: `gradcpt_perf_z_age_sex` for GradCPT and
`flanker_efficiency_z_age_sex` for Flanker. The Flanker column name is kept for
pipeline compatibility, but the cdrv9 diagnostic selected `flanker_simple_score`
as its source because it correlated at least 0.95 with the efficiency split
score. We refit the production scoring recipe on first valid sittings, applied
the same transforms/loadings/age-sex residualization parameters to later valid
sittings, and verified exact agreement with the saved pipeline columns
(`r = 1.0`, max absolute difference below `5e-16`).

Retest reliability was estimated among final GWAS-cohort samples with a first
valid sitting and a first later valid sitting more than 30 days and no more
than four years after the first:

| Task score | Source score | Retest pairs | Pearson reliability | Spearman reliability | Mean retest - first z | Median gap days |
|---|---|---:|---:|---:|---:|---:|
| `gradcpt_perf_z_age_sex` | `gradcpt_perf_factor` | 10,051 | 0.744942 | 0.736733 | 0.044443 | 306.7 |
| `flanker_efficiency_z_age_sex` | `flanker_simple_score` | 9,199 | 0.669197 | 0.675637 | 0.076173 | 300.2 |

A stricter sensitivity restricted to participants with exactly two scoreable
sittings gave similar reliabilities: `0.739571` for GradCPT and `0.663849` for
Flanker.

Among participants with exactly two scoreable sittings in this same
`>30 days` and `<=4 years` retest window, the gap between first and second
sitting was concentrated around the next annual measurement wave:

| Gap percentile | GradCPT days | Flanker days |
|---:|---:|---:|
| 10th | 75.9 | 73.7 |
| 20th | 157.2 | 154.1 |
| 30th | 244.9 | 236.9 |
| 40th | 285.9 | 280.0 |
| 50th | 329.5 | 329.3 |
| 60th | 336.4 | 336.4 |
| 70th | 338.9 | 338.9 |
| 80th | 341.7 | 341.9 |
| 90th | 345.4 | 346.0 |

For the proxy-task correlations, the complete-case task cohorts were weighted
back to the full final-phenotype EA-response distribution using deterministic
inverse stratum weights. EA-response categories that did not pass the `>20`
reporting rule in a task cohort were suppressed and excluded from that task's
weighted comparison (`One Through Four` for GradCPT; `Never Attended` and
`One Through Four` for Flanker).

| Task score | N | Unweighted Pearson r | EA-weighted Pearson r |
|---|---:|---:|---:|
| `gradcpt_perf_z_age_sex` | 48,769 | 0.357281 | 0.383354 |
| `flanker_efficiency_z_age_sex` | 45,474 | 0.343408 | 0.368800 |

The equivalent fixed-seed EA-matched resample check, using one task at a time
rather than requiring both tasks, gave:

| Task score | EA-matched resample r | Reliability | Resample lower bound |
|---|---:|---:|---:|
| `gradcpt_perf_z_age_sex` | 0.383032 | 0.744942 | 0.443787 |
| `flanker_efficiency_z_age_sex` | 0.370720 | 0.669197 | 0.453178 |

Dividing the proxy-task correlations by `sqrt(reliability)` gives these
task-specific lower bounds:

| Task score | `sqrt(reliability)` | Unweighted lower bound | EA-weighted lower bound | Bootstrap 95% CI for EA-weighted bound |
|---|---:|---:|---:|---:|
| `gradcpt_perf_z_age_sex` | 0.863100 | 0.413951 | 0.444160 | 0.433428 to 0.454687 |
| `flanker_efficiency_z_age_sex` | 0.818045 | 0.419791 | 0.450831 | 0.436475 to 0.464468 |

The two task-specific estimates are close. The tighter lower bound is the
Flanker-based estimate, about `0.451` on the EA-weighted scale. Using the
exactly-two-sittings reliability instead gives `0.445769` for GradCPT and
`0.452643` for Flanker; shifting the primary reliability by +/-0.03 keeps the
EA-weighted lower bounds in the same range (`0.435477` to `0.453383` for
GradCPT, `0.441053` to `0.461289` for Flanker).

Taken together, the final no-teacher GWAS phenotype has an estimated
lower-bound correlation of about `0.44` to `0.45` with latent cognitive ability
from these GradCPT/Flanker validation checks, with bootstrap uncertainty
roughly spanning `0.43` to `0.46`.

The local analysis outputs are in:

```text
data/tmp/g_proxy_latent_g_lower_bound_v9/

teacher_baseline_gradcpt_flanker_mean_correlations.tsv
no_teacher_task_specific_resampled_lower_bounds.tsv
production_task_exactly_two_gap_deciles.tsv
```

The calibration-cohort scale check was:

| Variable | N | Mean | SD | Skew |
|---|---:|---:|---:|---:|
| `ses_ea_proxy_z` | 55,528 | 0.338 | 0.815 | -0.705 |
| `gradcpt_flanker_finetuned_ea_proxy_z` | 55,528 | 0.330 | 0.770 | -0.796 |
| `gradcpt_flanker_direct_xgb_proxy_z` | 55,528 | 0.310 | 0.953 | -0.807 |
| `gradcpt_flanker_factor_z` | 55,528 | -0.050 | 0.964 | -0.133 |
| Final no-teacher calibrated phenotype | 55,528 | 0.333 | 0.897 | -0.813 |

The full-cohort phenotype distribution was standardized after prediction:

| Variable | N | Mean | SD | Skew |
|---|---:|---:|---:|---:|
| `gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z` | 280,101 | 0.000 | 1.000 | -0.636 |

The six fold-safe linear calibration coefficients were:

| Fit | Train N | Predict N | Intercept | `ses_ea_proxy_z` | Fine-tuned XGB | Direct XGB |
|---|---:|---:|---:|---:|---:|---:|
| OOF fold 0 | 39,567 | 50,817 | -0.190793 | 0.014060 | 0.138828 | 0.292954 |
| OOF fold 1 | 39,715 | 50,288 | -0.191847 | 0.002897 | 0.141663 | 0.304754 |
| OOF fold 2 | 39,453 | 50,765 | -0.186495 | 0.004757 | 0.136234 | 0.299020 |
| OOF fold 3 | 39,647 | 50,495 | -0.188028 | 0.005466 | 0.133069 | 0.300163 |
| OOF fold 4 | 39,586 | 50,409 | -0.191318 | 0.009536 | 0.134071 | 0.292828 |
| Applied kinholdout | 46,744 | 27,327 | -0.187863 | 0.004780 | 0.140821 | 0.297797 |

The multivariable coefficient significance check showed that the residual
SES-EA term was not significant after including the two GradCPT/Flanker proxy
terms, while both GradCPT/Flanker terms were highly significant in every fit:

| Fit | `ses_ea_proxy_z` p | Fine-tuned XGB p | Direct XGB p |
|---|---:|---:|---:|
| OOF fold 0 | 0.13817493 | 1.3041346e-31 | 8.2972005e-301 |
| OOF fold 1 | 0.76077458 | 3.8469774e-33 | <1e-300 |
| OOF fold 2 | 0.61955678 | 1.8937197e-30 | 1.5232672e-309 |
| OOF fold 3 | 0.56879236 | 4.1985292e-29 | <1e-300 |
| OOF fold 4 | 0.31694220 | 5.9990936e-30 | 3.1707073e-302 |
| Applied kinholdout | 0.58651081 | 2.7751468e-38 | <1e-300 |

This means the final no-teacher linear combiner is driven by the direct XGBoost
proxy and the fine-tuned proxy. The original SES-EA proxy contributes little
additional independent signal once those two survey-derived cognitive proxies
are in the same regression.

The final GWAS is configured to use:

```text
Phenotype:  gradcpt_flanker_factor18_no_teacher_calibrated_proxy_z
Samples:    280,101
Covariates: sex_c + PC1_AVG ... PC10_AVG
RINT:       enabled
Step 1:     494,816 direct variants
Step 2:     7,252,333 WGS variants across 22 chromosomes
```

After completion, lightweight GWAS outputs are written under:

```text
regenie_output/g_ea_proxy_sbayesrc7m_gwas/lightweight/
```

Diagnostic files for the final calibration are written under:

```text
regenie_input/g_ea_proxy_sbayesrc7m/diagnostics/

factor18_no_teacher_calibration_coefficients.tsv
factor18_no_teacher_calibration_correlations.tsv
factor18_no_teacher_calibration_distributions.tsv
```
