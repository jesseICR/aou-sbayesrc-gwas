# pan_aou_gwas — Pan-UKB-style All of Us HapMap3 residualize-first GWAS

A phenotype-wide GWAS factory over the All of Us survey, physical-measurement,
and ZIP3 socioeconomic-context data, run on the pre-built HapMap3 HQ bfile with
PLINK2 covariate-free linear regression on **pre-residualized** phenotypes. The full design is in
[`SPECSHEET.md`](SPECSHEET.md).

**This run produces the GWAS summary-statistic files only.** Downstream LDSC
h²/rg is deferred.

## What is here

```text
SPECSHEET.md                       the finalized specsheet (method + decisions)
scripts/parse_codebooks.py         parse the v9 codebook workbook -> item inventory
scripts/ordinal_rules.py           the ordinal-mapping knowledge base (templates + overrides)
scripts/build_manifests.py         inventory -> the four phenotype manifests
scripts/build_pan_aou_sex_covar.py build the pan-AoU binary sex covariate
scripts/pan_aou_gwas.py            build residualized phenotypes + run PLINK2
run_pan_aou_gwas.sh                on-platform orchestrator (env, extract, keep-list, run)
metadata/                          generated manifests (see below)
```

## Metadata files

| File | One row per | Purpose |
| --- | --- | --- |
| `survey_item_inventory.tsv` | codebook item | audit source: field type, options, validation |
| `survey_question_manifest.tsv` | question | disposition (ordinal/binary/nominal/numeric/excluded/flagged), rule, sensitivity |
| `ordinal_answer_templates.tsv` | (rule, answer phrase) | the shared ordinal rule library, machine-readable |
| `ordinal_mapping_manifest.tsv` | (ordinal question, answer) | the exact answer→value mapping applied at runtime |
| `flagged_questions.tsv` | flagged question | sensitive + medium-confidence + uncertain-ordinal items for review |
| `pfhh_self_allowlist.tsv` | condition | the 33 PFHH self-history conditions (binary + burden sumscore) |
| `external_scores.tsv` | cognitive/proxy score | registry of ETM task + EA-proxy scores to GWAS (verify paths on platform) |
| `sex_specific_items.tsv` | sex-specific item | female/male analysis-sample filters for reproductive, anatomy, and sex-ploidy phenotypes |
| `ea_proxy_feature_sources.tsv` | EA-proxy source question | supplemental live v9 question IDs used by the SES-EA/direct-XGB feature contract but absent from the codebook-derived matcher |
| `composite_items_manifest.tsv` | (instrument, item, answer) | validated composite scores (GAD-7, PHQ-9, PSS, BFI-2, ...) item scoring |
| `COMPOSITE_SCORES.md` | instrument | human-readable definition of each composite (items, scoring, combination) |

Regenerate the manifests from the codebook (safe to run off-platform):

```bash
python3 scripts/parse_codebooks.py
python3 scripts/build_manifests.py
```

## Scope of diagnosis phenotypes (important)

Broad disease / condition diagnoses are **not** run. The Personal Medical
History, Family Health History, and broad Personal-and-Family-Health-History
items (1,174 items) are all excluded. The **only** diagnosis phenotypes are the
33 PFHH **self-history** items in `pfhh_self_allowlist.tsv` — neuro/nervous-system,
mental-health/substance-use, plus fibromyalgia and recent fractured/broken bones.
Each condition yields two phenotypes: a binary `pfhh_self_has_<condition>` (case =
"Self"; control = completed the category screen and did not self-report) **and** a
quantitative `pfhh_burden_<condition>` genetic-relatedness-weighted family-history
sumscore (self = 1, first-degree = 0.5, grandparent = 0.25, summed) as an aggregate
liability proxy. See SPECSHEET §9.1 and §11.

## Run on the All of Us Researcher Workbench

Run from inside an AoU Verily Jupyter terminal, after the main
`get_genotypes.sh` pipeline has produced the European keep-list, genetic sex,
PCs, and the HapMap3 HQ bfile.

```bash
cd aou-sbayesrc-gwas/pan_aou_gwas
bash run_pan_aou_gwas.sh --setup-only   # extract + build residualized phenotypes, no GWAS
bash run_pan_aou_gwas.sh --smoke        # a handful of phenotypes end-to-end (smoke test)
bash run_pan_aou_gwas.sh                 # full run
```

The orchestrator resolves inputs from `${DX_OUTPUT_DIR}` (= workspace-bucket
`sbayesrc_genotypes/`): the HapMap3 bfile, `pca_eur/fit_pca_iids.txt`,
`pca_eur/aou_projected.sscore`, `genetic_sex/sex_covar.txt`,
`genetic_sex/sex_ploidy_qc.tsv`, and the identical-component exclusion list.
The run builds `data/pan_aou_gwas_work/sample_qc/pan_aou_sex_covar.txt` from
the strict main-pipeline sex covariate plus explicit pan-AoU imputation rules
for assigned-male DRAGEN `X0`/`XO` and skipped/prefer-not-to-answer sex-at-birth
rows with DRAGEN `XX`/`XY`. Set `PAN_AOU_SKIP_MHWB=1` if the off-cycle
Mental Health / Well-Being CDR is unavailable.

FYI, the cdrv9/v9 audit rows motivating those explicit pan-AoU additions were:

| Sex at birth response | Pan-AoU sex code | DRAGEN sex ploidy | Inclusion rule | Count |
| --- | ---: | --- | --- | ---: |
| Male | 1 | X0 | assigned-male noncanonical WGS sex ploidy coded male | 1,093 |
| PMI: Skip | missing | XX | missing/nonbinary sex at birth coded from DRAGEN ploidy | 1,071 |
| PMI: Skip | missing | XY | missing/nonbinary sex at birth coded from DRAGEN ploidy | 961 |
| Male | 1 | XO | assigned-male noncanonical WGS sex ploidy coded male | 140 |
| Prefer not to answer | missing | XY | missing/nonbinary sex at birth coded from DRAGEN ploidy | 83 |
| Prefer not to answer | missing | XX | missing/nonbinary sex at birth coded from DRAGEN ploidy | 64 |

The run also extracts a small person-level age covariate and includes
`dragen_x0_xo_male`, a male-only binary GWAS for DRAGEN `X0`/`XO` sex ploidy
versus DRAGEN `XY`, intended as a candidate mosaic loss-of-Y phenotype. This
sex-stratified GWAS is residualized on age and the first 10 genetic PCs.

Female reproductive/anatomy survey phenotypes listed in
`metadata/sex_specific_items.tsv` are run only in pan-AoU female-coded samples.
These sex-stratified GWAS are also residualized on age and the first 10 PCs.

Survey phenotypes use the latest valid response per participant/question, so a
later skip/prefer-not-to-answer response does not mask an earlier valid response.
Two-answer single-select questions emit only one binary GWAS because the two
one-vs-rest encodings are exact complements; the omitted side is listed as
`redundant_binary_complement` in `skipped_phenotypes.tsv`.
Same-survey reused codebook item labels are disambiguated with `_q<question_id>`
in the phenotype ID when the shared label refers to distinct follow-up questions.
PHQ-9 and GAD-7 item and sumscore GWAS pool EHHWB with COPE using EHHWB priority.
PSS-10 item and sumscore GWAS pool SDOH with COPE using SDOH priority. COPE
fill-in is included as a `from_cope` residualization covariate.

The run also extracts `ds_zip_code_socioeconomic` and GWASes the seven numeric
ZIP3 context traits: deprivation index, median income, poverty, assisted income,
no health insurance, vacant housing, and high-school education fraction. Raw
ZIP3 codes and ACS vintage are retained only in the extract for auditability.

## Method (one line)

Covariates (`age_c`, `sex_c`, `age_c:sex_c`, `PC1–PC10`) are regressed out of the
phenotype first (LPM residual for binary, IRNT→residual for quantitative); PLINK2
then runs one covariate-free `--glm allow-no-covars` linear pass over the whole
autosomal bfile — no `--mac`/`--geno`/`--hwe`/`--threads`. This is genetically
identical to covariates-in-PLINK2 and ~30× faster (validated in
`~/projects/ukgwas/covariate_experiment`).

Exceptions are explicit in `phenotype_manifest.tsv`: already age/sex-normalized
external scores use sex+PC covariates, and sex-stratified phenotypes use age+PC
covariates. Pooled PHQ-9/GAD-7 and PSS-10 phenotypes also include `from_cope`.

## Local validation

The manifest builders and the phenotype/residualization logic in
`pan_aou_gwas.py` run and were validated locally on a synthetic fixture with
`--skip-gwas`. The BigQuery extraction and PLINK2 association require the AoU
platform; smoke-test them there with `--smoke` before the full run.
```
