# pan_aou_gwas — Pan-UKB-style All of Us HapMap3 residualize-first GWAS

A phenotype-wide GWAS factory over the All of Us survey + physical-measurement
data, run on the pre-built HapMap3 HQ bfile with PLINK2 covariate-free linear
regression on **pre-residualized** phenotypes. The full design is in
[`SPECSHEET.md`](SPECSHEET.md).

**This run produces the GWAS summary-statistic files only.** Downstream LDSC
h²/rg is deferred.

## What is here

```text
SPECSHEET.md                       the finalized specsheet (method + decisions)
scripts/parse_codebooks.py         parse the v9 codebook workbook -> item inventory
scripts/ordinal_rules.py           the ordinal-mapping knowledge base (templates + overrides)
scripts/build_manifests.py         inventory -> the four phenotype manifests
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
| `pfhh_self_allowlist.tsv` | condition | the 33 `self_has_<condition>` PFHH self-history phenotypes |

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
mental-health/substance-use, plus fibromyalgia and recent fractured/broken bones —
run as binary `self_has_<condition>` (case = "Self"; control = completed the
category screen and did not self-report). See SPECSHEET §9.1 and §11.

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
`pca_eur/aou_projected.sscore`, `genetic_sex/sex_covar.txt`, and the
identical-component exclusion list. Set `PAN_AOU_SKIP_MHWB=1` if the off-cycle
Mental Health / Well-Being CDR is unavailable.

## Method (one line)

Covariates (`age_c`, `sex_c`, `age_c:sex_c`, `PC1–PC10`) are regressed out of the
phenotype first (LPM residual for binary, IRNT→residual for quantitative); PLINK2
then runs one covariate-free `--glm allow-no-covars` linear pass over the whole
autosomal bfile — no `--mac`/`--geno`/`--hwe`/`--threads`. This is genetically
identical to covariates-in-PLINK2 and ~30× faster (validated in
`~/projects/ukgwas/covariate_experiment`).

## Local validation

The manifest builders and the phenotype/residualization logic in
`pan_aou_gwas.py` run and were validated locally on a synthetic fixture with
`--skip-gwas`. The BigQuery extraction and PLINK2 association require the AoU
platform; smoke-test them there with `--smoke` before the full run.
```
