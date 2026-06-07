# AoU SBayesRC GWAS pipeline

A reproducible, idempotent pipeline that builds per-chromosome PLINK2 pfiles
holding the ~7.35 M SBayesRC SNPs from the All of Us v8 ACAF WGS callset, with
the variant `ID` column populated with SBayesRC rsids. It also builds the
~501k direct-SNP bfile, plus a higher-quality direct-SNP bfile filtered with
AoU EUR frequency/missingness metrics. It runs ADMIXTURE K=6 ancestry
projection, compares the resulting ancestry fractions to AoU-provided ancestry
calls/fractions, runs KING kinship from the high-quality direct bfile, compares
those estimates to AoU's provided relatedness table, classifies close
relationships, and selects the unrelated European sample set used to fit PCA.
It then builds the PCA-ready SNP bfile from that sample set, fits 20 European
ancestry PCs, and projects those PCs onto all samples in the high-quality
direct bfile. Before GWAS setup, it builds a conservative sample-QC exclusion
list for identical-genotype components of size three or larger and then builds
final REGENIE Step 1/Step 2 genotype inputs with GWAS-specific variant filters.
The final optional example builds a height phenotype/covariate set and runs a
continuous-trait REGENIE GWAS.

The pipeline **must be run from inside an AoU Verily Jupyter session
terminal** (the standard interactive analysis environment on the All of Us
Researcher Workbench). It cannot be run from an off-platform `wb` CLI
session — VPC SC and IAM together require the orchestrator to run as the
pod's pet service account, which is what the Jupyter pod gives you for free.

From there, the per-chromosome compute runs as **parallel Google Batch jobs
via dsub**, one worker per autosome in `us-central1` (the same region as the
workspace bucket, so output delocalization is intra-region and effectively
line-rate). On a 4-vCPU Jupyter pod a serial run takes ~2–3 days; the dsub
fan-out finishes in ~2 hours.

The pod we developed and tested this on (the minimum-cost default that
worked):

| Component | Value |
|---|---|
| Pod machine type | `n1-highmem-4` |
| vCPUs | 4 |
| RAM | 25 GB |
| Zone | `us-east4-a` (workspace's `terra-default-location`) |
| Local overlay disk | ~492 GB |
| Workspace bucket region | `us-central1` |

The Jupyter pod only orchestrates and submits the dsub tasks — the actual
plink2 work happens on Batch worker VMs we provision in `us-central1`, each
with `--min-cores 4 --min-ram 32 --disk-size 300` by default. You can leave
the Jupyter session open or close the browser (background the run with
`nohup`); the Batch jobs run independently.

## What this pipeline produces

```
${WORKSPACE_BUCKET}/sbayesrc_genotypes/wgs_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/wgs_pfiles/chr{1..22}.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile/chr1_22_merged.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.filter_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.sample_missingness_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_admixture_k6.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/europeans/classified_european_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_vs_ours/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_subset_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_missingness_threshold_counts.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/aou_hq_direct_rel.kin0
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/qc/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/close_relations.csv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/select_pca_europeans.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ready.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_projected.sscore
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.{eigenval,eigenvec,eigenvec.allele}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_eur_counts.acount
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_project_pca.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/sex_covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/genetic_sex_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_components.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_sample_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/gwas_genotype_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/gwas_step1_direct.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/gwas_step2_wgs.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/fit_pca_af/gwas_step2_fit_pca_alt_freqs_passing.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/{phen.txt,covar.txt,training_iids.txt}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/height_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step1/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step2/chr{1..22}/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/lightweight/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/regenie_gwas.summary.tsv

# gwas_dev branch optional command:
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/ea_gwas/{phen.txt,covar.txt,training_iids.txt}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/ea_gwas/ea_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/ea_gwas/ea_answer_counts.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/ea_gwas/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/ea_gwas/lightweight/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/income_gwas/{phen.txt,covar.txt,training_iids.txt}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/income_gwas/income_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/income_gwas/income_answer_counts.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/income_gwas/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/income_gwas/lightweight/
```

Each `summary.tsv` reports `requested / src_variants / src_samples /
biallelic_total / extracted / out_samples / missing / remapped / unmapped`
for that chromosome. A combined summary across all processed chromosomes is
also written to `logs/sbayesrc_extract_summary.tsv` at the end of each run.

## Pipeline steps

### Step 1 — `store_sbayesrc_ids.py` (local)

Reads `data/support/sbayesrc_hg38.csv` (chrom, pos, ref, alt, rsid; downloaded
on first run from a public GitHub release) and writes per-chromosome
`data/sbayesrc_ids/chr{N}.extract.txt` (one `chr{N}:pos:REF:ALT` ID per line)
and `chr{N}.idmap.txt` (`chr{N}:pos:REF:ALT \t rsid`). The ID format matches
the AoU `acaf_threshold/plink_bed/chr{N}.bim` second column exactly, so
`plink2 --extract` is a direct ID match.

### Step 2 — `wgs_extract_variants.sh` (dsub fan-out)

Submits **one dsub task per autosome** to Google Batch in `us-central1`.
Each task runs `dsub_extract_worker.sh` on its own Batch worker, with the
per-chromosome AoU pgen/pvar/psam staged in via `dsub --input` and the final
extracted pgen/pvar/psam/summary.tsv delocalized via `dsub --output-recursive`
to `${WORKSPACE_BUCKET}/sbayesrc_genotypes/wgs_pfiles/`.

Per-worker pipeline (four passes, same as the previous serial-on-Jupyter
implementation):

1. **Pass 1** — `plink2 --pfile acaf_threshold.chr{N} --make-pgen 'multiallelics=-' --output-chr chrM --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10000 --out chr{N}.split`. Splits multi-allelics and assigns `chr:pos:REF:ALT` IDs.

2. **Pass 1.5** (inline awk in `dsub_extract_worker.sh`) — right-trims shared REF/ALT suffix in the split-pvar. plink2's split keeps the longest indel's REF as the anchor for every split row, so SNPs at a multi-allelic site come out over-padded (e.g. `REF=TATG ALT=CATG` for what is really `T:C`). Trimming recovers the minimal form and re-stamps the ID. Purely a text rewrite — pgen indexes by row order with allele codes (0=REF, 1=ALT), so genotypes don't move. On chr22 this raises recovery from 95,130 → 97,960 of 98,065 SBayesRC variants. The canonical, testable reference for this logic is `normalize_pvar_alleles.py`; the worker uses an inline awk version because the default `marketplace.gcr.io/google/ubuntu2204` Batch worker image ships without python3 and VPC SC blocks `apt install`.

3. **Pass 2** — `plink2 --pfile chr{N}.split --extract chr{N}.extract.txt --no-pheno --output-chr chrM --make-pgen --out chr{N}`. Default `pvar-cols` preserves the INFO column (AoU acaf_threshold INFO includes `AC`, `AF`, `AN`, `AS_QUALapprox`, `CALIBRATION_SENSITIVITY`, `QUALapprox`, `SCORE`).

4. **Pass 3** (inline awk) — rewrites the ID column of `chr{N}.pvar` from `chr:pos:REF:ALT` to the SBayesRC rsid via `chr{N}.idmap.txt`.

No QC filters applied — all variants are kept regardless of FILTER/QUAL/
SCORE. ACAF threshold already enforces MAC ≥100 per ancestry. Downstream
code can apply additional filters (e.g. `--extract-if-info
"CALIBRATION_SENSITIVITY < 0.99"`) against the preserved INFO column.

Per-task stdout/stderr lives at
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/logs/dsub/<job-id>.<task>-{stdout,stderr}.log`.

### Step 3 — Direct-SNP bfile for REGENIE step 1

`get_genotypes.sh` downloads the UKBB-derived direct-SNP rsid list to
`data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt` if it is not
already cached. `prepare_direct_snps.py` cross-references that list against
`sbayesrc_hg38.csv` and the extracted WGS pvars, writing per-chromosome rsid
extract lists and a missing-SNP report for direct SNPs absent from the
AoU/SBayesRC pfiles.

`extract_direct_snps.sh` then builds
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_pfiles/chr{N}.{pgen,pvar,psam}`.
Present SNPs are extracted from the WGS pfiles. Direct SNPs absent from the
AoU/SBayesRC pfiles are tracked in `data/direct_snps/missing_direct_snps.tsv`.

`make_direct_bfile.sh` submits one Google Batch/dsub worker in `us-central1`
with an SSD data disk (default: 8 vCPU, 32 GB RAM, 300 GB `pd-ssd`) to merge
the 22 direct pfiles into
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile/chr1_22_merged.{bed,bim,fam}`.
The merged bfile contains the direct SNPs present in the AoU extracted pfiles;
absent SNPs are not encoded as all-missing variants. The merged bfile is
intended for REGENIE step 1.

### Step 4 — Higher-quality direct-SNP bfile

`get_genotypes.sh` downloads `sbayesrc_liftover_results.csv` from the public
SBayesRC liftover release if it is not already cached locally. This file
contains `A1_hg38`, `A2_hg38`, and `A1Freq`, which Step 4 converts to the
SBayesRC hg38 ALT-allele frequency after matching against the AoU direct bfile
REF/ALT alleles.

`make_hq_direct_bfile.sh` builds an AoU EUR keep-list from AoU's computed
ancestry predictions (`ancestry_pred == eur`) and the direct-bfile `.fam`.
It then submits a dsub metrics job that runs plink2 on the raw direct bfile
with that EUR keep-list:

```text
plink2 --bfile chr1_22_merged --keep aou_eur.keep --freq --missing variant-only
```

The local filter builder (`filter_hq_direct_snps.py`) then applies these
filters in order:

```text
1. original UKBB direct-SNP rsid list
2. present in the AoU direct bfile
3. SBayesRC liftover ALT frequency is available and alleles match
4. abs(AoU EUR ALT frequency - SBayesRC ALT frequency) <= 0.04
5. AoU EUR MAF >= 0.007
6. AoU EUR variant missingness <= 0.05
```

The ordered count log is written to:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.filter_summary.tsv
```

The per-variant QC table and final rsid extract list are written alongside it:

```text
chr1_22_merged_hq.variant_qc.tsv
chr1_22_merged_hq.extract.txt
```

Finally, a second dsub job extracts the passing variants from the raw direct
bfile and writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.{bed,bim,fam}
```

The final high-quality bfile keeps the same samples as the raw direct bfile;
only variants are filtered. The worker also computes sample missingness over
the final variant set for all samples and for the AoU EUR keep-list:

```text
chr1_22_merged_hq.sample_missingness_all.smiss
chr1_22_merged_hq.sample_missingness_eur.smiss
chr1_22_merged_hq.sample_missingness_summary.tsv
```

### Step 5 — ADMIXTURE K=6 projection

ADMIXTURE projection starts from the high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

`admixture_prep.sh` downloads two public inputs if they are not already cached
locally under gitignored paths: the ADMIXTURE Linux binary and the K=6 global
reference allele-frequency TSV from `public-statgen`. It stages those files to
the workspace bucket for Batch workers.

The prep worker then:

```text
1. computes all-sample variant missingness on the HQ direct bfile
2. keeps SNPs with all-sample missingness <= 0.05
3. intersects those SNPs with the K=6 reference SNPs
4. extracts the overlap into an ADMIXTURE-specific bfile
5. aligns PLINK BIM alleles to the reference frequencies
6. writes the aligned bfile plus ref_aligned.P
```

The all-sample missingness filter is specific to ADMIXTURE projection. It does
not modify the REGENIE Step-1 HQ direct bfile on disk.

Outputs from prep are written to:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/aou_admixture_aligned.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/ref_aligned.P
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/admixture_prep_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/admixture_align_log.txt
```

`admixture_split_batches.sh` splits the aligned bfile into batches of 20,000
people and copies `ref_aligned.P` to each batch as `batch_NNN.6.P.in`, which is
the filename ADMIXTURE expects in projection mode.

`admixture_run_projection.sh` submits one dsub task per missing batch Q file,
runs:

```text
admixture -j$(nproc) -P batch_NNN.bed 6
```

and concatenates the per-batch Q files into:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_admixture_k6.tsv
```

The final TSV columns are:

```text
FID IID European East_Asian American African South_Asian Oceanian
```

### Step 6 — AoU-vs-ours ancestry comparison and European classifier

`compare_aou_ancestry.sh` compares this pipeline's ADMIXTURE K=6 fractions to
AoU's v8 auxiliary ancestry outputs:

```text
aux/ancestry/echo_v4_r2.ancestry_preds.tsv
aux/admixture_estimates/aou_admixture_estimates_rye_v8.Q
```

AoU's RYE fractions include `mid` for Middle Eastern ancestry. Our K=6
projection uses the public-statgen global model and has `Oceanian`
instead of `mid`, so Step 6 treats MID as an AoU-specific component rather
than forcing a one-to-one mapping.

The European classifier uses fixed ancestry-fraction thresholds:

```text
European >= 0.8
African <= 0.1
American <= 0.1
East_Asian <= 0.1
Oceanian <= 0.1
South_Asian: no cap
```

The PLINK keep-list for downstream EUR-only analyses is written to:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/europeans/classified_european_iids.txt
```

The AoU-vs-ours comparison directory is:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_vs_ours/
```

Key summary outputs include:

```text
aou_vs_ours.summary.tsv
component_pair_metrics.tsv
component_correlation_matrix.tsv
fraction_distribution_summary.tsv
aou_pred_counts.tsv
aou_pred_vs_ours_european_summary.tsv
european_set_overlap_summary.tsv
european_set_group_ancestry_summary.tsv
discordant_set_component_summary.tsv
european_set_group_counts_by_aou_pred.tsv
aou_mid_threshold_summary.tsv
aou_mid_threshold_component_summary.tsv
aou_mid_bin_counts.tsv
aou_mid_bin_component_summary.tsv
```

The comparison also writes plots under `aou_vs_ours/plots/`, including
component scatter plots, ancestry-fraction distributions, a hard-call vs
dominant-component heatmap, European set-overlap counts, discordant European
call composition plots, and AoU MID-threshold composition plots.

### Step 7 — KING kinship and close relationship classification

Kinship starts from the high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

`subset_kinship_snps.sh` downloads the public UK Biobank SNP QC file
`ukb_snp_qc.txt` if it is not already cached locally, keeps SNPs where
`in_Relatedness == 1`, intersects those rsids with the HQ direct bfile, and
then submits a dsub worker to compute all-sample variant missingness on that
intersection. The final KING SNP list keeps variants with missingness strictly
less than `KINSHIP_MISSING_MAX` (default: `0.01`).

This is a **variant missingness** filter, not a sample missingness filter. No
`--mind`/sample filtering is applied before KING; the KING run uses all samples
present in `direct_bfile_hq/chr1_22_merged_hq`. The only sample-level output at
this stage is the Step 4 missingness report, which is diagnostic unless a later
pipeline step explicitly consumes it.

The subset stage writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_subset_all_sample_missingness.vmiss
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_missingness_threshold_counts.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_subset_summary.tsv
```

The summary reports the key set-operation counts:

```text
ukb_relatedness_snps
n_intersection_hq_direct
n_intersection_and_missing_lt_0.01
```

The threshold-count table reports how many SNPs would pass common missingness
cutoffs (`0.05`, `0.04`, `0.03`, `0.02`, `0.01`) from the same `.vmiss` file.
If the `.vmiss` file already exists and only `KINSHIP_MISSING_MAX` changes,
`subset_kinship_snps.sh` re-filters locally instead of re-submitting the dsub
job and re-localizing the large HQ direct bfile.

`get_genotypes.sh` proceeds directly from the SNP subset to the KING run. KING
is the largest single compute step in the pipeline, but it is idempotent: once
matching KING outputs exist, later runs skip it. To tighten the all-sample SNP
missingness threshold before KING, override `KINSHIP_MISSING_MAX` before
launching the pipeline:

```bash
# Override the default all-sample SNP missingness threshold.
KINSHIP_MISSING_MAX=0.02 bash get_genotypes.sh 2>&1
```

`run_king_kinship.sh` runs:

```text
plink2 --make-king-table --king-table-filter 0.035
```

on the HQ direct bfile with the final SNP extract list. It does not
materialize a separate kinship-only bfile; the exact SNP set is captured by
`ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt` and the parameter
files. KING writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/aou_hq_direct_rel.kin0
```

`kinship_qc.sh` compares our KING kinship coefficients against AoU's provided
KING-style relatedness table:

```text
v8/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv
```

AoU's table provides the pairwise kinship coefficient, so this comparison is
kinship-only; IBS0 diagnostics are available only from our PLINK2/KING output.
Summary tables and plots are written under:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/qc/
```

`classify_relations.sh` then applies fixed KING kinship/IBS0 thresholds on our
KING table:

```text
identical:    kinship >= 0.3535 and IBS0 < 0.0012
parent_child: 0.1767 <= kinship < 0.3535 and IBS0 < 0.0012
sibling:      0.1767 <= kinship < 0.3535 and IBS0 >= 0.0012
```

The AoU version deliberately does not apply a birth-year/month sibling age-gap
filter because no portable AoU phenotype dependency is part of this genotype
pipeline.

Validation run summary from the first completed AoU v8 run in this workspace
with `KINSHIP_MISSING_MAX=0.01` and `KING_TABLE_FILTER=0.035`:

```text
UKBB in_Relatedness SNPs:                         93,511
Intersection with HQ direct bfile:                85,223
Final KING SNPs after all-sample missingness <1%: 84,550

KING pairs reported at kinship >=0.035:           78,142
AoU provided relatedness pairs:                   39,681
Overlapping pair set:                             39,575
Pearson r vs AoU kinship:                         0.9894
Mean absolute kinship difference:                 0.00893

Close relationships, kinship >=0.1767:            26,215
  sibling:                                         7,826
  parent_child:                                   16,140
  identical/twin/duplicate:                        2,249
```

These numbers are a validation/accounting record for that run; they are not
hardcoded into the pipeline. New AoU releases or changed thresholds should be
summarized from the files in `sbayesrc_genotypes/kinship/`.

### Step 8 — Select unrelated European IIDs for PCA fitting

`select_pca_europeans.sh` selects the European subset used to fit ancestry
PCA. It starts from this pipeline's European keep-list:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/europeans/classified_european_iids.txt
```

It then identifies European samples involved in `sibling` or `identical`
relationships from:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/close_relations.csv
```

Those samples are the seed exclusion set. The script expands that seed set to
include everyone directly related to those seeds at a third-degree kinship
threshold:

```text
KINSHIP >= 0.0441941
```

using:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/aou_hq_direct_rel.kin0
```

The candidate PCA set is:

```text
classified Europeans - expanded sibling/identical-relative exclusion set
```

Finally, the script runs PLINK2's maximal unrelated selector on the candidate
set:

```text
plink2 --psam candidate_pca_europeans.psam \
       --king-cutoff-table aou_hq_direct_rel.kin0 0.0441941
```

This final cutoff guarantees that no retained PCA-fitting pair has KING
kinship at or above the third-degree threshold. The step writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/select_pca_europeans.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_eur_log.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/seed_sibling_identical_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/expanded_sibling_identical_relatives_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/candidate_pca_europeans_iids.txt
```

The summary reports the number of classified Europeans, seed samples,
expanded exclusions, candidate PCA Europeans before the final KING cutoff,
samples removed by the final cutoff, and final `fit_pca_iids`.

Validation run summary from the first completed AoU v8 run in this workspace
with `PCA_KINSHIP_THRESHOLD=0.0441941` and
`PCA_SEED_RELATIONSHIPS=sibling,identical`:

```text
Classified Europeans:                         234,889
Seed sibling/identical relationship rows:      10,075
Samples in seed relationship rows:             17,788
European seed samples:                          7,771
KING graph edges at kinship >=0.0441941:       60,895
KING graph nodes:                              82,914

Expanded exclusions total:                      8,941
Expanded exclusions in Europeans:               8,886
Expanded exclusions not in Europeans:              55
Candidate PCA Europeans before KING cutoff:   226,003
Removed by final KING cutoff:                  11,339
Final PCA fitting IIDs:                       214,664

Verification non-European retained:                 0
Verification expanded exclusions retained:          0
Verification related pairs retained:                0
```

As with the kinship validation counts, these are accounting numbers from that
run, not hardcoded expectations. New AoU releases or changed thresholds should
be summarized from `pca_eur/select_pca_europeans.summary.tsv`.

### Step 9 — QC SNPs for PCA

`pca_snp_qc.sh` builds the bfile used to fit European ancestry PCs. It starts
from the high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

and keeps only the final PCA-fitting IIDs from Step 8:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
```

The PCA-specific filters are applied in this order:

```text
1. keep final PCA-fitting IIDs
2. abs(AoU EUR ALT frequency - SBayesRC ALT frequency) <= 0.03
3. MAF >= 0.01 in the PCA-fitting samples
4. variant missingness <= 0.01
5. sample missingness <= 0.01
6. exclude long-range LD regions, hg38
7. LD prune with --indep-pairwise 1000 80 0.1
```

The tighter ALT-frequency agreement filter reuses Step 4's
`chr1_22_merged_hq.variant_qc.tsv`, which already contains
`abs_alt_freq_diff` after aligning AoU and SBayesRC hg38 alleles. This avoids
recomputing frequency concordance and keeps the Step 9 accounting directly
tied to the Step 4 QC table.

The long-range LD regions are hg38 regions downloaded from the public `plinkQC`
resource by the orchestrator and staged to the private Batch worker.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ready.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_afdiff_0.03.extract.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ld_prune.prune.in
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ld_prune.prune.out
```

`pca_snp_qc.filter_steps.tsv` reports, for every sequential filter, the input
SNP/sample count, dropped SNP/sample count, and remaining SNP/sample count.
`pca_snp_qc.summary.tsv` records the thresholds and final `pca_ready` counts.

Observed Step 9 accounting from the current v8 run:

| Filter | Dropped SNPs | Dropped samples | Remaining SNPs | Remaining samples |
|---|---:|---:|---:|---:|
| Source high-quality direct bfile | 0 | 0 | 498,890 | 414,830 |
| Keep Step 8 PCA-fitting IIDs | 0 | 200,166 | 498,890 | 214,664 |
| ALT-frequency difference <= 0.03 | 622 | 0 | 498,268 | 214,664 |
| MAF >= 0.01 | 3,579 | 0 | 494,689 | 214,664 |
| Variant missingness <= 0.01 | 3,313 | 0 | 491,376 | 214,664 |
| Sample missingness <= 0.01 | 0 | 83 | 491,376 | 214,581 |
| Exclude hg38 long-range LD regions | 8,026 | 0 | 483,350 | 214,581 |
| LD prune `1000 80 0.1` | 348,370 | 0 | 134,980 | 214,581 |

These numbers are run accounting, not hardcoded expectations. Recompute them
from `pca_snp_qc.filter_steps.tsv` after changing thresholds or moving to a
new AoU release.

### Step 10 — Fit PCA and project all samples

`fit_project_pca.sh` fits European ancestry principal components and projects
them onto the full high-quality direct bfile. The fit set is the Step 9
PCA-ready bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ready
```

That bfile contains the unrelated European PCA-fitting samples from Step 8
after the Step 9 SNP/sample QC and LD pruning. PCA is fit with PLINK2's
approximate PCA algorithm:

```text
plink2 --bfile pca_ready --pca allele-wts 20 approx --seed 0
```

The projection target is the all-sample high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

The worker computes allele counts in the PCA fit set, verifies that every PCA
SNP is present in the projection bfile, parses the `A1` and PC-weight columns
from `aou_pcs.eigenvec.allele`, and projects all samples with PLINK2
`--score` using `no-mean-imputation` and `variance-standardize`.

PLINK2 writes `aou_pcs.eigenvec.allele` with one score row per allele. For
biallelic SNPs, that means two allele-weight rows per SNP; the pipeline
therefore validates the number of unique SNP IDs in that file, not the raw row
count.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.eigenval
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.eigenvec
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.eigenvec.allele
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_eur_counts.acount
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_projected.sscore
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_project_pca.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_project_pca.params.tsv
```

The summary reports fitted variants, PCA-fit samples, projected samples, PCs,
seed, score-file columns, and PCA SNPs missing from the projection bfile.
Missing projection SNPs should be zero because Step 9 starts from
`direct_bfile_hq`.

Observed Step 10 accounting from the current v8 run:

| Metric | Value |
|---|---:|
| PCA SNPs used for fitting/projection | 134,980 |
| PCA-fitting samples | 214,581 |
| Allele-weight rows | 269,960 |
| Unique allele-weight SNP IDs | 134,980 |
| Projected samples | 414,830 |
| Missing PCA SNPs in projection bfile | 0 |
| PCs | 20 |

Eigenvalues from that run:

```text
437.348
89.0169
70.7415
30.5331
24.0492
14.2525
10.5199
8.24073
8.02872
7.47943
7.19182
6.86026
6.7948
6.75767
6.74678
6.71023
6.69641
6.68407
6.6745
6.66575
```

### Step 11 — Sex covariate and sex/ploidy QC

`get_genetic_sex.sh` builds the binary sex covariate used by the height GWAS
example. It queries sex at birth from the AoU CDR `person` table, joins WGS
sex-ploidy metrics from AoU genomic QC, and keeps only samples with confident
binary sex:

```text
sex at birth is Female or Male
AND WGS sex ploidy is XX or XY
AND sex at birth and WGS ploidy are concordant
```

The output covariate uses FIDs from `direct_bfile_hq/chr1_22_merged_hq.fam`
and `sex_01` coding `0=female, 1=male`. The height setup step centers this as
`sex_c = sex_01 - 0.5`, so female is `-0.5` and male is `0.5`.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/sex_covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/genetic_sex_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/sex_ploidy_crosstab.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/genetic_sex_log.txt
```

Observed Step 11 accounting from the current v8 run:

| Metric | Value |
|---|---:|
| Sample universe | 414,830 |
| Binary sex at birth | 410,445 |
| Missing/non-binary sex at birth | 4,385 |
| Missing WGS sex ploidy | 1 |
| Non-canonical WGS sex ploidy | 1,451 |
| Concordant confident binary samples | 408,993 |
| Confident female | 249,842 |
| Confident male | 159,151 |
| Binary sex/ploidy discordances | 0 |
| Confident sample percent | 98.592918% |

### Step 12 — Sample-QC exclusions for anomalous identical components

`build_identical_component_sample_qc.sh` builds a conservative GWAS exclusion
list from the close-relationship calls generated in Step 7. It treats
`relationship == identical` pairs as graph edges and computes connected
components. Ordinary monozygotic twin relationships are expected to appear as
simple size-2 components; components of size three or larger are not plausible
ordinary twin pairs and are excluded from downstream GWAS sample sets.

The default rule is:

```text
exclude all samples in identical-genotype connected components with size >= 3
```

The threshold can be changed with:

```bash
IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE=3
```

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_components.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_sample_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_sample_qc.log
```

Observed Step 12 accounting from the current v8 run:

| Metric | Value |
|---|---:|
| Identical KING pairs | 2,249 |
| Unique samples in identical graph | 3,860 |
| Identical components | 1,896 |
| Size-2 components | 1,854 |
| Size-3 components | 36 |
| Size-4 components | 4 |
| Size-5 components | 1 |
| Size-23 components | 1 |
| Components size >=3 | 42 |
| Samples excluded | 152 |

### Step 13 — Final GWAS genotype inputs

`make_gwas_genotype_inputs.sh` builds the final genotype files consumed by the
height REGENIE example. This step is separate from the earlier
`direct_bfile_hq` build because the final GWAS inputs use stricter,
purpose-specific QC:

```text
GWAS Step 1 source:
  direct_bfile_hq/chr1_22_merged_hq

GWAS Step 1 filters:
  variant missingness <= 0.01 in our classified European samples
  abs(fit_pca ALT frequency - SBayesRC/snp.info ALT frequency) <= 0.03
  MAF >= 0.007 in fit_pca_iids

GWAS Step 2 source:
  wgs_pfiles/chr{1..22}

GWAS Step 2 filters:
  variant missingness <= 0.03 in our classified European samples
  abs(fit_pca ALT frequency - SBayesRC/snp.info ALT frequency) <= 0.04
  MAF >= 0.007 in fit_pca_iids
```

The sample panels are intentionally different:

```text
classified European samples:
  europeans/classified_european_iids.txt
  used only for the genotype-missingness filter

fit_pca_iids:
  pca_eur/fit_pca_iids.txt
  used for ALT-frequency concordance and MAF
```

The SBayesRC/snp.info ALT frequency comes from
`data/support/sbayesrc_liftover_results.csv`. The filter helper matches
`A1_hg38`/`A2_hg38` to each AoU REF/ALT pair, converts `A1Freq` to the AoU ALT
allele, and drops a variant if the liftover row is missing or the alleles do
not match unambiguously.

The step first computes metrics with PLINK2:

```text
plink2 --bfile direct_bfile_hq/chr1_22_merged_hq \
  --keep pca_eur/fit_pca_iids.txt --freq counts
plink2 --bfile direct_bfile_hq/chr1_22_merged_hq \
  --keep europeans/classified_european_iids.txt --missing variant-only

plink2 --pfile wgs_pfiles/chrN \
  --keep pca_eur/fit_pca_iids.txt --freq counts
plink2 --pfile wgs_pfiles/chrN \
  --keep europeans/classified_european_iids.txt --missing variant-only
```

`build_gwas_genotype_filters.py` then applies the ordered filters, writes the
extract lists, writes per-variant QC tables, and writes a combined allele
frequency file for the passing 7M WGS variants. The MAF threshold replaces a
separate MAC threshold; the summary still reports the minimum fit-pca minor
allele count among passing variants as a sanity check.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/gwas_genotype_qc.params.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/gwas_genotype_qc.summary.tsv

${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.variant_qc.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.fit_pca_alt_freqs.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/gwas_step1_direct.filter_steps.tsv

${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/qc/chr{1..22}.variant_qc.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/extracts/chr{1..22}.extract.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/fit_pca_af/chr{1..22}.fit_pca_alt_freqs_passing.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/fit_pca_af/gwas_step2_fit_pca_alt_freqs_passing.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/gwas_step2_wgs.filter_steps.tsv
```

Observed Step 13 accounting from the current v8 run:

| Metric | Step 1 direct | Step 2 WGS |
|---|---:|---:|
| Source variants | 498,890 | 7,349,435 |
| Dropped: liftover missing / allele mismatch | 0 | 0 |
| Dropped: ALT-frequency difference above threshold | 567 | 24,009 |
| Dropped: fit-pca MAF < 0.007 | 1 | 246 |
| Dropped: classified-European missingness above threshold | 3,338 | 73,787 |
| Final variants | 494,984 | 7,251,393 |
| Final samples | 414,830 | 414,830 |
| Minimum fit-pca MAC among retained variants | 3,027 | 2,991 |
| Minimum fit-pca MAF among retained variants | 0.00705088 | 0.00700198 |
| Maximum classified-European missingness among retained variants | 0.00999621 | 0.0299971 |
| Maximum retained absolute ALT-frequency difference | 0.02998992 | 0.03999598 |

The drop counts are ordered/conditional: a variant removed by an earlier rule is
not counted again by later rules.

### Step 14 — Height GWAS input setup

`setup_height_gwas.sh` builds the phenotype, covariate, and keep files for a
height GWAS in all samples classified as European by this pipeline. It uses
AoU program-collected height measurements:

```text
measurement_concept_id = 3036277
measurement_source_concept_id = 903133
measurement_type_concept_id = 44818701
unit_concept_id = 8582
minimum height = 140 cm
```

The BigQuery result is exported to the workspace bucket, processed inside the
AoU environment, and not downloaded off-platform. For each participant, the
helper takes the median valid height and the mean age at height measurement.
It then intersects with:

```text
1. our European ancestry keep-list
2. Step 12 sample-QC exclusion list
3. Step 11 confident sex covariates
4. the high-quality direct-bfile .fam sample IDs
5. Step 10 projected PC scores
```

Covariates are:

```text
age_c                  = age_at_height - mean(age_at_height)
sex_c                  = sex_01 - 0.5
age_c_sex_c_inter      = age_c * sex_c
PC1_AVG ... PC10_AVG   = first 10 projected European PCs
```

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/training_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/phen.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/base_covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/height_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/height_gwas_log.txt
```

Observed Step 14 accounting from the current v8 run:

| Metric | Value |
|---|---:|
| Our classified Europeans | 234,889 |
| Samples in sample-QC exclusion list | 152 |
| Europeans removed by sample QC | 37 |
| Europeans after sample QC | 234,852 |
| Confident sex covariate rows | 408,993 |
| Height query rows after source/min-height filters | 439,858 |
| Projected PC rows | 414,830 |
| Europeans missing a height row | 33,090 |
| Height candidates missing confident sex | 2,819 |
| Height+sex candidates missing fam row | 0 |
| Height+sex+fam candidates missing PCs | 0 |
| Final GWAS samples | 198,943 |
| GWAS female | 118,723 |
| GWAS male | 80,220 |
| Mean height | 169.1081199 cm |
| Median height | 168.3 cm |
| Mean age at height measurement | 55.65723437 years |
| PCs included | 10 |

### Step 15 — Height GWAS with REGENIE

`run_continuous_regenie_gwas.sh` runs the optional continuous-trait REGENIE
height example. It is gated by `RUN_HEIGHT_GWAS=1` because it launches one
Step 1 Batch job plus a 22-task Step 2 Batch array. The default run applies
rank-inverse normal transformation (`--apply-rint`) and uses the Step 14
covariates.

Step 1 uses the final GWAS direct-SNP bfile from Step 13:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1
```

Step 2 tests the final Step 13 WGS pfiles:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}
```

The orchestrator stages a small REGENIE runtime bundle to the workspace bucket
from the Workbench-provided REGENIE binary. The bundle includes the MKL shared
libraries needed on the default Ubuntu Batch worker. For AoU WGS pfiles, the
Step 2 worker also writes a local REGENIE-compatible `.psam` with `FID=0`
because AoU pfiles use a `#IID`-only psam header. It rewrites the localized
Step 1 prediction list so Step 2 reads the localized LOCO prediction file.
Neither change modifies the stored pfiles.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step1/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step2/chr{1..22}/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/lightweight/chr{1..22}.height_example.regenie_lite.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/lightweight/regenie_lite.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/regenie_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/regenie_gwas.params.tsv
```

After Step 2 verification, the shared runner writes compact per-chromosome
association summaries with columns `rsid`, `allele1`, `a1freq`, `n`, `beta`,
`se`, and `log10p`. Set `REGENIE_MAKE_LIGHTWEIGHT_OUTPUTS=0` or pass
`--no-lightweight` to skip these compact files.

The previous full height GWAS output predates Step 13 final genotype
filtering. Re-run with `RUN_HEIGHT_GWAS=1` after Step 13 completes to generate
current REGENIE counts for the stricter final genotype inputs.

### gwas_dev — Additional exploratory GWAS commands

The `gwas_dev` branch adds standalone commands for additional exploratory
continuous-trait REGENIE GWAS runs. These commands are deliberately separate
from `get_genotypes.sh`. They assume the main pipeline has completed through
Step 13, but they do not require the optional height GWAS to have run.

Cheap setup-only checks:

```bash
bash run_ea_gwas.sh --setup-only
bash run_income_gwas.sh --setup-only
```

Full GWAS submissions:

```bash
bash run_ea_gwas.sh
bash run_income_gwas.sh
```

Both commands use:

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
| `1585941` | Highest Grade: Never Attended | 1.0 |
| `1585942` | Highest Grade: One Through Four | 2.5 |
| `1585943` | Highest Grade: Five Through Eight | 6.5 |
| `1585944` | Highest Grade: Nine Through Eleven | 10.0 |
| `1585945` | Highest Grade: Twelve Or GED | 13.0 |
| `1585946` | Highest Grade: College One to Three | 15.0 |
| `1585947` | Highest Grade: College Graduate | 17.0 |
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

The setup scripts write the current workspace-specific sample counts and answer
counts to `{ea,income}_gwas.summary.tsv` and `{ea,income}_answer_counts.tsv`.

### ses_ea_proxy primary setup and scoring

`run_ses_ea_proxy_gwas.sh` builds the primary SES-EA proxy scores and the
matching REGENIE input files. Setup/scoring is the default behavior; it does
not submit REGENIE unless `--run-gwas` is passed explicitly.

Recommended run order:

```bash
# 1. Build genotype/sample-QC/PCA/sex/genotype inputs.
nohup bash get_genotypes.sh > logs/run.log 2>&1 &

# 2. Build SES-EA proxy scores and REGENIE input files.
bash run_ses_ea_proxy_gwas.sh --setup-only

# 3. Build ETM cognitive scores from the proxy-score outputs.
bash run_etm_cog_task_factors.sh --stage-aggregate

# 4. Combine the task scores into one ETM general cognitive/performance score.
bash run_etm_g_from_task_scores.sh --stage-aggregate
```

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

The goal is to produce a non-genetic proxy for the education-attainment
teacher label, then review out-of-sample model performance and covariate
correlations before deciding whether to run GWAS. The setup therefore stops
after score generation by default. REGENIE is only an opt-in follow-up.

Samples are restricted to classified European IIDs, confirmed genetic sex in
`genetic_sex/sex_covar.txt`, samples not in the identical-component size `>=3`
sample-QC exclusion list, codeable EA from The Basics question `1585940`, and
age at that The Basics response `>=26`. Confirmed genetic sex is also included
as an XGBoost model feature.

The score uses a cross-fit design:

```text
1. Split eligible fit_pca_iids into 5 seeded folds.
2. For each fold, fit the EA residualization OLS and XGBoost model on the
   other four folds only, then predict the held-out fold.
3. Fit a sixth model on all eligible fit_pca_iids.
4. Apply the sixth model to eligible classified-European samples that were not
   in fit_pca_iids.
```

The EA teacher label is mapped to years from answers `1585941` through
`1585948`. In each cross-fit fold, EA years are residualized on
`yob_c`, `sex_c`, and `yob_c * sex_c` using only the four-fifths training
pool, then z-scored using that training-pool residual mean and SD. The sixth
model uses the same residualization procedure fit on all eligible
`fit_pca_iids`.

Survey features come from The Basics, Lifestyle, Overall Health, Healthcare
Access & Utilization, Personal and Family Health History, Social Determinants
of Health, and Behavioral Health & Personality. BHP is read from the
off-cycle Mental Health / Well-Being CDR dataset; override
`WORKSPACE_MHWB_CDR` if the dataset name differs. The Washington Group
disability items are sourced from The Basics. ZIP3-derived socioeconomic
features come from `ds_zip_code_socioeconomic`; raw ZIP codes are not used.

Feature extraction is multi-select safe. For each person/question, the setup
selects the latest survey timestamp, keeps all answer rows from that timestamp,
and one-hot encodes every retained `answer_concept_id` for nominal and
multi-select fields. Ordered Likert-style fields are encoded as one ordinal
numeric feature; continuous survey values are parsed as numeric features.

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
and sixth-model applied proxy scores are included. This command does not run a
GWAS.

```bash
bash run_etm_cog_task_factors.sh --stage-aggregate
```

The scorer reads Flanker, GradCPT, Delay Discounting, and Emotional Recognition
from the ETM off-cycle dataset. Set `WORKSPACE_ETM_CDR` or pass
`--etm-dataset PROJECT.DATASET` if the dataset name differs. It computes age at
each task from ETM
`test_start_date_time` and the main CDR `person.birth_datetime`, and it uses the
confirmed genetic-sex `sex_c` already present in
`regenie_input/ses_ea_proxy/base_covar.txt`.

For each task, the scoring order is deliberately fixed:

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
before FA/PCA. Age and confirmed genetic sex are removed only from the final
task score, so the factor/PCA loadings are learned from the task measurements
themselves rather than from age/sex-residualized inputs.

The primary valid-sitting filters are:

| Task | Excluded when any primary condition is true |
|---|---|
| Delay Discounting | `flag_median_rt != 0`, `flag_catch_trials != 0`, or `test_restarted` is true |
| GradCPT | `flag_trial_flags != 0`, `flag_non_response != 0`, `flag_omission_error_rate != 0`, or `test_restarted` is true |
| Flanker | `flag_accuracy != 0`, `flag_trial_flags != 0`, or `test_restarted` is true |
| Emotional Recognition | `flag_median_rtc != 0`, `flag_same_response != 0`, `flag_trial_flags != 0`, or `test_restarted` is true |

`any_timeouts` is retained as a diagnostic rather than an exclusion by default.

The command writes a long diagnostic score table plus a one-row-per-sample
recommended-score table. The recommended scores are:

```text
dd_patience_z_age_sex          # sourced from official -lnk
gradcpt_perf_z_age_sex         # PC1 of dprime + RT consistency/speed
flanker_efficiency_z_age_sex   # Flanker efficiency source selected by diagnostics
emorecog_perf_z_age_sex        # Emotional Recognition PC1 of score + RT consistency/speed
```

The current rationale is:

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
  Use the efficiency score. The one-factor efficiency/interference blend failed
  the loading rule, and the interference split was unstable. The efficiency
  score is selected against the official Flanker score by the predeclared
  simple-score rule.

Emotional Recognition:
  Score the task with the GradCPT-analog PC1 of score, -log(cv_rtc), and
  -log(median_rtc), after fixed transforms, winsorization, and z-scoring. The
  combined PC1 is then residualized on sex_c + age_at_test + age_at_test^2 and
  z-scored. The older per-emotion rate-correct efficiency factor, per-emotion
  accuracy factor, and simple score are still written as diagnostics.
```

Current aggregate diagnostics from the scored cohort are below. These are not
pipeline constants; they are regenerated into the diagnostic TSVs whenever the
command is run on a new cohort/CDR.

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

Repeat valid sittings are uncommon, but they provide a useful test-retest
diagnostic. The diagnostic fits the production scoring recipe on first valid
sittings, applies those same transforms/loadings/age-sex residualization
parameters unchanged to each person's second valid sitting, and then correlates
first-vs-second task scores.

| Task score | Repeat pairs | Pearson r | Spearman r | Mean second - first z |
|---|---:|---:|---:|---:|
| `dd_patience_z_age_sex` | 223 | 0.781 | 0.771 | 0.036 |
| `gradcpt_perf_z_age_sex` | 64 | 0.892 | 0.904 | 0.150 |
| `flanker_efficiency_z_age_sex` | 173 | 0.776 | 0.772 | 0.121 |
| `emorecog_perf_z_age_sex` | 209 | 0.802 | 0.758 | 0.070 |

The repeat gaps are mostly same-day retries rather than year-scale retests. The
right tail is long for DD, Flanker, and Emotional Recognition, but the median
second sitting is only minutes after the first:

| Task | Repeat pairs | Median gap | P10 gap | P75 gap | P90 gap | P95 gap | Same-day repeats |
|---|---:|---:|---:|---:|---:|---:|---:|
| DD | 223 | 0.003 days, about 4.9 min | 0.000 days | 47.3 days | 157.5 days | 204.1 days | 65.0% |
| GradCPT | 64 | 0.007 days, about 10.6 min | 0.000 days | 0.015 days | 2.1 days | 7.0 days | 84.4% |
| Flanker | 173 | 0.017 days, about 23.8 min | 0.000 days | 15.4 days | 153.4 days | 204.1 days | 57.2% |
| Emotional Recognition | 209 | 0.004 days, about 5.7 min | 0.000 days | 20.1 days | 157.4 days | 209.9 days | 69.4% |

These correlations are encouraging, but they should be read as short-interval
repeatability more than long-term stability because most repeats happened on
the same day.

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
the task-score output and the SES-EA proxy cohort files only; it does not query
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

One-task scores are therefore intentionally shrunk toward zero; the command does
not divide by the sum of available weights or z-score within missingness
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

The four-domain score is accepted by the loading checks and adds coverage for
people with only Emotional Recognition observed, but in the current run it is
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

## Prerequisites

You are inside an AoU Verily Jupyter session terminal (the standard
interactive analysis environment for the All of Us Researcher Workbench).
The session provides everything the pipeline needs:

- **Env vars** set automatically by the Workbench: `GOOGLE_PROJECT`,
  `CDR_STORAGE_PATH`. (`$WORKSPACE_BUCKET` is unreliable on some pods —
  `get_genotypes.sh` derives the writable bucket URI from the FUSE mount
  table instead, which is portable across users.)
- **Controlled-tier dataset bucket** FUSE-mounted (read-only) at
  `/home/jupyter/workspace/data_controlled/vwb-aou-datasets-controlled/`.
- **Workspace bucket** FUSE-mounted (read-write) at
  `/home/jupyter/workspace/workspace-bucket/`.
- **`plink2`** preinstalled at `/opt/workbench-tools/binaries/bin/plink2`.
  The orchestrator stages this binary to the workspace bucket once, and each
  Batch worker `--input`s it back.
- **`regenie`** preinstalled at `/opt/workbench-tools/binaries/bin/regenie`
  for the optional height GWAS example.
- **`dsub`** preinstalled at `/opt/conda/envs/jupyter/bin/dsub`.
- **Internet access from the Jupyter pod** for public reference downloads.
  Batch workers use private networking, so the orchestrator downloads public
  ADMIXTURE inputs and stages them to the workspace bucket.
- **gcloud Application Default Credentials** (run `gcloud auth
  application-default login` once if `gcloud storage` or dsub return auth
  errors).
- **Pet service account** of the pod, used as the Batch worker identity
  (auto-detected via `gcloud config get-value account`).

Python dependencies for the local helper scripts are listed in
`requirements.txt` and auto-installed by `get_genotypes.sh`. `dsub` is provided
by the Workbench image.

## Usage

Smoke-test one chromosome of the WGS extraction step (one dsub task, ~15–25 min;
the direct-bfile step is skipped because it requires all 22 chromosomes):
```bash
SBAYESRC_TEST_CHROM=22 bash get_genotypes.sh 2>&1
```

Full all-22 run (recommended in the background):
```bash
nohup bash get_genotypes.sh > logs/run.log 2>&1 &
```

Run the optional height GWAS example as well:
```bash
RUN_HEIGHT_GWAS=1 nohup bash get_genotypes.sh > logs/run_height_gwas.log 2>&1 &
```

The default `bash get_genotypes.sh` command is designed to run the full
pipeline end to end. Expensive steps, including KING and the optional REGENIE
height GWAS, have parameter/count-based idempotency checks so later runs skip
matching outputs instead of re-submitting work.

The run logs to `logs/run_YYYYMMDD_HHMMSS.log` (timestamped) and tees through
to the foreground if attached. Each Batch worker's stdout/stderr is uploaded
to `${WORKSPACE_BUCKET}/sbayesrc_genotypes/logs/dsub/`.

Steps that query the AoU CDR use a temporary BigQuery dataset in the user's
workspace project. The pipeline uses `SBAYESRC_BQ_TMP_DATASET` if set;
otherwise it chooses an existing dataset in this order: `sbayesrc_tmp`,
`high_quality_cohort`, then the first dataset returned by `bq ls`. The AoU
pet service account may not be allowed to create new BigQuery datasets, so set
`SBAYESRC_BQ_TMP_DATASET` to an existing writable dataset if auto-selection is
not appropriate for your workspace.

## Portability

No user-specific identifiers are hardcoded. The pipeline derives:

- the **workspace bucket URI** from the FUSE mount table (not `$WORKSPACE_BUCKET`, which on some pods points to a non-existent cloned-bucket name);
- the **GCP project** from `$GOOGLE_PROJECT` (set by the Workbench);
- the **pet service account** from `gcloud config get-value account` (specific to each pod);
- the **VPC network/subnetwork** as `projects/$GOOGLE_PROJECT/global/networks/network` and `projects/$GOOGLE_PROJECT/regions/us-central1/subnetworks/subnetwork` — these are the canonical AoU project VPC names and are the same for every user;
- the **AoU controlled-tier path** as a stable gs:// URI shared across users.

Any AoU researcher with controlled-tier access should be able to clone the
repo and run it as-is.

## Idempotency

- Step 1 skips chromosomes whose local `chr{N}.{extract,idmap}.txt` already
  exist non-empty.
- Step 2 skips chromosomes whose `chr{N}.pgen` already exists at the
  destination (cheap gcsfuse metadata-only lookup via `test -f` on the
  workspace bucket mount — no API call needed for the skip check).
- Step 3 skips direct pfiles and the merged direct bfile when they already
  exist with the expected variant counts.
- Step 4 skips the high-quality direct bfile when its parameter file, summary,
  bfile, and sample-missingness summary already exist with matching thresholds
  and expected variant counts. The EUR frequency/missingness metrics are also
  skipped when their output line counts match the raw direct bfile.
- Step 5 skips ADMIXTURE prep, split, projection batches, and final concat
  independently when their parameter files and expected row counts match.
  Projection only submits missing or stale batch `.Q` files.
- Step 6 skips the AoU-vs-ours ancestry comparison when its parameter file,
  summary tables, key plots, and European keep-list already exist with matching
  input sizes and thresholds.
- Step 7 skips the kinship SNP subset, KING run, AoU comparison, and close
  relationship classifier independently when their parameter files, thresholds,
  and expected counts match.
- Step 8 skips the PCA European selector when its parameter file, summary, and
  `fit_pca_iids.txt` exist with matching input sizes, threshold, seed
  relationships, and output count.
- Step 9 skips PCA SNP QC when its parameter file, summary, filter-step table,
  and `pca_ready.{bed,bim,fam}` exist with matching input sizes, thresholds,
  high-LD-region checksum, and final counts.
- Step 10 skips PCA fitting/projection when its parameter file, summary,
  eigenvalue/eigenvector/allele-weight files, allele counts, and projected
  score file exist with matching input sizes, PC settings, and output counts.
- Step 11 skips the sex covariate/QC build when its parameter file, summary,
  crosstab, log, and `sex_covar.txt` exist with matching input sizes and
  concordance settings.
- Step 12 skips identical-component sample QC when its parameter file, summary,
  component tables, exclusion list, and log exist with matching input sizes
  and exclusion threshold.
- Step 13 skips final GWAS genotype construction when its parameter file,
  summary, extract lists, Step 1 bfile, and Step 2 per-chromosome pfiles exist
  with matching input sizes, thresholds, and expected variant counts. It also
  skips metric-scan dsub jobs independently when their output line counts match
  the source genotype files and sample keep-list sizes.
- Step 14 skips height GWAS input setup when its parameter file, summary,
  phenotype, covariate, keep, and log files exist with matching concept IDs,
  minimum-height threshold, PC count, input sizes, and sample-QC exclusion
  list.
- Step 15 is skipped unless `RUN_HEIGHT_GWAS=1`. When enabled, it skips
  REGENIE Step 1 if the prediction list, params, and summary match the desired
  run parameters. It then submits only Step 2 chromosomes missing matching
  per-chromosome summary/params files, and writes a final GWAS summary only
  after all selected chromosomes validate.

A re-run of `get_genotypes.sh` after a successful run should skip every step
and submit zero dsub tasks.

## Files

| File | Role |
|---|---|
| `get_genotypes.sh` | Orchestrator. Sets paths/env (including the dsub knobs) and calls each step script. |
| `store_sbayesrc_ids.py` | Step 1 — local SBayesRC ID + idmap generation. |
| `wgs_extract_variants.sh` | Step 2 — submits the dsub `--tasks` TSV (one row per missing chrom), waits, polls dstat until every task is terminal, verifies bucket outputs, writes combined summary. |
| `dsub_extract_worker.sh` | Per-Batch-worker script — the four-pass extract for one chromosome. |
| `normalize_pvar_alleles.py` | Canonical, testable reference for the pass-1.5 right-trim logic that's inlined as awk in the worker. Not used at runtime — kept for unit-testability and local debugging. |
| `prepare_direct_snps.py` | Step 3a — prepares per-chromosome direct-SNP extraction lists and missing-SNP reports. |
| `extract_direct_snps.sh` | Step 3b — extracts present direct SNPs from WGS pfiles. |
| `make_direct_bfile.sh` | Step 3c — submits/verifies the direct-bfile dsub merge. |
| `dsub_direct_bfile_worker.sh` | Per-Batch-worker script — merges direct pfiles into a sorted intermediate pgen, converts to `chr1_22_merged.{bed,bim,fam}`, and writes a summary. |
| `make_hq_direct_bfile.sh` | Step 4 — creates the AoU EUR keep-list, submits EUR frequency/missingness metrics, applies high-quality direct-SNP filters, builds the final filtered bfile, and verifies outputs. |
| `dsub_hq_direct_metrics_worker.sh` | Step 4 metrics worker — computes EUR allele frequencies and variant missingness on the raw direct bfile. |
| `filter_hq_direct_snps.py` | Step 4 local filter builder — joins EUR metrics to SBayesRC liftover frequencies and writes the ordered filter summary, variant QC table, and final extract list. |
| `dsub_hq_direct_bfile_worker.sh` | Step 4 bfile worker — extracts the passing variants, writes `chr1_22_merged_hq.{bed,bim,fam}`, and computes sample missingness over the final variant set. |
| `admixture_prep.sh` | Step 5a — stages ADMIXTURE inputs and submits the prep/alignment worker. |
| `dsub_admixture_prep_worker.sh` | Step 5a worker — applies all-sample `geno <= 0.05`, intersects reference SNPs, aligns alleles, and writes `aou_admixture_aligned` plus `ref_aligned.P`. |
| `admixture_align_alleles.py` | Step 5a reference helper — readable/testable Python version of the allele-alignment logic; the Batch worker uses an inline awk implementation to avoid depending on python inside the worker image. |
| `admixture_split_batches.sh` | Step 5b — submits the batch-splitting worker. |
| `dsub_admixture_split_worker.sh` | Step 5b worker — creates 20,000-person batch bfiles and per-batch `.P.in` files. |
| `admixture_run_projection.sh` | Step 5c — submits projection tasks and the final concat job. |
| `dsub_admixture_project_worker.sh` | Step 5c worker — runs ADMIXTURE `-P` for one batch and verifies Q row counts. |
| `dsub_admixture_concat_worker.sh` | Step 5c worker — concatenates batch Q files and writes `aou_admixture_k6.tsv`. |
| `compare_aou_ancestry.sh` | Step 6 — idempotent wrapper for the AoU-vs-ours ancestry comparison and European keep-list. |
| `compare_aou_ancestry.py` | Step 6 helper — joins AoU hard ancestry calls, AoU RYE fractions, and our ADMIXTURE fractions; writes aggregate tables and plots. |
| `subset_kinship_snps.sh` | Step 7a — intersects UKBB `in_Relatedness == 1` SNPs with the HQ direct bfile and applies all-sample missingness `< KINSHIP_MISSING_MAX`. |
| `dsub_kinship_subset_worker.sh` | Step 7a worker — computes variant missingness and writes the final KING SNP extract list plus count summary. |
| `run_king_kinship.sh` | Step 7b — submits/verifies the PLINK2 KING run from the HQ direct SNP subset. |
| `dsub_king_kinship_worker.sh` | Step 7b worker — runs `plink2 --make-king-table --king-table-filter 0.035`. |
| `kinship_qc.sh` | Step 7c — idempotent wrapper comparing our KING output to AoU's provided relatedness table. |
| `kinship_qc.py` | Step 7c helper — writes kinship comparison summaries, pair-level overlap data, scatter plots, and Bland-Altman plots. |
| `classify_relations.sh` | Step 7d — idempotent wrapper for close-relationship classification. |
| `classify_relations.py` | Step 7d helper — classifies sibling, parent_child, and identical pairs from KING kinship/IBS0 thresholds. |
| `select_pca_europeans.sh` | Step 8 — idempotent wrapper selecting unrelated European IIDs for fitting PCA. |
| `select_pca_europeans.py` | Step 8 helper — expands European sibling/identical seeds to their third-degree relatives, removes them, and runs PLINK2 `--king-cutoff-table` for the final unrelated set. |
| `pca_snp_qc.sh` | Step 9 — stages PCA SNP QC inputs and submits/verifies the Batch worker that builds `pca_ready.{bed,bim,fam}`. |
| `dsub_pca_snp_qc_worker.sh` | Step 9 worker — applies the PCA-fitting sample keep-list, tighter frequency concordance, MAF/geno/mind filters, long-range-LD exclusion, and LD pruning. |
| `fit_project_pca.sh` | Step 10 — stages PCA inputs and submits/verifies the Batch worker that fits European PCs and projects them to all samples. |
| `dsub_fit_project_pca_worker.sh` | Step 10 worker — fits PLINK2 PCA with allele weights, computes fit-set allele counts, verifies projection SNP coverage, and writes all-sample PC scores. |
| `get_genetic_sex.sh` | Step 11 — queries AoU sex at birth, joins WGS sex ploidy, and writes the confident binary sex covariate plus QC summaries. |
| `get_genetic_sex.py` | Step 11 helper — builds `sex_covar.txt`, the sex/ploidy crosstab, summary, and verification log. |
| `build_identical_component_sample_qc.sh` | Step 12 — builds the sample-QC exclusion list for identical-genotype components of size `IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE` or larger. |
| `build_identical_component_sample_qc.py` | Step 12 helper — computes identical-pair connected components, writes component tables, and writes the exclusion FID/IID file. |
| `make_gwas_genotype_inputs.sh` | Step 13 — computes GWAS-specific allele frequency/missingness metrics, builds final Step 1/Step 2 extract lists, and submits the final genotype-build workers. |
| `build_gwas_genotype_filters.py` | Step 13 helper — joins fit-pca allele counts, classified-EUR missingness, and SBayesRC/snp.info frequencies; writes filter summaries, variant QC tables, extract lists, and passing Step 2 fit-pca ALT frequencies. |
| `dsub_gwas_direct_metrics_worker.sh` | Step 13 worker — computes fit-pca allele counts and classified-EUR missingness for `direct_bfile_hq`. |
| `dsub_gwas_wgs_metrics_worker.sh` | Step 13 worker — computes fit-pca allele counts and classified-EUR missingness for one WGS pfile chromosome. |
| `dsub_gwas_step1_direct_worker.sh` | Step 13 worker — extracts the final REGENIE Step 1 bfile from `direct_bfile_hq`. |
| `dsub_gwas_step2_wgs_worker.sh` | Step 13 worker — extracts one final REGENIE Step 2 WGS pfile chromosome. |
| `setup_height_gwas.sh` | Step 14 — queries program-collected AoU height, exports the result inside the workspace bucket, and builds REGENIE phenotype/covariate/keep files after sample-QC exclusions. |
| `setup_height_gwas.py` | Step 14 helper — intersects Europeans, sample-QC exclusions, height, confident sex, genotype IDs, and projected PCs; centers covariates; writes summaries and verification checks. |
| `run_continuous_regenie_gwas.sh` | Step 15 — optional continuous-trait REGENIE runner using the final Step 13 genotype inputs. |
| `make_lightweight_regenie_outputs.py` | Step 15 helper — converts full per-chromosome REGENIE outputs into compact `rsid/allele/frequency/effect/SE/log10p` TSVs. |
| `dsub_regenie_step1_worker.sh` | Step 15 worker — runs REGENIE Step 1, writes LOCO predictions, and verifies sample/variant counts. |
| `dsub_regenie_step2_worker.sh` | Step 15 worker — runs one REGENIE Step 2 chromosome, adapting AoU `#IID` psam headers and localized Step 1 prediction paths for REGENIE. |
| `run_ses_ea_proxy_gwas.sh` | Optional phenotype/GWAS command downstream of `get_genotypes.sh`; builds SES-EA proxy scores by default and runs REGENIE only with `--run-gwas`. |
| `setup_ses_ea_proxy_gwas.sh` | SES-EA proxy setup wrapper; checks required `get_genotypes.sh` outputs, extracts survey features, trains/saves XGBoost models, and writes REGENIE inputs. |
| `setup_ses_ea_proxy_gwas.py` | SES-EA proxy helper; implements multi-select-safe feature encoding, missing-data handling, cross-fit teacher residualization, XGBoost training, and diagnostics. |
| `run_etm_cog_task_factors.sh` | Optional cognitive-score command downstream of SES-EA proxy setup; creates recommended ETM cognitive scores and aggregate diagnostics without running GWAS. |
| `score_etm_cog_task_factors.py` | ETM cognitive scoring helper; implements task QC, first-valid sitting selection, final age/sex norming, DD `-lnk`, GradCPT PC1, Flanker efficiency, and Emotional Recognition score/RT PC1 outputs. |
| `cognitive_task_factor_score_spec.md` | Specification for the ETM cognitive scoring method and final recommended-score contract. |
| `run_etm_g_from_task_scores.sh` | Optional downstream command that combines recommended ETM task scores into one ETM general cognitive/performance factor and aggregate diagnostics. |
| `score_etm_general_factor.py` | ETM-g helper; fits the complete-case one-factor model, scores observed-task patterns with regression weights, and writes diagnostics. |
| `etm_general_factor_score_spec.md` | Specification for the downstream ETM general cognitive/performance factor score. |
| `requirements.txt` | Python dependencies for the local helper scripts. |
| `CLAUDE.md` | Project conventions, AoU platform notes, portability rules, dsub-from-Jupyter recipe. Gitignored — local developer reference. |
