# AoU SBayesRC GWAS pipeline

A reproducible, idempotent pipeline that builds per-chromosome PLINK2 pfiles
holding the ~7.35 M SBayesRC SNPs from the All of Us v8 ACAF WGS callset, with
the variant `ID` column populated with SBayesRC rsids. It also builds the
~501k direct-SNP bfile used as the REGENIE step-1 marker backbone, plus a
higher-quality direct-SNP bfile filtered with AoU EUR frequency/missingness
metrics, and runs ADMIXTURE K=6 ancestry projection.

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
```

Each `summary.tsv` reports `requested / src_variants / src_samples /
biallelic_total / extracted / out_samples / missing / remapped / unmapped`
for that chromosome. A combined summary across all processed chromosomes is
also written to `logs/sbayesrc_extract_summary.tsv` at the end of each run.

Downstream UKBB-analog steps after ADMIXTURE (kinship, PCA, REGENIE phenotype
setup) are not yet ported — see
`reference/ukbb-sbayesrc-gwas/get_genotypes.sh` for the full target pipeline
(gitignored locally; cloned for reference).

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

The run logs to `logs/run_YYYYMMDD_HHMMSS.log` (timestamped) and tees through
to the foreground if attached. Each Batch worker's stdout/stderr is uploaded
to `${WORKSPACE_BUCKET}/sbayesrc_genotypes/logs/dsub/`.

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
| `requirements.txt` | Python dependencies for the local helper scripts. |
| `CLAUDE.md` | Project conventions, AoU platform notes, portability rules, dsub-from-Jupyter recipe. Gitignored — local developer reference. |
| `reference/ukbb-sbayesrc-gwas/` | UKBB analog this pipeline mirrors. Gitignored — clone locally for reference only. |
