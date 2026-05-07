# AoU SBayesRC GWAS pipeline

A reproducible, idempotent pipeline that builds per-chromosome PLINK2 pfiles
holding the ~7.35M SBayesRC SNPs from the All of Us v8 ACAF WGS callset, with
the variant `ID` column populated with SBayesRC rsids. Each pipeline step is a
separate sub-script orchestrated by `get_genotypes.sh`. This is the AoU analog
of the early phase of the UKBB SBayesRC GWAS pipeline in
`reference/ukbb-sbayesrc-gwas/`.

## What this pipeline currently produces

`${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/wgs_pfiles/chr{1..22}.{pgen,pvar,psam}`
plus per-chromosome `chr{N}.summary.tsv` files with extraction counts.
`WORKSPACE_BUCKET_URI` is auto-derived from
`wb resource resolve --id=workspace-bucket` (e.g.
`gs://workspace-bucket-<workspace-project-id>`).

Subsequent steps from the UKBB analog (kinship, ADMIXTURE, PCA, REGENIE) are
not yet ported — see `reference/ukbb-sbayesrc-gwas/get_genotypes.sh` for the
full target pipeline.

## Pipeline steps

1. **`store_sbayesrc_ids.py`** — Read `data/support/sbayesrc_hg38.csv`
   (chrom, pos, ref, alt, rsid) and write per-chromosome
   `data/sbayesrc_ids/chr{N}.extract.txt` (one `chr{N}:pos:REF:ALT` ID per
   line) and `chr{N}.idmap.txt` (`chr{N}:pos:REF:ALT \t rsid`). The ID format
   matches AoU's `acaf_threshold/plink_bed/chr{N}.bim` second column exactly,
   so `plink2 --extract` is a direct match.

2. **`upload_sbayesrc_ids.sh`** — Upload both files per chromosome to
   `${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/sbayesrc_ids/` via
   `gcloud storage cp` with `--billing-project` for requester-pays.

3. **`wgs_extract_variants.sh`** — For each autosome, submit a dsub
   (Google Batch) job that runs:
   ```
   plink2 \
     --bed/--bim/--fam <AoU acaf_threshold/plink_bed/chr{N}.{bed,bim,fam}> \
     --extract chr{N}.extract.txt \
     --update-name chr{N}.idmap.txt \
     --no-pheno --make-pgen --out chr{N}
   ```
   Outputs `chr{N}.{pgen,pvar,psam}` + `chr{N}.summary.tsv` to
   `${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/wgs_pfiles/`. All 22 jobs run
   in parallel; the script blocks until all finish, then streams the summary
   TSVs back into a combined local report at
   `logs/sbayesrc_extract_summary.tsv`.

   No QC filters are applied here. ACAF threshold already enforces MAC ≥100
   per ancestry; missingness QC can be layered later.

## Prerequisites

- `wb` CLI authenticated (`wb auth login`).
- Google Application Default Credentials set up
  (`gcloud auth application-default login`).
- An active workspace selected (`wb workspace set --id=<your-workspace-id>`).
- IAM permission to submit Google Batch jobs in the workspace's GCP project.

Python dependencies (`pandas`, `dsub`) are installed automatically by
`get_genotypes.sh` via `pip install -r requirements.txt`.

## Usage

```
bash get_genotypes.sh 2>&1
```

Run in the background — it submits jobs and waits for them, which can take a
while. Logs go to `logs/run_<timestamp>.log` and per-chromosome
`logs/dsub/chr{N}.log`.

## Portability

No user-specific identifiers are hardcoded. The orchestrator derives:

- `GOOGLE_PROJECT` — from `wb workspace describe --format=json`
- `WORKSPACE_BUCKET_URI` — from `wb resource resolve --id=workspace-bucket`

Override either via environment variable if auto-detection fails.

## Idempotency

- Step 1 skips chromosomes whose local `chr{N}.{extract,idmap}.txt` already
  exist.
- Step 2 skips uploads whose destination object already exists
  (`gcloud storage ls`).
- Step 3 skips chromosomes whose `chr{N}.pgen` already exists at the
  destination — never submits a dsub job just to check.

A re-run of `get_genotypes.sh` after a successful run should skip every step
and submit zero jobs.

## Files

- `get_genotypes.sh` — orchestrator; sets env vars and calls the three step
  scripts.
- `store_sbayesrc_ids.py` — local SBayesRC ID + idmap generation.
- `upload_sbayesrc_ids.sh` — local→GCS upload with idempotency.
- `wgs_extract_variants.sh` — per-chromosome dsub extraction submission +
  wait + summary.
- `requirements.txt` — pandas (for `store_sbayesrc_ids.py`) + dsub (for
  `wgs_extract_variants.sh`); installed by `get_genotypes.sh` via
  `pip install -r`.
- `CLAUDE.md` — project conventions, AoU platform notes, portability rules.
- `reference/ukbb-sbayesrc-gwas/` — UKBB analog this pipeline mirrors.
- `reference/fileheads.txt` — heads of the AoU v8 input files (used to
  confirm column ordering of the `plink_bed` bim).
