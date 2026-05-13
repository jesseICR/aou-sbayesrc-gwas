# AoU SBayesRC GWAS pipeline

A reproducible, idempotent pipeline that builds per-chromosome PLINK2 pfiles
holding the ~7.35 M SBayesRC SNPs from the All of Us v8 ACAF WGS callset, with
the variant `ID` column populated with SBayesRC rsids. The AoU analog of the
early phase of the UKBB SBayesRC GWAS pipeline (UKBB steps 1–3: generate
variant ID list, upload, extract).

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
```

Each `summary.tsv` reports `requested / src_variants / src_samples /
biallelic_total / extracted / out_samples / missing / remapped / unmapped`
for that chromosome. A combined summary across all processed chromosomes is
also written to `logs/sbayesrc_extract_summary.tsv` at the end of each run.

Downstream UKBB-analog steps (kinship, ADMIXTURE, PCA, REGENIE) are not yet
ported — see `reference/ukbb-sbayesrc-gwas/get_genotypes.sh` for the full
target pipeline (gitignored locally; cloned for reference).

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
- **gcloud Application Default Credentials** (run `gcloud auth
  application-default login` once if `gcloud storage` or dsub return auth
  errors).
- **Pet service account** of the pod, used as the Batch worker identity
  (auto-detected via `gcloud config get-value account`).

Python dependencies (`pandas`, `dsub`) are pinned in `requirements.txt` and
auto-installed by `get_genotypes.sh`.

## Usage

Smoke-test a single chromosome end-to-end (one dsub task, ~15–25 min):
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
| `requirements.txt` | `pandas` (step 1), `dsub` (step 2 orchestration). |
| `CLAUDE.md` | Project conventions, AoU platform notes, portability rules, dsub-from-Jupyter recipe. Gitignored — local developer reference. |
| `reference/ukbb-sbayesrc-gwas/` | UKBB analog this pipeline mirrors. Gitignored — clone locally for reference only. |
