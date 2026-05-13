#!/bin/bash
# dsub_extract_worker.sh — Runs inside one Batch worker for one chromosome.
#
# Submitted by wgs_extract_variants.sh as a dsub task (one task per chrom).
# Same four passes as the previous serial-on-Jupyter version of the extract:
#
#   Pass 1   — plink2 split multi-allelics + assign chr:pos:REF:ALT IDs
#   Pass 1.5 — right-trim shared REF/ALT suffix in the split pvar so SNPs
#              co-located with indels match parsimonious SBayesRC IDs.
#              Logic is identical to normalize_pvar_alleles.py (kept in the
#              repo as the canonical, testable reference); here it's inlined
#              in awk so the worker doesn't need python3 (the default
#              marketplace.gcr.io/google/ubuntu2204 image ships without it,
#              and VPC SC blocks installing it).
#   Pass 2   — plink2 --extract SBayesRC variants, preserve INFO
#   Pass 3   — awk: remap pvar ID from chr:pos:REF:ALT to rsid via idmap
#
# Inputs (env vars set by dsub --input / --env):
#   CHROM         e.g. "chr3"
#   PLINK2        local path to statically-linked plink2 binary
#   PGEN PVAR PSAM   local paths to AoU acaf_threshold.chrN.{pgen,pvar,psam}
#   EXTRACT       local path to chrN.extract.txt (one ID per line)
#   IDMAP         local path to chrN.idmap.txt (chr:pos:REF:ALT \t rsid)
#
# Output (dsub --output-recursive OUTDIR):
#   ${OUTDIR}/${CHROM}.pgen          — final extracted pgen
#   ${OUTDIR}/${CHROM}.pvar          — final pvar with rsids
#   ${OUTDIR}/${CHROM}.psam          — sample file (unchanged from source)
#   ${OUTDIR}/${CHROM}.summary.tsv   — one-row TSV with key counts

set -euo pipefail

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[${CHROM} $(ts)] $*"; }
free_root() { df -BG --output=avail / | tail -1 | tr -d ' G'; }

# Make plink2 executable BEFORE anything that tries to invoke it.
chmod +x "${PLINK2}"

log "=== starting on $(hostname) ==="
log "cpu=$(nproc) mem=$(awk '/MemTotal/{printf "%.1f", $2/1024/1024}' /proc/meminfo)G disk-root=$(free_root)G"
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"

mkdir -p "${OUTDIR}"

# pgen/pvar/psam all come from the same gs:// dir so they share a local prefix
# after dsub localization. Derive it from PGEN.
SRC_PREFIX="${PGEN%.pgen}"
if [[ ! -f "${SRC_PREFIX}.pvar" || ! -f "${SRC_PREFIX}.psam" ]]; then
    log "ERROR: pgen/pvar/psam don't share local prefix '${SRC_PREFIX}'"
    ls -la "${PGEN}" "${PVAR}" "${PSAM}" 2>&1 | sed "s/^/  /"
    exit 1
fi

# Scratch on the data disk (where dsub mounts /mnt/data).
SCRATCH=/mnt/data/scratch
mkdir -p "${SCRATCH}"
SPLIT_PREFIX="${SCRATCH}/${CHROM}.split"
FINAL_PREFIX="${OUTDIR}/${CHROM}"

# Source counts
REQUESTED=$(wc -l < "${EXTRACT}")
SRC_VARIANTS=$(grep -vc '^#' "${SRC_PREFIX}.pvar")
SRC_SAMPLES=$(grep -vc '^#' "${SRC_PREFIX}.psam")
log "requested SBayesRC variants:      ${REQUESTED}"
log "AoU pgen source variants (multi): ${SRC_VARIANTS}"
log "AoU pgen source samples:          ${SRC_SAMPLES}"

# --------- Pass 1: split + ID assign ----------------------------------------
log "pass 1: split multi + assign IDs ..."
"${PLINK2}" \
    --pfile "${SRC_PREFIX}" \
    --output-chr chrM \
    --make-pgen 'multiallelics=-' \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 10000 \
    --no-pheno \
    --threads "$(nproc)" \
    --out "${SPLIT_PREFIX}"
BIALLELIC_TOTAL=$(grep -vc '^#' "${SPLIT_PREFIX}.pvar")
log "post-split biallelic variants:    ${BIALLELIC_TOTAL}"

# Free input pgen now that split is materialized — saves ~80 GB on chr1.
rm -f "${SRC_PREFIX}.pgen" "${SRC_PREFIX}.pvar" "${SRC_PREFIX}.psam"
log "free / after input cleanup: $(free_root)G"

# --------- Pass 1.5: normalize pvar alleles (right-trim shared suffix) ------
# Mirrors normalize_pvar_alleles.py: while ref and alt both have >1 base and
# end in the same base, drop the trailing base. If the (REF, ALT) pair changes,
# rebuild the ID column as chr:pos:REF:ALT. Genotypes don't move — pgen indexes
# by row order with allele codes (0=REF, 1=ALT), only the text labels change.
log "pass 1.5: right-trim split pvar alleles (inline awk) ..."
awk -v OFS='\t' '
    BEGIN { FS = "\t"; total = 0; rewritten = 0 }
    /^#/  { print; next }
    NF < 5 { print; next }
    {
        ref = $4; alt = $5
        rl = length(ref); al = length(alt)
        while (rl > 1 && al > 1 && substr(ref, rl, 1) == substr(alt, al, 1)) {
            ref = substr(ref, 1, rl - 1)
            alt = substr(alt, 1, al - 1)
            rl--; al--
        }
        if (ref != $4 || alt != $5) {
            $4 = ref; $5 = alt
            $3 = $1 ":" $2 ":" ref ":" alt
            rewritten++
        }
        total++
        print
    }
    END {
        printf "[normalize_pvar_alleles: total=%d rewritten=%d (%.2f%%)]\n", \
            total, rewritten, (total > 0 ? 100.0 * rewritten / total : 0.0) > "/dev/stderr"
    }
' "${SPLIT_PREFIX}.pvar" > "${SPLIT_PREFIX}.pvar.tmp"
mv -f "${SPLIT_PREFIX}.pvar.tmp" "${SPLIT_PREFIX}.pvar"

# --------- Pass 2: extract SBayesRC subset ----------------------------------
log "pass 2: extract SBayesRC subset ..."
"${PLINK2}" \
    --pfile "${SPLIT_PREFIX}" \
    --extract "${EXTRACT}" \
    --no-pheno \
    --output-chr chrM \
    --make-pgen \
    --threads "$(nproc)" \
    --out "${FINAL_PREFIX}"

# Free split intermediate (~115 GB on chr1)
rm -f "${SPLIT_PREFIX}".{pgen,pvar,psam,log}
log "free / after split cleanup: $(free_root)G"

# --------- Pass 3: rsid remap -----------------------------------------------
log "pass 3: rsid remap in pvar ..."
awk -v IDMAP="${IDMAP}" '
    BEGIN {
        FS = OFS = "\t"
        while ((getline line < IDMAP) > 0) {
            split(line, a, "\t"); m[a[1]] = a[2]
        }
        close(IDMAP)
    }
    /^#/ { print; next }
    { if ($3 in m) $3 = m[$3]; print }
' "${FINAL_PREFIX}.pvar" > "${FINAL_PREFIX}.pvar.tmp"
mv -f "${FINAL_PREFIX}.pvar.tmp" "${FINAL_PREFIX}.pvar"

# --------- Counts + summary -------------------------------------------------
EXTRACTED=$(grep -vc '^#' "${FINAL_PREFIX}.pvar")
OUT_SAMPLES=$(grep -vc '^#' "${FINAL_PREFIX}.psam")
MISSING=$(( REQUESTED - EXTRACTED ))
REMAPPED=$(awk -F'\t' '!/^#/ { if ($3 !~ /^chr[0-9XYM]+:[0-9]+:/) c++ } END { print c+0 }' "${FINAL_PREFIX}.pvar")
UNMAPPED=$(( EXTRACTED - REMAPPED ))

log "extracted variants:                 ${EXTRACTED}"
log "output samples:                     ${OUT_SAMPLES}"
log "SBayesRC variants not found in AoU: ${MISSING}"
log "pvar IDs remapped to rsid:          ${REMAPPED}"
log "pvar IDs left in chr:pos form:      ${UNMAPPED}"

{
    printf 'chrom\trequested\tsrc_variants\tsrc_samples\tbiallelic_total\textracted\tout_samples\tmissing\tremapped\tunmapped\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${CHROM}" "${REQUESTED}" "${SRC_VARIANTS}" "${SRC_SAMPLES}" \
        "${BIALLELIC_TOTAL}" "${EXTRACTED}" "${OUT_SAMPLES}" \
        "${MISSING}" "${REMAPPED}" "${UNMAPPED}"
} > "${OUTDIR}/${CHROM}.summary.tsv"

# Drop plink2's progress-bar log; we have our own structured stdout.
rm -f "${FINAL_PREFIX}.log"

log "=== done ==="
ls -lh "${OUTDIR}" | sed "s/^/  /"
