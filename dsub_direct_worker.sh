#!/bin/bash
# dsub_direct_worker.sh — Build one chromosome's present direct-SNP pfile on Batch.

set -euo pipefail

: "${CHROM:?CHROM not set}"
: "${DESIRED:?DESIRED not set}"
: "${AVAILABLE:?AVAILABLE not set}"
: "${MISSING:?MISSING not set}"
: "${PLINK2:?PLINK2 not set}"
: "${PGEN:?PGEN not set}"
: "${PVAR:?PVAR not set}"
: "${PSAM:?PSAM not set}"
: "${EXTRACT:?EXTRACT not set}"
: "${OUTDIR:?OUTDIR not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[${CHROM} $(ts)] $*"; }
free_data() { df -BG --output=avail /mnt/data | tail -1 | tr -d ' G'; }

chmod +x "${PLINK2}"

mkdir -p "${OUTDIR}"
SCRATCH=/mnt/data/scratch
mkdir -p "${SCRATCH}"

SRC_PREFIX="${PGEN%.pgen}"
if [[ ! -f "${SRC_PREFIX}.pvar" || ! -f "${SRC_PREFIX}.psam" ]]; then
    log "ERROR: pgen/pvar/psam don't share local prefix '${SRC_PREFIX}'"
    ls -la "${PGEN}" "${PVAR}" "${PSAM}" 2>&1 | sed "s/^/  /"
    exit 1
fi
if [[ ! -f "${PVAR}" || ! -f "${PSAM}" ]]; then
    log "ERROR: missing localized pvar/psam inputs"
    ls -la "${PGEN}" "${PVAR}" "${PSAM}" 2>&1 | sed "s/^/  /"
    exit 1
fi

FINAL_PREFIX="${OUTDIR}/${CHROM}"
PRESENT_PREFIX="${SCRATCH}/${CHROM}.present"

sample_count=$(grep -vc '^#' "${SRC_PREFIX}.psam")
available=$(grep -vc '^#' "${EXTRACT}" || true)
if [[ "${available}" -ne "${AVAILABLE}" ]]; then
    log "ERROR: extract list has ${available} variants, expected AVAILABLE=${AVAILABLE}"
    exit 1
fi
desired="${DESIRED}"
missing="${MISSING}"

log "=== starting on $(hostname) ==="
log "plink2 = $("${PLINK2}" --version 2>&1 | head -1 || true)"
log "samples=${sample_count} desired=${desired} available=${available} missing=${missing} free-data=$(free_data)G"

if [[ "${available}" -gt 0 ]]; then
    log "extracting present direct SNPs ..."
    "${PLINK2}" \
        --pfile "${SRC_PREFIX}" \
        --extract "${EXTRACT}" \
        --make-pgen \
        --sort-vars \
        --no-pheno \
        --threads "$(nproc)" \
        --out "${PRESENT_PREFIX}"
fi

if [[ "${available}" -le 0 ]]; then
    log "ERROR: no direct SNPs requested for ${CHROM}"
    exit 1
fi

final_count=$(grep -vc '^#' "${PRESENT_PREFIX}.pvar")
if [[ "${final_count}" -ne "${available}" ]]; then
    log "ERROR: final direct pfile has ${final_count} variants, expected available=${available}"
    exit 1
fi

mv -f "${PRESENT_PREFIX}.pgen" "${FINAL_PREFIX}.pgen"
mv -f "${PRESENT_PREFIX}.pvar" "${FINAL_PREFIX}.pvar"
mv -f "${PRESENT_PREFIX}.psam" "${FINAL_PREFIX}.psam"
[[ -f "${PRESENT_PREFIX}.log" ]] && mv -f "${PRESENT_PREFIX}.log" "${FINAL_PREFIX}.log"

{
    printf 'chrom\tdesired\tavailable_in_wgs_pfiles\tmissing_from_wgs_pfiles\tfinal_variants\tsamples\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${CHROM}" "${desired}" "${available}" "${missing}" "${final_count}" "${sample_count}"
} > "${OUTDIR}/${CHROM}.summary.tsv"

log "=== done ==="
ls -lh "${OUTDIR}" | sed "s/^/  /"
