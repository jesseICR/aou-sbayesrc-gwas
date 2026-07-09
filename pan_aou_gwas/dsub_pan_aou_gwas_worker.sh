#!/bin/bash
# dsub_pan_aou_gwas_worker.sh - Run one pan-AoU HapMap3 PLINK2 GWAS batch.

set -euo pipefail

: "${PLINK2:?PLINK2 not set}"
: "${BED:?BED not set}"
: "${BIM:?BIM not set}"
: "${FAM:?FAM not set}"
: "${PHENO:?PHENO not set}"
: "${KEEP:?KEEP not set}"
: "${BATCH_MANIFEST:?BATCH_MANIFEST not set}"
: "${OUTDIR:?OUTDIR not set}"
: "${BATCH_INDEX:?BATCH_INDEX not set}"

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[pan-aou-gwas $(ts)] $*"; }

chmod +x "${PLINK2}"

SCRATCH="/mnt/data/scratch/pan_aou_gwas_${BATCH_INDEX}"
mkdir -p "${SCRATCH}" "${OUTDIR}"

BFILE="${SCRATCH}/hapmap3_bfile_hq"
ln -sf "${BED}" "${BFILE}.bed"
ln -sf "${BIM}" "${BFILE}.bim"
ln -sf "${FAM}" "${BFILE}.fam"

PREFIX="${SCRATCH}/batch_${BATCH_INDEX}"
PHENO_NAMES="$(awk 'NR > 1 {print $2}' "${BATCH_MANIFEST}" | paste -sd, -)"
if [[ -z "${PHENO_NAMES}" ]]; then
    log "ERROR: no phenotype names found in ${BATCH_MANIFEST}"
    exit 1
fi
DUP_PHENO_NAMES="$(awk 'NR > 1 {count[$2]++} END {for (name in count) if (count[name] > 1) print name}' "${BATCH_MANIFEST}" | sort | paste -sd, -)"
if [[ -n "${DUP_PHENO_NAMES}" ]]; then
    log "ERROR: duplicate phenotype names in batch manifest: ${DUP_PHENO_NAMES}"
    exit 1
fi

log "=== starting on $(hostname) ==="
log "batch_index=${BATCH_INDEX}"
log "phenotypes=$(awk 'NR > 1 {c++} END {print c+0}' "${BATCH_MANIFEST}")"
log "plink2=$("${PLINK2}" --version 2>&1 | head -1 || true)"
df -h /mnt/data | sed 's/^/  /'

"${PLINK2}" \
    --bfile "${BFILE}" \
    --keep "${KEEP}" \
    --pheno "${PHENO}" \
    --pheno-name "${PHENO_NAMES}" \
    --glm allow-no-covars cols=chrom,pos,a1freq,nobs,beta,se,p \
    --no-input-missing-phenotype \
    --threads "$(nproc)" \
    --out "${PREFIX}"

[[ -s "${PREFIX}.log" ]] || { log "ERROR: missing PLINK log ${PREFIX}.log"; exit 1; }

write_lightweight() {
    local glm="$1"
    local out="$2"
    awk -v OFS='\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                h[$i] = i
            }
            chrom = ("#CHROM" in h) ? h["#CHROM"] : h["CHROM"]
            pos = h["POS"]
            id = h["ID"]
            a1 = h["A1"]
            a1freq = ("A1_FREQ" in h) ? h["A1_FREQ"] : h["A1FREQ"]
            nobs = h["OBS_CT"]
            beta = h["BETA"]
            se = h["SE"]
            pcol = h["P"]
            print "rsid", "chrom", "pos", "allele1", "a1freq", "n", "beta", "se", "p", "log10p"
            next
        }
        {
            p = (pcol ? $pcol : "")
            log10p = ""
            if (p != "" && p != "NA" && p + 0 > 0) {
                log10p = -log(p + 0) / log(10)
            }
            print (id ? $id : ""), (chrom ? $chrom : ""), (pos ? $pos : ""), (a1 ? $a1 : ""), \
                  (a1freq ? $a1freq : ""), (nobs ? $nobs : ""), (beta ? $beta : ""), \
                  (se ? $se : ""), p, log10p
        }
    ' "${glm}" | gzip -c > "${out}"
}

tail -n +2 "${BATCH_MANIFEST}" | while IFS=$'\t' read -r pheno_id pheno_name glm_path sumstats_path gwas_params trait_type kind n n_cases n_controls covar_mode sex_filter extra_covariates construction_id; do
    [[ -n "${pheno_id}" ]] || continue
    local_glm="${PREFIX}.${pheno_name}.glm.linear"
    if [[ ! -s "${local_glm}" ]]; then
        log "ERROR: PLINK did not write expected output ${local_glm}"
        exit 1
    fi
    pheno_dir="${OUTDIR}/${pheno_id}"
    mkdir -p "${pheno_dir}"
    cp "${local_glm}" "${pheno_dir}/${pheno_id}.${pheno_name}.glm.linear"
    cp "${PREFIX}.log" "${pheno_dir}/${pheno_id}.plink2.log"
    write_lightweight "${local_glm}" "${pheno_dir}/${pheno_id}.sumstats.tsv.gz"
    {
        printf 'parameter\tvalue\n'
        printf 'pheno_id\t%s\n' "${pheno_id}"
        printf 'pheno_name\t%s\n' "${pheno_name}"
        printf 'trait_type\t%s\n' "${trait_type}"
        printf 'kind\t%s\n' "${kind}"
        printf 'n\t%s\n' "${n}"
        printf 'n_cases\t%s\n' "${n_cases}"
        printf 'n_controls\t%s\n' "${n_controls}"
        printf 'covar_mode\t%s\n' "${covar_mode:-full}"
        printf 'sex_filter\t%s\n' "${sex_filter:-all}"
        printf 'extra_covariates\t%s\n' "${extra_covariates:-}"
        printf 'construction_id\t%s\n' "${construction_id:-}"
    } > "${pheno_dir}/${pheno_id}.gwas.params.tsv"
done

log "=== done batch ${BATCH_INDEX} ==="
