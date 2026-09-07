#!/bin/bash
set -euo pipefail

# ==============================================================================
# Pipeline Step 3: MarkDuplicates (Standard GATK Engine)
# Concurrency: 8 parallel workers, 24 GB heap per worker (~192 GB total RAM)
# Inputs:
#   /data/OSCC2/fastq/WGS/normal/bam/*.bam
#   /data/OSCC2/fastq/WGS/tumour/bam/*.bam
# Outputs:
#   /data/OSCC2/fastq/WGS/{normal,tumour}/bam_markdup/*.markdup.bam
# ==============================================================================

MAX_JOBS=8
JAVA_MEM="24G"

WGS_ROOT="/data/OSCC2/fastq/WGS"
GATK_EXEC="gatk"

process_markdup() {
    local in_bam="$1"
    local out_dir="$2"
    local sample_id
    sample_id=$(basename "${in_bam}" .bam)

    local out_bam="${out_dir}/${sample_id}.markdup.bam"
    local metrics="${out_dir}/metrics/${sample_id}.metrics.txt"
    local tmp_dir="${out_dir}/tmp/${sample_id}"

    if [ -f "${out_bam}" ] && [ -f "${out_bam}.bai" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SKIP] Deliverable exists: ${out_bam}"
        return 0
    fi

    mkdir -p "${tmp_dir}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [START] MarkDuplicates on ${sample_id}"

    "${GATK_EXEC}" --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${tmp_dir}" MarkDuplicates \
        -I "${in_bam}" \
        -O "${out_bam}" \
        -M "${metrics}" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT

    rm -rf "${tmp_dir}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FINISH] Completed: ${out_bam}"
}

export -f process_markdup
export GATK_EXEC JAVA_MEM

# ------------------------------------------------------------------------------
# 1. Normal Cohort
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing Normal cohort..."
mkdir -p "${WGS_ROOT}/normal/bam_markdup/metrics" "${WGS_ROOT}/normal/bam_markdup/tmp"

find "${WGS_ROOT}/normal/bam" -maxdepth 1 -name "*.bam" | sort | while read -r bam; do
    while [ "$(jobs -r -p | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
    process_markdup "${bam}" "${WGS_ROOT}/normal/bam_markdup" &
done
wait

# ------------------------------------------------------------------------------
# 2. Tumour Cohort
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing Tumour cohort..."
mkdir -p "${WGS_ROOT}/tumour/bam_markdup/metrics" "${WGS_ROOT}/tumour/bam_markdup/tmp"

find "${WGS_ROOT}/tumour/bam" -maxdepth 1 -name "*.bam" | sort | while read -r bam; do
    while [ "$(jobs -r -p | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
    process_markdup "${bam}" "${WGS_ROOT}/tumour/bam_markdup" &
done
wait

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All MarkDuplicates tasks completed."
