#!/bin/bash
set -euo pipefail

# ==============================================================================
# Pipeline Step 2: High-Throughput Parallel BWA-MEM Alignment & Sorting
# Purpose: Execute 8-sample concurrent alignment and coordinate sorting.
# Destination paths:
#   Tumor:  /data/OSCC2/fastq/WGS/tumour/bam/${sample_id}.bam
#   Normal: /data/OSCC2/fastq/WGS/normal/bam/${sample_id}.bam
# Allocation per job: 16 BWA threads, 4 sort threads, 8 GB sort memory
# ==============================================================================

MAX_JOBS=8
BWA_THREADS=16
SORT_THREADS=4
SORT_MEM="8G"

BWA_EXEC="/home/dyu/miniconda3/envs/tumour_env_clean/bin/bwa"
BWA_INDEX="/data/yudan_ref/genomics/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fa"

WGS_ROOT="/data/OSCC2/fastq/WGS"
TUMOUR_DIR="${WGS_ROOT}/tumour"
NORMAL_DIR="${WGS_ROOT}/normal"

mkdir -p "${TUMOUR_DIR}/bam" "${TUMOUR_DIR}/logs" "${TUMOUR_DIR}/tmp"
mkdir -p "${NORMAL_DIR}/bam" "${NORMAL_DIR}/logs" "${NORMAL_DIR}/tmp"

run_bwa() {
    local fq1="$1"
    local sample_root="$2"

    local bam_out="${sample_root}/bam"
    local log_out="${sample_root}/logs"
    local tmp_dir="${sample_root}/tmp"

    local filename
    filename=$(basename "${fq1}")
    local raw_name
    raw_name=$(echo "${filename}" | sed -E 's/(_combined)?_1\.(fq|fastq)\.gz//')
    local sample_id
    sample_id=$(echo "${raw_name}" | awk -F'_' '{print $1}')
    [[ "${raw_name}" =~ _M ]] && sample_id="${sample_id}_M"

    local fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    fq2="${fq2/_1.combined.fq.gz/_2.combined.fq.gz}"
    fq2="${fq2/_combined_1.fq.gz/_combined_2.fq.gz}"

    local final_bam="${bam_out}/${sample_id}.bam"
    local tmp_prefix="${tmp_dir}/${sample_id}.sort"

    if [ -f "${final_bam}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SKIP] Deliverable exists: ${final_bam}"
        return 0
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [START] Processing sample: ${sample_id}"

    local rg="@RG\tID:${sample_id}\tPU:Novogene\tSM:${sample_id}\tLB:WGS\tPL:ILLUMINA"

    "${BWA_EXEC}" mem -t "${BWA_THREADS}" -M -R "${rg}" "${BWA_INDEX}" "${fq1}" "${fq2}" 2> "${log_out}/${sample_id}.log" | \
    samtools sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -T "${tmp_prefix}" -o "${final_bam}" -

    samtools index -@ "${SORT_THREADS}" "${final_bam}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FINISH] Generated: ${final_bam}"
}

export -f run_bwa
export BWA_EXEC BWA_INDEX BWA_THREADS SORT_THREADS SORT_MEM

# ------------------------------------------------------------------------------
# 1. Tumor Cohort Alignment Queue (8 concurrent workers)
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Initializing Tumor alignment queue..."
find "${TUMOUR_DIR}" -maxdepth 1 -name "*_1.fq.gz" | sort | while read -r fq1; do
    while [ "$(jobs -r -p | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
    run_bwa "${fq1}" "${TUMOUR_DIR}" &
done
wait

# ------------------------------------------------------------------------------
# 2. Normal Cohort Alignment Queue (8 concurrent workers)
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Initializing Normal alignment queue..."
find "${NORMAL_DIR}" -maxdepth 1 -name "*_1.fq.gz" | sort | while read -r fq1; do
    while [ "$(jobs -r -p | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
    run_bwa "${fq1}" "${NORMAL_DIR}" &
done
wait

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All alignments finished."
