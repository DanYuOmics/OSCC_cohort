#!/bin/bash
BAM_DIR="/data/OSCC2/fastq/RNA/hisat2_out"
REF_GTF="/data/yudan_ref/RNA/genome.gtf"
OUT_DIR="/data/OSCC2/fastq/RNA/stringtie_assemble"

mkdir -p "${OUT_DIR}"

run_stringtie() {
    bam_file="$1"
    ref_gtf="$2"
    out_dir="$3"
    
    sample=$(basename "$bam_file" .bam)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Assembling transcript for sample: ${sample}"
    
    stringtie "$bam_file" \
        -G "$ref_gtf" \
        -l "${sample}" \
        -o "${out_dir}/${sample}.gtf" \
        -p 8
}

export -f run_stringtie

# 并行执行：7 个任务同时跑，每个任务 8 线程
ls -1 ${BAM_DIR}/*.bam | xargs -n 1 -P 7 -I {} bash -c 'run_stringtie "$@"' _ {} "${REF_GTF}" "${OUT_DIR}"

echo "All 48 samples transcript assembly completed."
