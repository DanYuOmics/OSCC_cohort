#!/bin/bash
FASTQ_PATH="/data/OSCC2/fastq/RNA"
SAMPLE_LIST="${FASTQ_PATH}/samples.txt"
REF_PATH="/data/yudan_ref/RNA/grch38_snp_tran/genome_snp_tran"
OUT_DIR="/data/OSCC2/fastq/RNA/hisat2_out"

mkdir -p "${OUT_DIR}"

run_hisat() {
    sample_prefix=$1
    REF=$2
    OUT=$3
    
    # 提取纯样本名（去掉目录路径）
    sample_name=$(basename "${sample_prefix}")
    
    # 执行 HISAT2 比对并通过管道流直接由 samtools 排序并输出 BAM
    hisat2 -p 12 --dta -x "${REF}" \
      -1 "${sample_prefix}_1.fq.gz" \
      -2 "${sample_prefix}_2.fq.gz" | \
    samtools sort -@ 4 -m 3G -o "${OUT}/${sample_name}.bam" -

    samtools index -@ 4 "${OUT}/${sample_name}.bam"
}

export -f run_hisat

# -P 14 启动 14 个样本并行计算 (14 * 16 = 224 核心)
cat "${SAMPLE_LIST}" | xargs -n 1 -P 14 -I {} bash -c 'run_hisat "$@" '"${REF_PATH}"' '"${OUT_DIR}" _ {}
