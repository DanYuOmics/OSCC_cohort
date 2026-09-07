#!/bin/bash
BAM_DIR="/data/OSCC2/fastq/RNA/hisat2_out"
MERGED_GTF="/data/OSCC2/fastq/RNA/stringtie_merged.gtf"
QUANT_OUT="/data/OSCC2/fastq/RNA/ballgown_out"

mkdir -p "${QUANT_OUT}"

quant_sample() {
    bam_file="$1"
    merged_gtf="$2"
    quant_out="$3"
    
    sample=$(basename "$bam_file" .bam)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quantifying sample: ${sample}"
    
    mkdir -p "${quant_out}/${sample}"
    
    # -e: 仅对 merged.gtf 中的已知/组装转录本定量
    # -B: 输出 Ballgown 所需的 5 个 ctab 文件
    # -A: 输出基因级别的丰度估算表 (gene_abundances.tsv)
    stringtie "$bam_file" \
        -e -B \
        -p 8 \
        -G "$merged_gtf" \
        -o "${quant_out}/${sample}/${sample}.gtf" \
        -A "${quant_out}/${sample}/gene_abundances.tsv"
}

export -f quant_sample

# 7 进程并行，每个样本分配 8 线程
ls -1 ${BAM_DIR}/*.bam | xargs -n 1 -P 7 -I {} bash -c 'quant_sample "$@"' _ {} "${MERGED_GTF}" "${QUANT_OUT}"

echo "All samples quantified successfully."
