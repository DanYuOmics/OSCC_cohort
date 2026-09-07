#!/bin/bash

# 设置适度线程，避免与 BWA 争抢 I/O
QC_THREADS=6

FASTQ_DIR=/data/OSCC2/fastq/WGS
QC_OUT_DIR=/data/OSCC2/fastqc_out/WGS

mkdir -p "$QC_OUT_DIR"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting FastQC for WGS dataset..."

# 1. 正常对照组 (Normal) 质控
find "$FASTQ_DIR/normal" -name "*.fq.gz" | xargs -P "$QC_THREADS" -I {} fastqc -o "$QC_OUT_DIR" --noextract {}

# 2. 肿瘤组 (Tumour) 质控 (排除 raw_split_lanes 归档目录)
find "$FASTQ_DIR/tumour" -maxdepth 1 -name "*.fq.gz" | xargs -P "$QC_THREADS" -I {} fastqc -o "$QC_OUT_DIR" --noextract {}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All WGS FastQC jobs finished."

# 3. 自动生成全队列汇总报告 (如果安装了 multiqc)
if command -v multiqc &> /dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating MultiQC report..."
    multiqc "$QC_OUT_DIR" -o "$QC_OUT_DIR/multiqc_report"
fi
