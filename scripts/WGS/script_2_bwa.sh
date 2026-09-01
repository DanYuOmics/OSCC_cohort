#!/bin/bash

NUM_CORES=32
BWA_EXEC=/home/dyu/miniconda3/envs/tumour_env_clean/bin/bwa
BWA_INDEX=/data/yudan_ref/genomics/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fa
FASTQ_DIR=/data/OSCC2/fastq/WGS
BAM_DIR=/data/OSCC2/BAM

mkdir -p $BAM_DIR/tumour $BAM_DIR/normal

# 1. Tumour 样本比对 (23 对)
for fq1 in $FASTQ_DIR/tumour/*_1.fq.gz; do
    raw_name=$(basename "$fq1" | sed 's/_1.fq.gz//; s/_combined_1.fq.gz//')
    sample_id=$(echo "$raw_name" | awk -F'_' '{print $1}')
    [[ "$raw_name" =~ _M ]] && sample_id="${sample_id}_M"
    
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    fq2="${fq2/_1.combined.fq.gz/_2.combined.fq.gz}"
    
    out_bam="$BAM_DIR/tumour/${sample_id}.bam"
    if [ -f "$out_bam" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping existing: $out_bam"
        continue
    fi
    
    RG="@RG\tID:${sample_id}\tPU:Novogene\tSM:${sample_id}\tLB:WGS\tPL:ILLUMINA"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aligning tumour: $sample_id"
    $BWA_EXEC mem -t $NUM_CORES -M -R "$RG" "$BWA_INDEX" "$fq1" "$fq2" \
    | samtools view -Shb - > "$out_bam"
done

# 2. Normal 样本比对 (25 对)
for fq1 in $FASTQ_DIR/normal/*_1.fq.gz; do
    raw_name=$(basename "$fq1" | sed 's/_1.fq.gz//')
    sample_id=$(echo "$raw_name" | awk -F'_' '{print $1}')
    
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    out_bam="$BAM_DIR/normal/${sample_id}.bam"
    if [ -f "$out_bam" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping existing: $out_bam"
        continue
    fi
    
    RG="@RG\tID:${sample_id}\tPU:Novogene\tSM:${sample_id}\tLB:WGS\tPL:ILLUMINA"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aligning normal: $sample_id"
    $BWA_EXEC mem -t $NUM_CORES -M -R "$RG" "$BWA_INDEX" "$fq1" "$fq2" \
    | samtools view -Shb - > "$out_bam"
done

echo "All BWA alignments completed successfully."
