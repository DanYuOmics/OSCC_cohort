#!/bin/bash
# ==============================================================================
# Script Name: script_6_extract_pass.sh
# Description: Batch extraction of high-confidence PASS somatic variants (SNVs & INDELs)
#              from Mutect2 filtered VCFs and generation of cohort statistics.
# ==============================================================================
set -euo pipefail

# 1. Directory and environment settings
PROJECT_DIR="/data/OSCC2/fastq/WGS"
INPUT_DIR="${PROJECT_DIR}/mutect2"
OUT_DIR="${INPUT_DIR}/pass_only_vcfs"
SUMMARY_TABLE="${OUT_DIR}/cohort_somatic_pass_stats.tsv"

mkdir -p "${OUT_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting PASS variant extraction..."
echo -e "Sample\tPASS_Total\tPASS_SNV\tPASS_INDEL" > "${SUMMARY_TABLE}"

# 2. Iterate through all completed filtered VCF files
for vcf in "${INPUT_DIR}"/*/*.filtered.vcf.gz; do
    sname=$(basename "${vcf}" .filtered.vcf.gz)
    out_vcf="${OUT_DIR}/${sname}.somatic.pass.vcf.gz"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${sname}..."
    
    # Extract strict PASS somatic variants and compress using bgzip
    bcftools view -f PASS -O z -o "${out_vcf}" "${vcf}"
    
    # Generate tabix index for random access
    bcftools index -t "${out_vcf}"
    
    # Calculate variant counts
    total_pass=$(bcftools view -H "${out_vcf}" | wc -l)
    snv_pass=$(bcftools view -v snps -H "${out_vcf}" | wc -l)
    indel_pass=$(bcftools view -v indels -H "${out_vcf}" | wc -l)
    
    echo -e "${sname}\t${total_pass}\t${snv_pass}\t${indel_pass}" >> "${SUMMARY_TABLE}"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extraction completed successfully."
echo "=== Cohort PASS Somatic Statistics ==="
column -t "${SUMMARY_TABLE}"
