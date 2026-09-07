#!/bin/bash
set -euo pipefail

# ==============================================================================
# Script: script_4_download_ref.sh
# Purpose: Download GATK Best Practices common variants resource for GetPileupSummaries
# Destination: /data/yudan_ref/genomics/
# ==============================================================================

REF_DIR="/data/yudan_ref/genomics"
BASE_URL="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"
VCF_FILE="small_exac_common_3.hg38.vcf.gz"
TBI_FILE="small_exac_common_3.hg38.vcf.gz.tbi"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting download for ${VCF_FILE} and ${TBI_FILE}..."

# Download VCF with resume support
curl -L -C - -O "${BASE_URL}/${VCF_FILE}"

# Download Tabix index with resume support
curl -L -C - -O "${BASE_URL}/${TBI_FILE}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Verifying file integrity..."
gzip -t "${VCF_FILE}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Download and verification completed successfully."
