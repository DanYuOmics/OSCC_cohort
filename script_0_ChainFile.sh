#!/bin/bash
set -euo pipefail

# ==============================================================================
# Script 0: Prepare Reference Sequences and Liftover Chain Files
# Purpose:  Download and index bi-directional chain files and reference FASTA
#           assets required for GRCh38 <-> T2T-CHM13v2.0 (hs1) CrossMap liftover.
# ==============================================================================

REF_BASE_DIR="/data/yudan_ref/genomics"
CHAIN_DIR="${REF_BASE_DIR}/liftover"
T2T_DIR="${REF_BASE_DIR}/t2t"
HG38_FA="${REF_BASE_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fa"

mkdir -p "${CHAIN_DIR}" "${T2T_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting environment and resource preparation..."

# ------------------------------------------------------------------------------
# 1. Download and decompress bi-directional chain files
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: Setting up UCSC chain files..."
cd "${CHAIN_DIR}"

if [ ! -f "hs1ToHg38.over.chain" ]; then
    echo "[INFO] Fetching hs1ToHg38.over.chain.gz..."
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz
    gzip -dc hs1ToHg38.over.chain.gz > hs1ToHg38.over.chain
else
    echo "[SKIP] hs1ToHg38.over.chain already present."
fi

if [ ! -f "hg38ToHs1.over.chain" ]; then
    echo "[INFO] Fetching hg38ToHs1.over.chain.gz..."
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz
    gzip -dc hg38ToHs1.over.chain.gz > hg38ToHs1.over.chain
else
    echo "[SKIP] hg38ToHs1.over.chain already present."
fi

# ------------------------------------------------------------------------------
# 2. Verify GRCh38 reference FASTA and index
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: Verifying GRCh38 FASTA index..."
if [ ! -f "${HG38_FA}" ]; then
    echo "[ERROR] GRCh38 reference FASTA not found at ${HG38_FA}" >&2
    exit 1
fi

if [ ! -f "${HG38_FA}.fai" ]; then
    echo "[INFO] Indexing GRCh38 reference FASTA with samtools..."
    samtools faidx "${HG38_FA}"
else
    echo "[SKIP] GRCh38 index (.fai) verified."
fi

# ------------------------------------------------------------------------------
# 3. Download and index T2T-CHM13v2.0 (hs1) reference genome
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: Setting up T2T-CHM13v2.0 (hs1) FASTA..."
cd "${T2T_DIR}"

if [ ! -f "hs1.fa" ]; then
    if [ ! -f "hs1.fa.gz" ]; then
        echo "[INFO] Downloading T2T hs1.fa.gz from UCSC..."
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
    fi
    echo "[INFO] Decompressing hs1.fa.gz..."
    gzip -d hs1.fa.gz
else
    echo "[SKIP] hs1.fa already exists."
fi

if [ ! -f "hs1.fa.fai" ]; then
    echo "[INFO] Indexing hs1.fa with samtools..."
    samtools faidx hs1.fa
else
    echo "[SKIP] hs1.fa index (.fai) verified."
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All liftover reference dependencies ready."
