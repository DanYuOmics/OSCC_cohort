#!/bin/bash
set -euo pipefail

# ==============================================================================
# Pipeline Step 5: High-Throughput Somatic Variant Calling via GATK Mutect2
# Strategy:
#   1. Mutect2 Somatic Calling (Tumour-Normal matched mode with PoN & gnomAD)
#   2. Read Orientation Artifact Modeling (LearnReadOrientationModel)
#   3. Pileup Summaries over Common Biallelic Sites (GetPileupSummaries)
#   4. Contamination & Segment Estimation (CalculateContamination)
#   5. High-Confidence Somatic Filtering (FilterMutectCalls)
# Hardware:
#   4 concurrent pairs, 16 native PairHMM threads per pair
# ==============================================================================

MAX_JOBS=4
HMM_THREADS=16

# Reference datasets and resource bundles
REF_FA="/data/yudan_ref/genomics/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fa"
PON="/data/yudan_ref/genomics/1000g_pon.hg38.vcf.gz"
GERM="/data/yudan_ref/genomics/af-only-gnomad.hg38.vcf.gz"
COMMON_SNPS="/data/yudan_ref/genomics/small_exac_common_3.hg38.vcf.gz"

# Working directories
WGS_ROOT="/data/OSCC2/fastq/WGS"
TUMOUR_BAM_DIR="${WGS_ROOT}/tumour/bam_markdup"
NORMAL_BAM_DIR="${WGS_ROOT}/normal/bam_markdup"

OUT_ROOT="${WGS_ROOT}/mutect2"
LOG_DIR="${OUT_ROOT}/logs"
TMP_DIR="${OUT_ROOT}/tmp"

mkdir -p "${OUT_ROOT}" "${LOG_DIR}" "${TMP_DIR}"

run_paired_mutect2() {
    local tumour_bam="$1"
    local sample_id
    sample_id=$(basename "${tumour_bam}" .markdup.bam)

    # Map matched normal sample ID: DNA{N}_M -> DNA{N}, S{N}TB -> S{N}BA
    local normal_id=""
    if [[ "${sample_id}" =~ _M$ ]]; then
        normal_id="${sample_id%_M}"
    elif [[ "${sample_id}" =~ TB$ ]]; then
        normal_id="${sample_id%TB}BA"
    fi

    local normal_bam="${NORMAL_BAM_DIR}/${normal_id}.markdup.bam"

    # Verify presence of matched normal BAM
    if [ ! -f "${normal_bam}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] Matched Normal BAM missing for ${sample_id}: ${normal_bam}"
        return 1
    fi

    local sample_out="${OUT_ROOT}/${sample_id}"
    local sample_tmp="${TMP_DIR}/${sample_id}"
    local sample_log="${LOG_DIR}/${sample_id}.log"
    mkdir -p "${sample_out}" "${sample_tmp}"

    local final_filtered_vcf="${sample_out}/${sample_id}.filtered.vcf.gz"
    if [ -f "${final_filtered_vcf}" ] && [ -f "${final_filtered_vcf}.tbi" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SKIP] Deliverable already exists: ${final_filtered_vcf}"
        return 0
    fi

    {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [START] Processing pair: Tumour=${sample_id} | Normal=${normal_id}"

        # Extract sample names (SM tags) from BAM headers
        local tumour_sm
        local normal_sm
        tumour_sm=$(samtools view -H "${tumour_bam}" | awk -F'\t' '/^@RG/ {for(i=1;i<=NF;i++) if ($i ~ /^SM:/){split($i,a,":"); print a[2]; exit}}')
        normal_sm=$(samtools view -H "${normal_bam}" | awk -F'\t' '/^@RG/ {for(i=1;i<=NF;i++) if ($i ~ /^SM:/){split($i,a,":"); print a[2]; exit}}')

        # 1. Mutect2 Somatic Calling
        gatk --java-options "-Xmx24G -Djava.io.tmpdir=${sample_tmp}" Mutect2 \
            -R "${REF_FA}" \
            -I "${tumour_bam}" -tumor "${tumour_sm}" \
            -I "${normal_bam}" -normal "${normal_sm}" \
            --panel-of-normals "${PON}" \
            --germline-resource "${GERM}" \
            --f1r2-tar-gz "${sample_out}/${sample_id}.f1r2.tar.gz" \
            --native-pair-hmm-threads "${HMM_THREADS}" \
            --tmp-dir "${sample_tmp}" \
            -O "${sample_out}/${sample_id}.unfiltered.vcf.gz"

        # 2. Learn Read Orientation Model
        gatk --java-options "-Xmx16G -Djava.io.tmpdir=${sample_tmp}" LearnReadOrientationModel \
            -I "${sample_out}/${sample_id}.f1r2.tar.gz" \
            -O "${sample_out}/${sample_id}.read-orientation-model.tar.gz"

        # 3. Pileup Summaries (Tumour & Normal)
        gatk --java-options "-Xmx16G -Djava.io.tmpdir=${sample_tmp}" GetPileupSummaries \
            -I "${tumour_bam}" \
            -V "${COMMON_SNPS}" \
            -L "${COMMON_SNPS}" \
            -R "${REF_FA}" \
            -O "${sample_out}/${sample_id}.tumour.pileups.table"

        gatk --java-options "-Xmx16G -Djava.io.tmpdir=${sample_tmp}" GetPileupSummaries \
            -I "${normal_bam}" \
            -V "${COMMON_SNPS}" \
            -L "${COMMON_SNPS}" \
            -R "${REF_FA}" \
            -O "${sample_out}/${sample_id}.normal.pileups.table"

        # 4. Calculate Contamination
        gatk --java-options "-Xmx16G -Djava.io.tmpdir=${sample_tmp}" CalculateContamination \
            -I "${sample_out}/${sample_id}.tumour.pileups.table" \
            -matched "${sample_out}/${sample_id}.normal.pileups.table" \
            -O "${sample_out}/${sample_id}.contamination.table" \
            --tumor-segmentation "${sample_out}/${sample_id}.segments.table"

        # 5. Filter Mutect Calls
        gatk --java-options "-Xmx16G -Djava.io.tmpdir=${sample_tmp}" FilterMutectCalls \
            -R "${REF_FA}" \
            -V "${sample_out}/${sample_id}.unfiltered.vcf.gz" \
            --contamination-table "${sample_out}/${sample_id}.contamination.table" \
            --tumor-segmentation "${sample_out}/${sample_id}.segments.table" \
            --ob-priors "${sample_out}/${sample_id}.read-orientation-model.tar.gz" \
            -O "${final_filtered_vcf}"

        rm -rf "${sample_tmp}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FINISH] Completed: ${final_filtered_vcf}"

    } > "${sample_log}" 2>&1
}

export -f run_paired_mutect2
export MAX_JOBS HMM_THREADS REF_FA PON GERM COMMON_SNPS WGS_ROOT TUMOUR_BAM_DIR NORMAL_BAM_DIR OUT_ROOT LOG_DIR TMP_DIR

# ------------------------------------------------------------------------------
# Launch Parallel Somatic Calling Pool
# ------------------------------------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Launching Mutect2 Somatic Calling queue..."

find "${TUMOUR_BAM_DIR}" -maxdepth 1 -name "*.markdup.bam" | sort | while read -r tb; do
    while [ "$(jobs -r -p | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 10
    done
    run_paired_mutect2 "${tb}" &
done
wait

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All Mutect2 Somatic variant calling tasks completed."
