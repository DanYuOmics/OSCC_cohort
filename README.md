# OSCC Multi-Omics Sequencing Pipeline (WGS & RNA-Seq)

This repository contains the end-to-end data processing pipelines, quality control protocols, alignment workflows, and somatic variant calling procedures for the Oral Squamous Cell Carcinoma (OSCC) multi-omics cohort. The dataset integrates whole-genome sequencing (WGS) and bulk RNA-sequencing (RNA-Seq) across two independent sequencing batches.

---

## 1. Directory Structure and Workspace Organization

All raw FASTQ files, intermediate alignments, somatic variant calls, and reference resources are structured under `/data/OSCC2/` and reference storage `/data/yudan_ref/`:

```text
/data/OSCC2/
├── fastq/
│   ├── RNA/                                # Bulk RNA-Seq module
│   │   ├── samples.txt                     # RNA-Seq sample prefix list (48 samples)
│   │   ├── script_1_fastqc.sh              # FastQC execution script
│   │   ├── script_2_hisat.sh               # HISAT2 alignment & coordinate sorting
│   │   ├── script_3_stringtie_assemble.sh  # Transcript assembly per sample
│   │   ├── script_4_merge.sh               # Transcriptome merge script
│   │   ├── script_5_stringtie_quant.sh     # Abundance estimation (-e -B)
│   │   ├── script_6_extraction.sh          # Expression matrix extraction (prepDE.py)
│   │   ├── hisat2_out/                     # Sorted alignments (*.sort.bam, *.sort.bam.bai)
│   │   ├── stringtie_assemble/             # Sample-level GTFs
│   │   ├── stringtie_merged.gtf            # Non-redundant transcriptome GTF
│   │   ├── ballgown_out/                   # Ballgown quantification tables (*.ctab)
│   │   ├── gene_count_matrix.csv           # Gene-level raw count matrix
│   │   └── transcript_count_matrix.csv     # Transcript-level raw count matrix
│   │
│   └── WGS/                                # WGS raw data, alignments & Mutect2 workspace
│       ├── normal/
│       │   ├── *.fq.gz                     # Germline / normal FASTQs (25 samples)
│       │   └── bam_markdup/                # Deduplicated normal BAMs (*.markdup.bam, *.bai)
│       ├── tumour/
│       │   ├── *.fq.gz                     # Primary tumor FASTQs (23 pairs)
│       │   ├── raw_split_lanes/            # Archived multi-lane FASTQs (S19, S21, S24, S25)
│       │   └── bam_markdup/                # Deduplicated tumor BAMs (*.markdup.bam, *.bai)
│       ├── mutect2/                        # GATK4 somatic variant calling workspace
│       │   ├── <Sample_ID>/                # Sample-level Mutect2 raw calls & models
│       │   │   ├── *.unfiltered.vcf.gz     # Raw somatic candidate calls
│       │   │   ├── *.filtered.vcf.gz       # FilterMutectCalls output
│       │   │   ├── *.f1r2.tar.gz           # Read orientation model data
│       │   │   └── *.contamination.table   # Cross-sample contamination metrics
│       │   ├── pass_only_vcfs/             # High-confidence somatic deliverables
│       │   │   ├── *.somatic.pass.vcf.gz   # Tabix-indexed PASS-only somatic VCFs (23 pairs)
│       │   │   └── cohort_somatic_pass_stats.tsv # Cohort somatic SNV/INDEL statistics
│       │   └── logs/                       # Mutect2 execution logs
│       └── scripts/WGS/                    # Pipeline scripts (tracked in Git)
│           ├── script_0_ChainFile.sh       # Resource verification & chain file setup
│           ├── script_1_merge_lane.sh      # Dual-lane FASTQ merging for top-up samples
│           ├── script_1b_qc.sh             # FastQC quality control on raw reads
│           ├── script_2_bwa.sh             # BWA-MEM batch alignment and BAM generation
│           ├── script_3_markdup.sh         # Picard MarkDuplicates & coordinate indexing
│           ├── script_4_download_ref.sh    # Reference genome & GATK resource retrieval
│           ├── script_5_mutect2.sh         # Mutect2 somatic calling & filtering engine
│           └── script_6_extract_pass.sh    # PASS-only extraction & cohort mutation summary
│
├── OSCC_cohort/                            # Git tracking repository (DanYuOmics/OSCC_cohort)
│   └── scripts/
│       ├── RNA/                            # Archived RNA-Seq scripts
│       └── WGS/                            # Synchronized WGS pipeline scripts
└── fastqc_out/                             # Global FastQC reports

/data/yudan_ref/
├── RNA/
│   └── grch38_snp_tran/                    # HISAT2 GRCh38 SNP-aware reference indices
│       ├── genome_snp_tran.1.ht2 ... .8.ht2
│       └── Homo_sapiens.GRCh38.109.gtf     # Ensembl annotation GTF
└── genomics/
    ├── Homo_sapiens_assembly38/
    │   ├── Homo_sapiens_assembly38.fa      # GRCh38 reference FASTA
    │   ├── Homo_sapiens_assembly38.fa.fai  # Reference FASTA index
    │   ├── Homo_sapiens_assembly38.dict    # Sequence dictionary
    │   └── Homo_sapiens_assembly38.fa.*    # BWA indices (.bwt, .sa, .pac, .ann, .amb)
    ├── 1000g_pon.hg38.vcf.gz               # Mutect2 Panel of Normals (PoN)
    └── af-only-gnomad.hg38.vcf.gz          # gnomAD germline allele frequency resource
```

---

## 2. Sample Matrix & Multi-Omics Availability

The cohort comprises **25 individual patients** profiled across two sequencing batches:
- **Batch 1 (Patients 1–12):** 10 complete pairs with both WGS and RNA-Seq, alongside 2 Normal-only controls (Patients 4 and 6).
- **Batch 2 (Patients 13–26, excluding 17):** 13 complete pairs with full WGS and RNA-Seq coverage.

### Batch 1 (10 Complete Pairs + 2 Normal Controls)

| Patient ID | WGS Tumor ID | WGS Normal ID | RNA Tumor ID | RNA Normal ID | Data Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Patient 1** | `DNA1_M` | `DNA1` | `STRNA1` | `SPRNA1` | Complete (WGS + RNA) |
| **Patient 2** | `DNA2_M` | `DNA2` | `STRNA2` | `SPRNA2` | Complete (WGS + RNA) |
| **Patient 3** | `DNA3_M` | `DNA3` | `STRNA3` | `SPRNA3` | Complete (WGS + RNA) |
| **Patient 4** | *Missing* | `DNA4` | *Missing* | `SPRNA4` | Normal Only |
| **Patient 5** | `DNA5_M` | `DNA5` | `STRNA5` | `SPRNA5` | Complete (WGS + RNA) |
| **Patient 6** | *Missing* | `DNA6` | *Missing* | `SPRNA6` | Normal Only |
| **Patient 7** | `DNA7_M` | `DNA7` | `STRNA7` | `SPRNA7` | Complete (WGS + RNA) |
| **Patient 8** | `DNA8_M` | `DNA8` | `STRNA8` | `SPRNA8` | Complete (WGS + RNA) |
| **Patient 9** | `DNA9_M` | `DNA9` | `STRNA9` | `SPRNA9` | Complete (WGS + RNA) |
| **Patient 10** | `DNA10_M` | `DNA10` | `STRNA10` | `SPRNA10` | Complete (WGS + RNA) |
| **Patient 11** | `DNA11_M` | `DNA11` | `STRNA11` | `SPRNA11` | Complete (WGS + RNA) |
| **Patient 12** | `DNA12_M` | `DNA12` | `STRNA12` | `SPRNA12` | Complete (WGS + RNA) |

### Batch 2 (13 Complete Pairs)

| Patient ID | WGS Tumor ID | WGS Normal ID | RNA Tumor ID | RNA Normal ID | Data Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Patient 13** | `S13TB` | `S13BA` | `S13TA` | `S13NA` | Complete (WGS + RNA) |
| **Patient 14** | `S14TB` | `S14BA` | `S14TA` | `S14NA` | Complete (WGS + RNA) |
| **Patient 15** | `S15TB` | `S15BA` | `S15TA` | `S15NA` | Complete (WGS + RNA) |
| **Patient 16** | `S16TB` | `S16BA` | `S16TA` | `S16NA` | Complete (WGS + RNA) |
| **Patient 18** | `S18TB` | `S18BA` | `S18TA` | `S18NA` | Complete (WGS + RNA) |
| **Patient 19** | `S19TB` | `S19BA` | `S19TA` | `S19NA` | Complete (WGS + RNA) |
| **Patient 20** | `S20TB` | `S20BA` | `S20TA` | `S20NA` | Complete (WGS + RNA) |
| **Patient 21** | `S21TB` | `S21BA` | `S21TA` | `S21NA` | Complete (WGS + RNA) |
| **Patient 22** | `S22TB` | `S22BA` | `S22TA` | `S22NA` | Complete (WGS + RNA) |
| **Patient 23** | `S23TB` | `S23BA` | `S23TA` | `S23NA` | Complete (WGS + RNA) |
| **Patient 24** | `S24TB` | `S24BA` | `S24TA` | `S24NA` | Complete (WGS + RNA) |
| **Patient 25** | `S25TB` | `S25BA` | `S25TA` | `S25NA` | Complete (WGS + RNA) |
| **Patient 26** | `S26TB` | `S26BA` | `S26TA` | `S26NA` | Complete (WGS + RNA) |

> **Note on Multi-Lane Sequencing (Batch 2 Tumour WGS):**  
> Samples `S19TB`, `S21TB`, `S24TB`, and `S25TB` contain supplementary sequencing lanes merged via `script_1_merge_lane.sh` into `*.combined.fq.gz` prior to BWA alignment.

---

## 3. Execution Pipeline

### RNA-Seq Workflow
1. **Quality Control:** Multi-threaded FastQC across paired FASTQ files.
   ```bash
   bash /data/OSCC2/fastq/RNA/script_1_fastqc.sh
   ```
2. **Alignment & Coordinate Sorting:** Spliced alignment with HISAT2 using the GRCh38 SNP-aware reference (`genome_snp_tran`), piped to `samtools sort`.
   ```bash
   bash /data/OSCC2/fastq/RNA/script_2_hisat.sh
   ```
3. **Assembly & Quantification:** Reference-guided transcriptome reconstruction using StringTie, followed by raw count extraction via `prepDE.py`.
   ```bash
   bash /data/OSCC2/fastq/RNA/script_3_stringtie_assemble.sh
   bash /data/OSCC2/fastq/RNA/script_4_merge.sh
   bash /data/OSCC2/fastq/RNA/script_5_stringtie_quant.sh
   bash /data/OSCC2/fastq/RNA/script_6_extraction.sh
   ```

### WGS Workflow
1. **Multi-Lane Merging & Quality Control:** Dual-lane merging for top-up samples and FastQC validation.
   ```bash
   bash /data/OSCC2/fastq/WGS/scripts/WGS/script_1_merge_lane.sh
   bash /data/OSCC2/fastq/WGS/scripts/WGS/script_1b_qc.sh
   ```
2. **Read Alignment & Deduplication:** BWA-MEM paired-end mapping against GRCh38, followed by Picard `MarkDuplicates` to generate analysis-ready BAM files under `normal/bam_markdup/` and `tumour/bam_markdup/`.
   ```bash
   bash /data/OSCC2/fastq/WGS/scripts/WGS/script_2_bwa.sh
   bash /data/OSCC2/fastq/WGS/scripts/WGS/script_3_markdup.sh
   ```
3. **Somatic Variant Discovery (Mutect2 Engine):** Tumor-Normal paired somatic mutation calling via GATK4 Mutect2 using `1000g_pon.hg38.vcf.gz` (PoN) and `af-only-gnomad.hg38.vcf.gz`. Incorporates `LearnReadOrientationModel` (F1R2 artifact filtering) and `CalculateContamination`.
   ```bash
   bash /data/OSCC2/fastq/WGS/scripts/WGS/script_5_mutect2.sh
   ```
4. **PASS-Only Variant Extraction & Cohort Profiling:** Automated extraction of high-confidence somatic SNVs and INDELs (`-f PASS`), tabix indexing, and generation of cohort mutation burden statistics.
   ```bash
   bash /data/OSCC2/fastq/WGS/scripts/WGS/script_6_extract_pass.sh
   ```

---

## 4. Cohort Somatic Variant Statistics (PASS-Only)

All 23 tumor-normal pairs have completed somatic calling and PASS filtering. Final datasets are archived under `/data/OSCC2/fastq/WGS/mutect2/pass_only_vcfs/`.

| Sample ID | PASS Total | PASS SNV | PASS INDEL | Cohort Batch |
| :--- | :--- | :--- | :--- | :--- |
| **DNA1_M** | 304,003 | 269,980 | 15,830 | Batch 1 (DNA series) |
| **DNA2_M** | 345,086 | 300,117 | 23,177 | Batch 1 (DNA series) |
| **DNA3_M** | 260,904 | 221,389 | 23,528 | Batch 1 (DNA series) |
| **DNA5_M** | 297,714 | 252,049 | 29,143 | Batch 1 (DNA series) |
| **DNA7_M** | 308,677 | 274,905 | 13,420 | Batch 1 (DNA series) |
| **DNA8_M** | 263,590 | 233,137 | 16,801 | Batch 1 (DNA series) |
| **DNA9_M** | 274,451 | 245,845 | 10,986 | Batch 1 (DNA series) |
| **DNA10_M** | 271,941 | 235,375 | 16,788 | Batch 1 (DNA series) |
| **DNA11_M** | 207,469 | 188,675 | 5,875 | Batch 1 (DNA series) |
| **DNA12_M** | 217,927 | 195,833 | 9,142 | Batch 1 (DNA series) |
| **S13TB** | 11,043 | 10,070 | 858 | Batch 2 (S series) |
| **S14TB** | 24,949 | 22,398 | 2,301 | Batch 2 (S series) |
| **S15TB** | 23,014 | 21,963 | 905 | Batch 2 (S series) |
| **S16TB** | 19,687 | 17,933 | 1,561 | Batch 2 (S series) |
| **S18TB** | 25,426 | 23,124 | 1,966 | Batch 2 (S series) |
| **S19TB** | 9,807 | 8,807 | 884 | Batch 2 (S series) |
| **S20TB** | 16,303 | 14,740 | 1,382 | Batch 2 (S series) |
| **S21TB** | 15,083 | 12,379 | 2,513 | Batch 2 (S series) |
| **S22TB** | 20,838 | 19,125 | 1,531 | Batch 2 (S series) |
| **S23TB** | 16,629 | 15,166 | 1,278 | Batch 2 (S series) |
| **S24TB** | 7,545 | 6,846 | 595 | Batch 2 (S series) |
| **S25TB** | 14,029 | 12,897 | 981 | Batch 2 (S series) |
| **S26TB** | 25,979 | 23,841 | 1,842 | Batch 2 (S series) |
