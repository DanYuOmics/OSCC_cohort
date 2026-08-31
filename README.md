# OSCC Multi-Omics Sequencing Pipeline (WGS & RNA-Seq)

This repository contains the end-to-end data processing pipelines, quality control protocols, and alignment scripts for the Oral Squamous Cell Carcinoma (OSCC) multi-omics cohort. The dataset integrates whole-genome sequencing (WGS) and bulk RNA-sequencing (RNA-Seq) across two independent sequencing batches.

---

## 1. Directory Structure and File Paths

All raw FASTQ files, reference databases, and downstream output files are organized systematically under the central storage environment:

```text
/data/OSCC2/
├── fastq/
│   ├── RNA/                               # RNA-Seq raw & symlinked FASTQs (28 libraries)
│   │   ├── samples.txt                    # RNA-Seq sample prefix list
│   │   ├── script_1_fastqc.sh             # FastQC execution script
│   │   ├── script_2_hisat.sh              # Parallel HISAT2 alignment script
│   │   └── hisat2_out/                    # Coordinate-sorted BAM outputs (*.bam, *.bam.bai)
│   └── WGS/                               # WGS raw & symlinked FASTQs
│       ├── normal/                        # Germline / normal control FASTQs (25 samples)
│       └── tumour/                        # Primary tumor FASTQs (23 samples / 27 sequencing runs)
├── fastqc_out/                            # FastQC HTML and zip quality reports
└── copy_wgs.log                           # WGS disk-to-server data transfer log

/data/yudan_ref/
└── RNA/
    └── grch38_snp_tran/                   # HISAT2 pre-built index with SNP and transcript models
        ├── genome_snp_tran.1.ht2
        ├── genome_snp_tran.2.ht2
        └── ...
```

---

## 2. Sample Matrix & Multi-Omics Availability

The cohort comprises **25 individual patients** profiled across two sequencing batches:
* **Batch 1 (Patients 1–12):** 10 complete pairs with both WGS and RNA-Seq, alongside 2 Normal-only controls (Patients 4 and 6).
* **Batch 2 (Patients 13–26, excluding 17):** 13 complete pairs with full WGS and RNA-Seq coverage.

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
> Samples `S19TB`, `S21TB`, `S24TB`, and `S25TB` contain supplementary sequencing lanes (e.g. `V350396689_L03` / `V350396689_L04`) to satisfy genomic coverage requirements. These are merged during alignment and coordinate sorting.

---

## 3. Execution Pipeline

### RNA-Seq Workflow
1. **Quality Control:** Run multi-threaded FastQC across all paired FASTQ files.
   ```bash
   bash /data/OSCC2/fastq/RNA/script_1_fastqc.sh
   ```
2. **Parallel Alignment & Indexing:** Spliced alignment with HISAT2 using the GRCh38 SNP-aware reference (`genome_snp_tran`), piped into `samtools sort` and `samtools index`.
   ```bash
   bash /data/OSCC2/fastq/RNA/script_2_hisat.sh
   ```

### WGS Workflow
1. **Data Consolidation:** Standardized symlinking of Batch 1 and direct staging of Batch 2 into unified `normal/` (25 samples) and `tumour/` (27 runs) directories.
2. **Read Alignment & BAM Processing:** BWA-MEM2 paired mapping against GRCh38, followed by Picard MarkDuplicates and GATK4 BaseRecalibrator.
3. **Somatic Variant Discovery:** Paired Tumor-Normal variant calling with GATK4 Mutect2.
