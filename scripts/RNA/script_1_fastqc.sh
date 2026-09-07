#!/bin/bash
FASTQ_PATH=/data/OSCC2/fastq/RNA
OUT_DIR=/data/OSCC2/fastqc_out

mkdir -p ${OUT_DIR}
fastqc -o ${OUT_DIR} --threads 64 ${FASTQ_PATH}/*.fq.gz
