# rnaseq_nextflow

## ❓ What
- I am using this repository as a training environment to learn NextFlow

## 🎯 Aim
- Learn how to use NextFlow

## Details
- This pipeline is designed to run on the QMUL HPC (Apocrita)
- Pipeline:
    1. Subset FASTQ files down to 1M reads
    2. Run FASTQC on subsetted FASTQ files
    3. Run TRIMGALORE on full FASTQ files
    4. Run FASTQC on trimming output
    5. Run MULTIQC on outputs from FASTQC on subsetted reads and trimmed reads
