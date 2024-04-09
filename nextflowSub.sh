#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=7.5G

module load nextflow

READS='/data/WHRI-GenomeCentre/shares/Projects_RandD/Reference_Datasets/PolyA/GC-JK-9047_copy/Data/201023_NS500784_0721_AHFCJHBGXG/Alignment_1/20201024_095434/Fastq/fastq_lane-merged/*_L001_R{1,2}_001.fastq.gz'
ANALYSISDIR='/data/WHRI-GenomeCentre/shares/Projects_RandD/Reference_Datasets/PolyA/GC-JK-9047_copy/Analysis/nextflow'

nextflow -DXmx=1G \
         -C nextflow.config \
         run nextflowScriptDev.nf -resume --reads ${READS} --analysisdir ${ANALYSISDIR}
