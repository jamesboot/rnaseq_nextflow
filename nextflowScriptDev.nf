/*
 * Pipeline input parameters
 */

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    analysisdir  : ${params.analysisdir}
    """
    .stripIndent(true)

/*
 * Sub-sample to 250k reads, 1M lines
 */
process SUBSAMPLE {
    tag "Sub-sample on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)
    path parentFolder

    output:
    tuple val(sample_id), path("${parentFolder}/1M-Subsample/1M-${sample_id}_R1.fastq")
    tuple val(sample_id), path("${parentFolder}/1M-Subsample/1M-${sample_id}_R2.fastq")

    script:
    """
    if [ ! -e ${parentFolder}/1M-Subsample ]; then mkdir -p ${parentFolder}/1M-Subsample; fi
    gzip -cd ${reads[0]} | head -4000000 > ${parentFolder}/1M-Subsample/1M-${sample_id}_R1.fastq
    gzip -cd ${reads[1]} | head -4000000 > ${parentFolder}/1M-Subsample/1M-${sample_id}_R2.fastq
    """
}

/*
 * Sub-sample to 250k reads, 1M lines
 */
process FASTQC {
    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(read1)
    tuple val(sample_id), path(read2)
    path parentFolder

    output:
    path "${parentFolder}/fastqc"

    script:
    """
    # Load module
    module load fastqc

    # Make output folder
    if [ ! -e ${parentFolder}/1M-fastqc ]; then mkdir -p ${parentFolder}/1M-fastqc; fi
    
    # FastQC files
    fastqc --outdir=${parentFolder}/1M-fastqc ${read1}
    fastqc --outdir=${parentFolder}/1M-fastqc ${read2}
    """
}

/*
 * Define workflow
 */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    SUBSAMPLE(read_pairs_ch, params.analysisdir)
    FASTQC(SUBSAMPLE.out[0], SUBSAMPLE.out[1], params.analysisdir)
}