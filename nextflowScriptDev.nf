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
    publishDir "${params.analysisdir}/1M_Subsample"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("1M-${sample_id}_R1.fastq"), path("1M-${sample_id}_R2.fastq"), emit: reads

    script:
    """
    gzip -cd ${reads[0]} | head -4000000 > 1M-${sample_id}_R1.fastq
    gzip -cd ${reads[1]} | head -4000000 > 1M-${sample_id}_R2.fastq
    """
}

/*
 * Perform FastQC
 */
process FASTQC {
    publishDir "${params.analysisdir}/1M_fastqc"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    file "*fastqc*" into ch_fastqc
    
    script:
    """
    # Load module
    module load fastqc

    # FastQC files
    fastqc ${read1} ${read2}
    """
}

/*
 * Perform trimming
 */
process TRIMMING {
    publishDir "${params.analysisdir}/trimgalore_outs"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*val*.fq.gz"), emit: reads
	path("*report.txt"), optional: true, emit: report

    script:
    """
    # Load module
    module load trimgalore/0.6.5
    
	# Run trimgalore
	# Options --length, -e, --stringency, --quality are set to default
	trim_galore --quality 20 \
	--length 20 \
	-j ${NSLOTS} \
	--paired \
	--stringency 1 \
	-e 0.1 \
	${reads}
    """
}

/*
 * Perform FastQC, post trimming
 */
process FASTQC_PT {
    publishDir "${params.analysisdir}/post_trim_fastqc"

    input:
    tuple val(sample_id), path(reads)

    output:
    file "*fastqc*" into ch_fastqc_trim

    script:
    """
    # Load module
    module load fastqc

    # FastQC files
    fastqc ${reads}
    """
}

/*
 * Perform MultiQC
 */
process MULTIQC {
    conda 'multiqc=1.21'

    input:
    file("*")
    file("*")

    output:
    path("${params.analysisdir}/multiqc_report.html")

    script:
    """
    # Run multiqc
    multiqc .
    """
}

/*
 * Define workflow
 */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    SUBSAMPLE(read_pairs_ch)
    FASTQC(SUBSAMPLE.out.reads)
    TRIMMING(read_pairs_ch)
    FASTQC_PT(TRIMMING.out.reads)
    MULTIQC(
        ch_fastqc.collect(), 
        ch_fastqc_trim.collect())
}