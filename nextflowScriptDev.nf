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
    input:
    tuple val(sample_id), path(reads)
    path parentFolder

    output:
    tuple val(sample_id), path("${parentFolder}/1M_Subsample/1M-${sample_id}_R1.fastq"), path("${parentFolder}/1M_Subsample/1M-${sample_id}_R2.fastq"), emit: reads

    script:
    """
    if [ ! -e ${parentFolder}/1M_Subsample ]; then mkdir -p ${parentFolder}/1M_Subsample; fi
    gzip -cd ${reads[0]} | head -4000000 > ${parentFolder}/1M_Subsample/1M-${sample_id}_R1.fastq
    gzip -cd ${reads[1]} | head -4000000 > ${parentFolder}/1M_Subsample/1M-${sample_id}_R2.fastq
    """
}

/*
 * Perform FastQC
 */
process FASTQC {
    input:
    tuple val(sample_id), path(read1), path(read2)
    path parentFolder

    output:
    file "${parentFolder}/1M_fastqc/*_pre_trim.{zip,html}"
    
    script:
    """
    # Load module
    module load fastqc
    
    # Create output folders
    if [ ! -e ${parentFolder}/1M_fastqc ]; then mkdir -p ${parentFolder}/1M_fastqc; fi

    # FastQC files
    fastqc -o ${parentFolder}/1M_fastqc ${read1} ${read2}
    """
}

/*
 * Perform trimming
 */
process TRIMMING {
    input:
    tuple val(sample_id), path(reads)
    path parentFolder

    output:
    tuple val(sample_id), path("${parentFolder}/trimgalore_outs/*val*.fq.gz"), emit: reads
	path("${parentFolder}/trimgalore_outs/*report.txt"), optional: true, emit: report

    script:
    """
    # Load module
    module load trimgalore/0.6.5

    # Create output folders
    if [ ! -e ${parentFolder}/trimgalore_outs ]; then mkdir -p ${parentFolder}/trimgalore_outs; fi
    
	# Run trimgalore
	# Options --length, -e, --stringency, --quality are set to default
	trim_galore --quality 20 \
	--length 20 \
	-j ${NSLOTS} \
	--paired \
	--stringency 1 \
	-e 0.1 \
	--output_dir ${parentFolder}/trimgalore_outs \
	${reads}
    """
}

/*
 * Perform FastQC, post trimming
 */
process FASTQC_PT {
    input:
    tuple val(sample_id), path(reads)
    path parentFolder

    output:
    file "${parentFolder}/post_trim_fastqc/*_post_trim.{zip,html}"

    script:
    """
    # Load module
    module load fastqc

    # Create output folders
    if [ ! -e ${parentFolder}/post_trim_fastqc ]; then mkdir -p ${parentFolder}/post_trim_fastqc; fi
    
    # FastQC files
    fastqc -o ${parentFolder}/post_trim_fastqc ${reads}
    """
}

/*
 * Perform MultiQC
 */
process MULTIQC {
    conda 'multiqc=1.21'

    input:
    file '*'

    output:
    path "$params.analysisdir/multiqc_report.html"

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
    SUBSAMPLE(read_pairs_ch, params.analysisdir)
    fastqc_ch = FASTQC(SUBSAMPLE.out.reads, params.analysisdir)
    TRIMMING(read_pairs_ch, params.analysisdir)
    fastqcPT_ch = FASTQC_PT(TRIMMING.out.reads, params.analysisdir)
    MULTIQC(fastqcPT_ch.mix(fastqc_ch).collect())
}