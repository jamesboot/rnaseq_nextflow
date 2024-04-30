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
    tuple val(sample_id), path("${parentFolder}/1M_Subsample/1M-${sample_id}_R1.fastq")
    tuple val(sample_id), path("${parentFolder}/1M_Subsample/1M-${sample_id}_R2.fastq")

    script:
    """
    if [ ! -e ${parentFolder}/1M-Subsample ]; then mkdir -p ${parentFolder}/1M-Subsample; fi
    gzip -cd ${reads[0]} | head -4000000 > ${parentFolder}/1M_Subsample/1M-${sample_id}_R1.fastq
    gzip -cd ${reads[1]} | head -4000000 > ${parentFolder}/1M_Subsample/1M-${sample_id}_R2.fastq
    """
}

/*
 * Perform FastQC
 */
process FASTQC {
    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(read1)
    tuple val(sample_id), path(read2)
    path parentFolder

    output:
    path "${parentFolder}/1M_fastqc"

    script:
    """
    # Load module
    module load fastqc

    # Make output folder
    if [ ! -e ${parentFolder}/1M_fastqc ]; then mkdir -p ${parentFolder}/1M_fastqc; fi
    
    # FastQC files
    fastqc --outdir=${parentFolder}/1M_fastqc ${read1}
    fastqc --outdir=${parentFolder}/1M_fastqc ${read2}
    """
}

/*
 * Perform MultiQC
 */
process MULTIQC {
    input:
    path fastqcFolder

    output:
    file "${fastqcFolder}/pre_trim_multiqc_report.html"

    script:
    """
    # Load gcenv containing multiqc
    source /data/WHRI-GenomeCentre/gcenv/bin/activate

    # Run multiqc
    cd ${fastqcFolder}
    multiqc .
    """
}

/*
 * Perform trimming
 */
process TRIMMING {
    tag "Trimming on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)
    path parentFolder

    output:
    path "${parentFolder}/trimgalore_outs"

    script:
    """
    # Load module
    module load trimgalore/0.6.5

    # Make output folder
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
	${reads[0]} ${reads[1]}
    """
}

/*
 * Perform FastQC, post trimming
 */
process FASTQC_PT {
    tag "Post-trimming FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)
    path trimFolder

    output:
    path "${trimFolder}/post_trim_fastqc"

    script:
    """
    # Load module
    module load fastqc

    # Make output folder
    if [ ! -e ${trimFolder}/post_trim_fastqc ]; then mkdir -p ${trimFolder}/post_trim_fastqc; fi
    
    # FastQC files
    fastqc --outdir=${trimFolder}/post_trim_fastqc ${reads[0]}
    fastqc --outdir=${trimFolder}/post_trim_fastqc ${reads[1]}
    """
}

/*
 * Perform MultiQC, post-trimming
 */
process MULTIQC_PT {
    input:
    path fastqcFolder

    output:
    file "${fastqcFolder}/post_trim_multiqc_report.html"

    script:
    """
    # Load gcenv containing multiqc
    source /data/WHRI-GenomeCentre/gcenv/bin/activate

    # Run multiqc
    cd ${fastqcFolder}
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
    FASTQC(SUBSAMPLE.out[0], SUBSAMPLE.out[1], params.analysisdir)
    MULTIQC(FASTQC.out)
    TRIMMING(read_pairs_ch, params.analysisdir)
    Channel
        .fromFilePairs(TRIMMING.out, checkIfExists: true)
        .set { trim_read_pairs_ch }
    FASTQC(trim_read_pairs_ch, TRIMMING.out)
    MULTIQC_PT(FASTQC.out)
}