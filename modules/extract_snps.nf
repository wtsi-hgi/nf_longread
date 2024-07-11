/*
###############################################
####        extracting snp pipeline        ####
###############################################
*/

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow extract_snps {
    take:
    ch_sample
    ch_bam

    main:
    samtools_filtering(ch_bam)
    ch_snps_bam = samtools_filtering.out.snps_bam
    ch_sample_bam = ch_sample.join(ch_snps_bam)

    get_snps(ch_sample_bam)
    ch_snp_cov = get_snps.out.ch_snp_cov

    emit:
    ch_snp_cov
}

/*
#~~~~~~~~~#
# process #
#~~~~~~~~~#
*/

process samtools_filtering {
    label 'process_medium'

    publishDir "${params.outdir}/bamFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(bam), path(bai)

    output:
    tuple val(group), path("${group}.snps.bam"), path("${group}.snps.bam.bai"), emit: snps_bam

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -F 2308 -q ${params.mapq} -o ${group}.snps.bam $bam
    samtools index -@ $task.cpus ${group}.snps.bam
    """
}


process get_snps {
    label 'process_single'

    publishDir "${params.outdir}/snpCov/${group}", mode: "copy", overwrite: true

    input:
    tuple val(group), path(fastq), path(fasta), 
          path(bam), path(bai)

    output:
    path "${group}.snp_cov.txt", emit: ch_snp_cov
    path "${group}.snp_cov.png", emit: ch_snp_png

    script:  
    """
    python ${projectDir}/scripts/extract_snps.py -i ${bam} -o .
    """ 
}