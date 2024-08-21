/*
###############################################
####        extracting snv pipeline        ####
###############################################
*/

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow extract_snvs {
    take:
    ch_sample
    ch_bam

    main:
    samtools_filtering(ch_bam)
    ch_snvs_bam = samtools_filtering.out.snvs_bam
    ch_sample_bam = ch_sample.join(ch_snvs_bam)

    get_snvs(ch_sample_bam)
    ch_snv_cov = get_snvs.out.ch_snv_cov

    emit:
    ch_snv_cov
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
    tuple val(group), path("${group}.snvs.bam"), path("${group}.snvs.bam.bai"), emit: snvs_bam

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -F 2308 -q ${params.mapq} -o ${group}.snvs.bam $bam
    samtools index -@ $task.cpus ${group}.snvs.bam
    """
}


process get_snvs {
    label 'process_single'

    publishDir "${params.outdir}/snvCov/${group}", mode: "copy", overwrite: true

    input:
    tuple val(group), path(fastq), path(fasta), 
          path(bam), path(bai)

    output:
    path "${group}.snv_cov.txt", emit: ch_snv_cov
    path "${group}.snv_cov.png", emit: ch_snv_png

    script:  
    """
    python ${projectDir}/scripts/extract_snvs.py -i ${bam} -o . -b ${params.basequal} -r ${params.region}
    """ 
}