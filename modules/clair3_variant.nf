/*
#######################################
####        clair3 pipeline        ####
#######################################
*/

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow clair3_variant
{
    take:
    ch_sample
    ch_bam

    main:
    ch_fai = index_fai(ch_sample)

    clair3_call(ch_sample, ch_bam, ch_fai)
    ch_vcf = clair3_call.out.vcf

    emit:
    ch_sample
    ch_vcf
}

/*
#~~~~~~~~~#
# process #
#~~~~~~~~~#
*/

process index_fai
{
    label 'process_single'

    input:
    tuple val(group), path(fastq), path(fasta)

    output:
    path("${group}.ref.fasta"), emit: ref_fai

    script:
    """
    samtools faidx ${fasta}
    """
}

process clair3_call
{
    label 'process_high'

    publishDir "${params.outdir}/vcfFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(fastq), path(fasta)
    path bam
    path reference

    output:
    tuple path("clair3/${group}.vcf.gz"), path("clair3/${group}.vcf.gz.tbi"), emit: vcf

    script:
    def platform = ''
    def model_path = ''

    switch (params.platform)
    {
        case 'nanopore':
            platform = "ont"
            break
        case 'pcabio':
            platform = "hifi"
            break
        case 'hifi':
            platform = "hifi"
            break
    }

    switch (params.model)
    {
        case 'ont_r10':
            model_path = "${projectDir}/resources/r1041_e82_400bps_sup_v420"
            break
        case 'ont_r9':
            model_path = "${projectDir}/resources/r941_prom_sup_g5014"
            break
        case 'hifi':
            model_path = "${projectDir}/resources/hifi"
            break
    }

    """
    run_clair3.sh -b ${bam} \
                  -f ${reference} \
                  -m ${model_path} \
                  -t $task.cpus \
                  -p ${platform} \
                  -o clair3_outs \
                  --include_all_ctgs

    mv clair3_outs/merge_output.vcf.gz clair3_outs/${group}.vcf.gz
    mv clair3_outs/merge_output.vcf.gz.tbi clair3_outs/${group}.vcf.gz.tbi
    """
}
