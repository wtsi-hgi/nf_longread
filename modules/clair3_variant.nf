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

workflow clair3_variant {
    take:
    ch_sample
    ch_bam

    main:
    ch_fasta = ch_sample.map { group, fastq, fasta -> [group, fasta] }
    index_fai(ch_fasta)
    ch_fai = index_fai.out.ref_fai

    ch_joined = ch_bam.join(ch_fai)

    clair3_call(ch_joined)
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

process index_fai {
    label 'process_single'

    input:
    tuple val(group), path(fasta)

    output:
    tuple val(group), path("${fasta}"), path("${fasta}.fai"), emit: ref_fai

    script:
    """
    samtools faidx ${fasta}
    """
}

process clair3_call {
    label 'process_high'

    publishDir "${params.outdir}/vcfFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(group), path("clair3_outs/${group}.vcf.gz"), emit: vcf

    script:
    def platform = ''
    def model_path = ''

    switch (params.platform) {
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

    switch (params.model) {
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
                  -f ${fasta} \
                  -m ${model_path} \
                  -t $task.cpus \
                  -p ${platform} \
                  -o clair3_outs \
                  --no_phasing_for_fa \
                  --include_all_ctgs

    if ls clair3_outs/merge_output.vcf.gz &> /dev/null
    then
        mv clair3_outs/merge_output.vcf.gz clair3_outs/${group}.vcf.gz
    else
        echo "No variant found!"
    fi
    """
}
