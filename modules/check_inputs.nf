/*
#########################################
####        checking pipeline        ####
#########################################
*/

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow check_inputs
{
    take:
    sample_sheet

    main:
    Channel
        .fromPath(sample_sheet, checkIfExists: true)
        .splitCsv(header:true, sep:",")
        .map { row -> tuple(row.group, row.input_file, row.fasta) }
        .set { ch_sample_sheet }

    cat_reads(ch_sample_sheet)
    ch_cat_out = cat_reads.out.cat_out

    emit:
    ch_cat_out
}

/*
#~~~~~~~~~#
# process #
#~~~~~~~~~#
*/

process cat_reads
{
    label 'process_single'

    publishDir "${params.outdir}/mergedReads", mode: "copy", overwrite: true

    input:
    tuple val(group), val(input_file), val(fasta)

    output:
    tuple val(group), path("${group}.fastq.gz"), path("${group}.ref.fasta"), emit: cat_out

    script:
    """
    if ls ${input_file}/*.fastq &> /dev/null
    then
        pigz -9 ${input_file}/*.fastq
    fi

    cat ${input_file}/*.fastq.gz > ${group}.fastq.gz
    cp ${fasta} ${group}.ref.fasta
    """
}
