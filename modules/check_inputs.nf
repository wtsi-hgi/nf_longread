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

workflow check_inputs {
    take:
    sample_sheet

    main:
    Channel
        .fromPath(sample_sheet, checkIfExists: true)
        .splitCsv(header:true, sep:",")
        .map { row -> tuple(row.group, row.fastq, row.fasta) }
        .set { ch_sample }

    cat_reads(ch_sample)
    ch_cat_out = cat_reads.out.cat_out

    Channel
        .fromPath(sample_sheet, checkIfExists: true)
        .splitCsv(header:true, sep:",")
        .map { row -> tuple(row.group.toString(), row.barcode_start.toInteger(), row.barcode_end.toInteger(), row.barcode_template.toString()) }
        .set { ch_barcode }

    emit:
    ch_cat_out
    ch_barcode
}

/*
#~~~~~~~~~#
# process #
#~~~~~~~~~#
*/

process cat_reads {
    label 'process_single'

    publishDir "${params.outdir}/mergedReads", mode: "copy", overwrite: true

    input:
    tuple val(group), val(fastq), val(fasta)

    output:
    tuple val(group), path("${group}.fastq.gz"), path("${group}.ref.fasta"), emit: cat_out

    script:
    """
    if ls ${fastq}/*.fastq &> /dev/null
    then
        pigz -9 ${fastq}/*.fastq
    fi

    cat ${fastq}/*.fastq.gz > ${group}.fastq.gz
    cp ${fasta} ${group}.ref.fasta
    """
}
