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
        .map { row -> 
            def group = row.group.toString()
            def barcode_start = row.barcode_start ? row.barcode_start.toInteger() : null
            def barcode_end = row.barcode_end ? row.barcode_end.toInteger() : null
            def barcode_template = row.barcode_template?.toString()
            tuple(group, barcode_start, barcode_end, barcode_template) 
        }
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
