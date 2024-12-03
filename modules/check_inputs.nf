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
        .map { row -> tuple(row.group, row.directory, row.reference, row.gene_info) }
        .set { ch_input }

    cat_reads(ch_input)
    ch_sample = cat_reads.out.ch_sample

    ch_sample.flatMap { group, directory, reference, gene_info ->
                if (!gene_info.exists()) {
                    throw new RuntimeException("The gene_info file does not exist for group: ${group} (${gene_info})!")
                }

                gene_info
                    .readLines()
                    .collect { line -> 
                        def fields = line.split('\t')
                        tuple(group, fields[0], fields[1] as Integer, fields[2] as Integer) } }
             .set { ch_gene }
    
    Channel
        .fromPath(sample_sheet, checkIfExists: true)
        .splitCsv(header:true, sep:",")
        .map { row -> 
            def group = row.group.toString()
            def barcode_start = row.barcode_start ? row.barcode_start.toInteger() : null
            def barcode_end = row.barcode_end ? row.barcode_end.toInteger() : null
            def barcode_template = row.barcode_template?.toString()
            tuple(group, barcode_start, barcode_end, barcode_template) }
        .set { ch_barcode }

    emit:
    ch_sample
    ch_gene
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
    tuple val(group), val(directory), val(reference), val(gene_info)

    output:
    tuple val(group), path("${group}.fastq.gz"), path("${group}.ref.fasta"), path("${group}.gene_info.txt"), emit: ch_sample

    script:
    """
    if ls ${directory}/*.fastq &> /dev/null
    then
        pigz -9 ${directory}/*.fastq
    fi

    cat ${directory}/*.fastq.gz > ${group}.fastq.gz
    cp ${directory}/${reference} ${group}.ref.fasta
    cp ${directory}/${gene_info} ${group}.gene_info.txt
    """
}
