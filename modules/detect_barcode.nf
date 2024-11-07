/*
########################################
####        barcode pipeline        ####
########################################
*/

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow detect_barcode {
    take:
    ch_sample
    ch_barcode
    ch_bam

    main:
    samtools_filtering(ch_bam)
    ch_filtered_bam = samtools_filtering.out.filtered_bam

    ch_sample_bam = ch_sample.join(ch_filtered_bam)
    ch_sample_bam_barcode = ch_sample_bam.join(ch_barcode)

    samtools_stat(ch_filtered_bam)
    ch_idxstats = samtools_stat.out.idxstats

    ch_idxstats_list = ch_idxstats.map { group, idxstatsFile -> 
        idxstatsFile.text.split('\n')
            .findAll { it.trim() }
            .collect { line ->
                def cols = line.split('\t')
                if (!cols[0].contains("*")) {
                    tuple(group, cols[0], cols[1], cols[2])
                }
            }
            .findAll { it != null }
    }

    ch_idxstats_flatten = ch_idxstats_list.flatMap { it }

    ch_target = ch_sample_bam_barcode.combine(ch_idxstats_flatten, by: 0)
    //ch_target.view()
    
    extract_barcode(ch_target)
    ch_target_barcode = extract_barcode.out.ch_target_barcode

    emit:
    ch_target_barcode
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
    tuple val(group), path("${group}.filtered.bam"), path("${group}.filtered.bam.bai"), emit: filtered_bam

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -F 0x800 -q ${params.mapq} -o ${group}.filtered.bam $bam
    samtools index -@ $task.cpus ${group}.filtered.bam
    """
}

process samtools_stat {
    label 'process_single'

    publishDir "${params.outdir}/bamFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(bam), path(bai)

    output:
    tuple val(group), path("${group}.filtered.flagstat.txt"), emit: flagstat
    tuple val(group), path("${group}.filtered.idxstats.txt"), emit: idxstats
    
    script:
    """
    samtools flagstat ${bam} > ${group}.filtered.flagstat.txt
    samtools idxstats ${bam} > ${group}.filtered.idxstats.txt
    """
}

process extract_barcode {
    label 'process_single'

    publishDir "${params.outdir}/barcodeList/${group}", mode: "copy", overwrite: true

    input:
    tuple val(group), path(fastq), path(fasta), 
          path(bam), path(bai), 
          val(start), val(end), val(sequence),
          val(target), val(length), val(reads)

    output:
    path "${target}.barcodes.txt", emit: ch_target_barcode
    tuple path("${target}.bam"), path("${target}.bam.bai"), emit: ch_target_bam

    script:
    def barcode_opt = ''

    if (sequence != '') {
        barcode_opt = "-b ${sequence}"
    }

    if (start == null || start.toString().trim() == "" || end == null || end.toString().trim() == "") {
        throw new IllegalArgumentException("Error: start or end is empty for group: ${group}")
    }
    
    """
    samtools view -h -O BAM -o ${target}.bam ${bam} ${target}
    samtools index ${target}.bam

    python ${projectDir}/scripts/extract_barcodes.py -i ${target}.bam -s ${start} -e ${end} -o . ${barcode_opt} -q ${params.qualcut} -n ${params.numcut} -c ${params.countcut}
    """ 
}