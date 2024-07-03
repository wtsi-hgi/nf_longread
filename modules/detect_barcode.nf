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
    ch_list = extract_barcode.out.ch_list

    emit:
    ch_list
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
    path "${target}.barcodes.txt", emit: ch_list

    script:
    def barcode_opt = ''

    if (sequence != '') {
        barcode_opt = "-b ${sequence}"
    }

    if (start == null || start.toString().trim() == "" || end == null || end.toString().trim() == "") {
        throw new IllegalArgumentException("Error: start or end is empty for group: ${group}")
    }
    
    """
    awk '{if(\$1==">${target}"){print \$0; getline; print \$0; exit}}' ${fasta} > ${target}.fa
    samtools faidx ${target}.fa

    echo -e "${target}\\t1\\t\$((${start}-1))\\tleft_barcode\\t60\\t+" > ${target}.left_barcode.bed
    echo -e "${target}\\t${end}\\t\$((${end}+30))\\tright_barcode\\t60\\t-" > ${target}.right_barcode.bed

    samtools view -h -O SAM -o ${target}.sam ${bam} ${target}

    awk -F'\\t' -v OFS='\\t' '{if(!(\$0~/^@/)){\$2=0}; print \$0}' ${target}.sam | samtools view -bS - > ${target}.plus.bam
    samtools ampliconclip -b ${target}.left_barcode.bed -o ${target}.left_tmp.bam -f ${target}.left.txt --hard-clip --strand --clipped -O bam --tolerance 1000000 ${target}.plus.bam
    samtools view -b -h -O BAM -F 4 -o ${target}.left.bam ${target}.left_tmp.bam

    samtools view -h ${target}.left.bam | awk -F'\\t' -v OFS='\\t' '{if(!(\$0~/^@/)){\$2=16}; print \$0}' | samtools view -bS - > ${target}.minus.bam
    samtools ampliconclip -b ${target}.right_barcode.bed -o ${target}.right_tmp.bam -f ${target}.right.txt --hard-clip --strand --clipped -O bam --tolerance 1000000 ${target}.minus.bam
    samtools view -b -h -O BAM -F 4 -o ${target}.right.bam ${target}.right_tmp.bam

    samtools view -h -O SAM -o ${target}.right.sam ${target}.right.bam
    rm ${target}.right.bam

    awk -F'\\t' -v OFS='\\t' 'NR==FNR{if(!(\$0~/^@/)){a[\$1]=\$2};next}{if(!(\$0~/^@/)){\$2=a[\$1]};print \$0}' ${target}.sam ${target}.right.sam > ${target}.right.tmp.sam
    rm ${target}.right.sam
    mv ${target}.right.tmp.sam ${target}.right.sam

    samtools view -bS ${target}.right.sam > ${target}.right.bam
    samtools index ${target}.right.bam
    
    rm *.sam

    samtools fixmate -O bam ${target}.right.bam ${target}.right.fixmate.bam
    samtools calmd ${target}.right.fixmate.bam ${target}.fa --output-fmt bam > ${target}.barcode.bam
    samtools index ${target}.barcode.bam

    python ${projectDir}/scripts/extract_barcodes.py -i ${target}.barcode.bam -s ${start} -e ${end} -o . ${barcode_opt}

    rm *.bam
    rm *.fa
    """ 
}