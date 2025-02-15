/*
#########################################
####        minimap2 pipeline        ####
#########################################
*/

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow minimap2_align {
    take:
    ch_sample

    main:
    align_reads(ch_sample)
    ch_sam = align_reads.out.sam

    samtools_view_bam(ch_sam)
    ch_bam = samtools_view_bam.out.bam

    samtools_flagstat(ch_bam)

    emit:
    ch_bam
}

/*
#~~~~~~~~~#
# process #
#~~~~~~~~~#
*/

process align_reads {
    label 'process_high'

    input:
    tuple val(group), path(fastq), path(fasta)

    output:
    tuple val(group), path("${group}.sam"), emit: sam

    script:
    def preset = ''
    def kmer = ''
    def stranded = ''
    def tag = ''

    if (params.protocol == 'DNA') {
        switch (params.platform) {
            case 'nanopore':
                preset = "-x map-ont -O 8,28 -B 2"
                break
            case 'pcabio':
                preset = "-x map-pb -O 8,28 -B 2"
                break
            case 'hifi':
                preset = "-x map-hifi -A 1 -B 2 -O 12,30 -E 8,4"
                break
        }
    } else {
        preset = "-x splice"
        kmer = (params.protocol == 'directRNA') ? "-k 14" : "-k 15"
        stranded = (params.protocol == 'directRNA') ? "-uf" : ""
    }

    if (!params.skip_barcode) {
        tag = "--cs"
    } else {
        tag = "--MD"
    }

    """
    minimap2 -a -t $task.cpus $preset $kmer $stranded $tag --secondary=no ${fasta} ${fastq} > ${group}.sam
    """
}

process samtools_view_bam {
    label 'process_medium'

    publishDir "${params.outdir}/bamFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(sam)

    output:
    tuple val(group), path("${group}.sorted.bam"), path("${group}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -o ${group}.bam ${sam}
    samtools sort -@ $task.cpus -o ${group}.sorted.bam ${group}.bam
    samtools index -@ $task.cpus ${group}.sorted.bam
    """
}

process samtools_flagstat {
    label 'process_single'

    publishDir "${params.outdir}/bamFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(bam), path(bai)

    output:
    path "${group}.flagstat.txt"
    
    script:
    """
    samtools flagstat ${bam} > ${group}.flagstat.txt
    """
}