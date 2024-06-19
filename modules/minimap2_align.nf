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

workflow minimap2_align
{
    take:
    ch_sample

    main:
    index_mmi(ch_sample)
    ch_ref_mmi = index_mmi.out.ref_mmi

    align_reads(ch_sample, ch_ref_mmi)
    ch_sam = align_reads.out.sam

    samtools_view_bam(ch_sample, ch_sam)
    ch_bam = samtools_view_bam.out.bam

    emit:
    ch_bam
}

/*
#~~~~~~~~~#
# process #
#~~~~~~~~~#
*/

process index_mmi
{
    label 'process_high'

    input:
    tuple val(group), path(fastq), path(fasta)

    output:
    path "${group}.ref.fasta.mmi", emit: ref_mmi

    script:
    def preset = ''
    def kmer = ''
    def stranded = ''

    if (params.protocol == 'DNA') 
    {
        switch (params.platform)
        {
            case 'nanopore':
                preset = "map-ont"
                break
            case 'pcabio':
                preset = "map-pb"
                break
            case 'hifi':
                preset = "map-hifi"
                break
        }
        kmer = "15"
    } else {
        preset = "splice"
        kmer = (params.protocol == 'directRNA') ? "14" : "15"
        stranded = (params.protocol == 'directRNA') ? "-uf" : ""
    }

    """
    minimap2 -x $preset -k $kmer $stranded -t $task.cpus -d ${group}.ref.fasta.mmi $fasta
    """
}

process align_reads
{
    label 'process_high'

    input:
    tuple val(group), path(fastq), path(fasta)
    path index

    output:
    path "${group}.sam", emit: sam

    script:
    def preset = ''
    def kmer = ''
    def stranded = ''
    def md = ''

    if (params.protocol == 'DNA') 
    {
        switch (params.platform)
        {
            case 'nanopore':
                preset = "-x map-ont"
                break
            case 'pcabio':
                preset = "-x map-pb"
                break
            case 'hifi':
                preset = "-x map-hifi"
                break
        }
        kmer = "-k 15"
        md = "--MD"
    } else {
        preset = "-x splice"
        kmer = (params.protocol == 'directRNA') ? "-k 14" : "-k 15"
        stranded = (params.protocol == 'directRNA') ? "-uf" : ""
    }

    """
    minimap2 -a -t $task.cpus $preset $kmer $stranded $md ${index} $fastq > ${group}.sam
    """
}

process samtools_view_bam
{
    label 'process_medium'

    publishDir "${params.outdir}/bamFiles", mode: "copy", overwrite: true

    input:
    tuple val(group), path(fastq), path(fasta)
    path sam

    output:
    tuple path("${group}.sorted.bam"), path("${group}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -o ${group}.bam $sam
    samtools sort -@ $task.cpus -o ${group}.sorted.bam ${group}.bam
    samtools index -@ $task.cpus ${group}.sorted.bam
    """
}

