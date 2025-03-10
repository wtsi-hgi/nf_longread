/*
#########################################
####        longread pipeline        ####
#########################################
*/

/*
#~~~~~~~~~~~~~~~~~~#
# define functions #
#~~~~~~~~~~~~~~~~~~#
*/

def helpMessage() {
    log.info """
    Usage:
    nextflow run nf_longread/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet        Path of the sample sheet
    
    Optional arguments:
    Basic:
        --outdir              the directory path of output results, default: the current directory
    
    Alignment:
        --protocol            DNA, cDNA, directRNA, default: DNA
        --platform            nanopore, pacbio, hifi, default: nanopore
    
    Barcode Detection:
        --mapq                the mapping quality for filtering, default: 1
        --qualcut             the base quality in the barcode for filtering , default: 10
        --numcut              the number of low-quality bases in the barcode for filtering, default: 3
        --countcut            the number of reads supporting the barcode for filtering, default: 5

    Extract SNVs:
        --basequal            the base quality for filtering, default: 30
        --region              the expected region of variants, eg: 100,200, default: 0,0

    Step arguments:
        --skip_align          skip alignment
        --skip_barcode        skip barcode detection
        --skip_snvcov         skip snv coverage extraction
    """
}

/*
#~~~~~~~~~~~~~~~~~~#
# check parameters #
#~~~~~~~~~~~~~~~~~~#
*/
params.help = null
if(params.help) {
    helpMessage()
    exit 0
}

if (params.sample_sheet) {
    ch_input = file(params.sample_sheet)
} else {
    helpMessage()
    println "========> Error!!! Please specify the full path of the sample sheet!\n"
    exit 1
}

/* output options */
params.outdir       = params.outdir       ?: "$PWD"

/* platform options */
params.protocol     = params.protocol     ?: "DNA"
params.platform     = params.platform     ?: "nanopore"

/* barcode detection options */
params.libtype      = params.libtype      ?: "muta"
params.mapq         = params.mapq         ?: 1
params.qualcut      = params.qualcut      ?: 10
params.numcut       = params.numcut       ?: 3
params.countcut     = params.countcut     ?: 5

/* snv plot options */
params.basequal     = params.basequal     ?: 30
params.region       = params.region       ?: "0,0"

/* step options */
params.skip_align   = params.skip_align   ?: false
params.skip_barcode = params.skip_barcode ?: false
params.skip_snvcov  = params.skip_snvcov  ?: false

if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
}

if (params.platform != 'nanopore' && params.platform != 'pacbio' && params.platform != 'hifi') {
    exit 1, "Invalid protocol option: ${params.platform}. Valid options: 'nanopore', 'pacbio', 'hifi'"
}

if (params.libtype != 'muta' && params.libtype != 'gene') {
    exit 1, "Invalid protocol option: ${params.libtype}. Valid options: 'muta', 'gene'"
}

if (!params.skip_barcode) {
    if(params.skip_align) {
        exit 1, "Cannot run barcode detection with --skip_align!"
    }
}

if (!params.skip_snvcov) {
    if(params.skip_align) {
        exit 1, "Cannot run snv coverage extraction with --skip_align!"
    }
}
/*
#~~~~~~~~~~~~~~#
# load modules #
#~~~~~~~~~~~~~~#
*/
include { check_inputs } from '../modules/check_inputs.nf'
include { minimap2_align } from '../modules/minimap2_align.nf'
include { detect_barcode } from '../modules/detect_barcode.nf'
include { extract_snvs } from '../modules/extract_snvs.nf'

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow longread {
    check_inputs(ch_input)
    ch_sample = check_inputs.out.ch_sample
    ch_gene = check_inputs.out.ch_gene
    ch_barcode = check_inputs.out.ch_barcode

    ch_sample_input = ch_sample.map { group, fastq, fasta, gene_info -> tuple(group, fastq, fasta) }

    if (!params.skip_align) {
        minimap2_align(ch_sample_input)
        ch_bam = minimap2_align.out.ch_bam   
    }

    if (!params.skip_barcode) {
        ch_validated_barcode = ch_barcode.map { group, barcode_start, barcode_end, barcode_template ->
            if (!barcode_start || !barcode_end) { error "Error: barcode_start or barcode_end is empty. Barcode detection cannot proceed." }
            return tuple(group, barcode_start, barcode_end, barcode_template)
        }

        detect_barcode(ch_sample_input, ch_gene, ch_validated_barcode, ch_bam)
    }

    if (!params.skip_snvcov) {
        extract_snvs(ch_sample_input, ch_bam)
    }
}
