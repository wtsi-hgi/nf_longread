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
    nextflow run check_inputs.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet        Path of the sample sheet
    
    Optional arguments:
    Basic:
        --outdir              the directory path of output results, default: the current directory
    
    Alignment:
        --protocol            DNA, cDNA, directRNA, default: DNA
        --platform            nanopore, pacbio, hifi, default: nanopore
    
    Variant Calling:
        --model               the trainning model of variant calling, default: ont_r10
    
    Barcode Detection:
        --mapq                the mapping quality for filtering, default: 1

    Extract SNVs:
        --basequal            the base quality for filtering, default: 30
        --region              the expected region of variants, eg: 100,200, default: null,null

    Step arguments:
        --skip_align          skip alignment
        --skip_variant        skip variant calling
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

params.outdir       = params.outdir       ?: "$PWD"
params.protocol     = params.protocol     ?: "DNA"
params.platform     = params.platform     ?: "nanopore"
params.model        = params.model        ?: "ont_r10"
params.mapq         = params.mapq         ?: 1
params.basequal     = params.basequal     ?: 30
params.region       = params.region       ?: "0,0"

params.skip_align   = params.skip_align   ?: false
params.skip_variant = params.skip_variant ?: false
params.skip_barcode = params.skip_barcode ?: false
params.skip_snvcov  = params.skip_snvcov  ?: false

if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
}

if (params.platform != 'nanopore' && params.platform != 'pacbio' && params.platform != 'hifi') {
    exit 1, "Invalid protocol option: ${params.platform}. Valid options: 'nanopore', 'pacbio', 'hifi'"
}

if (!params.skip_variant) {
    if(params.skip_align) {
        exit 1, "Cannot run variant calling with --skip_align!"
    }
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
include { clair3_variant } from '../modules/clair3_variant.nf'
include { detect_barcode } from '../modules/detect_barcode.nf'
include { extract_snvs } from '../modules/extract_snvs.nf'

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow longread {
    check_inputs(ch_input)
    ch_sample = check_inputs.out.ch_cat_out
    ch_barcode = check_inputs.out.ch_barcode

    if (!params.skip_align) {
        minimap2_align(ch_sample)
        ch_bam = minimap2_align.out.ch_bam

        if (!params.skip_variant) {
            clair3_variant(ch_sample, ch_bam)
        }

        if (!params.skip_barcode) {
            ch_validated_barcode = ch_barcode.map { group, barcode_start, barcode_end, barcode_template ->
                if (!barcode_start || !barcode_end) { error "Error: barcode_start or barcode_end is empty. Barcode detection cannot proceed." }
                return tuple(group, barcode_start, barcode_end, barcode_template)
            }

            detect_barcode(ch_sample, ch_validated_barcode, ch_bam)
        }

        if (!params.skip_snvcov) {
            extract_snvs(ch_sample, ch_bam)
        }
    }
}
