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

def helpMessage()
{
    log.info """
    Usage:
    nextflow run check_inputs.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet        Path of the sample sheet
    
    Optional arguments:
        --outdir              the directory path of output results, default: the current directory
        --protocol            DNA, cDNA, directRNA, default: DNA
        --platform            nanopore, pacbio, hifi, default: nanopore
        --model               the trainning model of variant calling, default: ont_r10

    Step arguments:
        --skip_align          skip alignment
        --skip_variant        skip variant calling
    """
}

/*
#~~~~~~~~~~~~~~~~~~#
# check parameters #
#~~~~~~~~~~~~~~~~~~#
*/
params.help = null
if(params.help)
{
    helpMessage()
    exit 0
}

if (params.sample_sheet)
{
    ch_input = file(params.sample_sheet)
} else {
    helpMessage()
    println "========> Error!!! Please specify the full path of the sample sheet!\n"
    exit 1
}

params.outdir = params.outdir ?: "$PWD"

params.protocol = params.protocol ?: 'DNA'
if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA')
{
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
}

params.platform = params.platform ?: 'nanopore'
if (params.platform != 'nanopore' && params.platform != 'pacbio' && params.platform != 'hifi')
{
    exit 1, "Invalid protocol option: ${params.platform}. Valid options: 'nanopore', 'pacbio', 'hifi'"
}

params.skip_align = params.skip_align ?: false

params.skip_variant = params.skip_variant ?: false
if (!params.skip_variant)
{
    if(params.skip_align)
    {
        exit 1, "Cannot run variant calling with skipping alignment!"
    }
}

params.model = params.model ?: 'ont_r10'

/*
#~~~~~~~~~~~~~~#
# load modules #
#~~~~~~~~~~~~~~#
*/
include { check_inputs } from '../modules/check_inputs.nf'
include { minimap2_align } from '../modules/minimap2_align.nf'
include { clair3_variant } from '../modules/clair3_variant.nf'

/*
#~~~~~~~~~~#
# workflow #
#~~~~~~~~~~#
*/

workflow longread
{
    check_inputs(ch_input)
        .set{ ch_sample }

    if (!params.skip_align)
    {
        minimap2_align(ch_sample)
            .set{ ch_bam }

        if (!params.skip_variant)
        {
            clair3_variant(ch_sample, ch_bam)
        }
    }
}
