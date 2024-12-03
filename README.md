<div align="center">
<h1 align="center">Long Read Processing Pipeline</h1>
  <p align="center">A nextflow pipeline of processing long reads</p>
</div>

## Table of Contents
<details open>
<summary><b>[Show or Hide]</b></summary>

1. [Dependencies](#dependencies)
2. [File Format](#file-format)
    - [Structure of input directories](#structure)
    - [Sample sheet](#samplesheet)
3. [Usage](#usage)
    - [Run job](#runjob)
    - [Usage options](#options)
</details>

<!-- Dependencies-->
## Dependencies
* nextflow
* minimap2
* samtools
* bamtools
* python packages
    - pysam
    - biopython
    - numpy
    - scipy
    - matplotlib
    - pandas 

<!-- File Format-->
## File Format

<a id="structure"></a>

### Structure of input directories
![example](./image/inputs.png)

<a id="samplesheet"></a>

### Sample Sheet -- csv
| group | barcode_start | barcode_end | barcode_template | directory | reference | gene_info |
| - | - | - | - | - | - | - |
| bc2081 | 5832 | 5869 | NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN | /path/of/directory/ | GUK1_PGJJ162.fa | bc2081.txt |
| bc2082 | 4819 | 4856 | NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN | /path/of/directory/ | GUK1_PTB198.fa | bc2082.txt |
| bc2083 | 4791 | 4828 | NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN | /path/of/directory/ | GUK1_PH003.fa | bc2083.txt |

> [!Note]  
> 1. The sample sheet must be a csv file and the header must be like below in the example
> 2. For barcode association, you need everything in the field
> 3. For variant coverer plot, you only need group, directory, reference

### Gene Info -- tsv

| | | |
| - | - | - |
| hGUK1_wt | 6256 | 591 |
| yGUK1_wt | 6256 | 561 |
| pGUK1_wt | 6256 | 597 |

> [!Note]  
> 1. The gene info file must be a tsv file and in the directory according to the sample sheet
> 2. The gene info file has no header
> 3. The first column is the reference ID corresponding to the gene (check your reference fasta file)
> 4. The second column is the start position of the gene in the reference
> 5. The third column is the length of the gene

<!-- Usage-->
## Usage

<a id="runjob"></a>

### Run job
submit the bash script below

```bash
#!/bin/bash
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -R "select[mem>1000] rusage[mem=1000]"
#BSUB -M 1000
#BSUB -q normal

# modules
module load HGI/common/nextflow/23.10.0
module load HGI/softpack/users/fs18/nf_longread

#--------------#
# user specify #
#--------------#
# LSF group
export LSB_DEFAULT_USERGROUP=hgi

# Paths
export INPUTSAMPLE=$PWD/inputs/samplesheet.csv
export OUTPUTRES=$PWD/outputs

#-----------#
# pipelines #
#-----------#
nextflow run -resume nf_longread/main.nf --sample_sheet $INPUTSAMPLE \
                                         --protocol DNA \
                                         --platform hifi \
                                         --outdir $OUTPUTRES \
                                         --skip_snvcov
```

<a id="options"></a>

### Usage options
```bash
nextflow run check_inputs.nf --sample_sheet "/path/of/sample/sheet"

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
```