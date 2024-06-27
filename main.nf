#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { longread } from './workflows/longread.nf'

workflow {
    longread()
}
