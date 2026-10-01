#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRAINING_WORKFLOW } from './workflows/training'

workflow {
    TRAINING_WORKFLOW()
}