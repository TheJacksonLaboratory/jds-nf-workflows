#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/mitochondria_variant_calling.nf"
include {param_log} from "${projectDir}/bin/log/mitochondria_variant_calling.nf"
include {extract_csv_bam} from "${projectDir}/bin/shared/extract_csv_bam.nf"
include {MT_VARIANT_CALLING} from "${projectDir}/subworkflows/mt_variant_calling"


// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

bam_ch = extract_csv_bam(file(params.csv_input, checkIfExists: true))
input_ch = bam_input_ch.map{it -> [it[0], file(it[2]), file(it[3])]}

// main workflow
workflow MITOCHONDRIA_VARIANT_CALLING {
    MT_VARIANT_CALLING(input_ch)
    // workflow found in: subworkflows/mt_variant_calling.nf
    // workflow run as subworkflow due to re-use in WGS and PTA (possible) workflows. 
}
