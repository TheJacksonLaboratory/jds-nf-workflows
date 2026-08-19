#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules

include {PACBIO} from "../subworkflows/pacbio"
include {ILLUMINA} from "../subworkflows/illumina"
include {ONT} from "../subworkflows/ont"

// main workflow
workflow GERMLINE_SV {
  if (params.data_type == 'pacbio') {
    PACBIO()
  } else if (params.data_type == 'illumina') {
    ILLUMINA()
  } else if (params.data_type == 'ont') {
    ONT()
  } else {
    // if it is not a validate data type, an acceptable string
    exit 1, "'--data_type': \"${params.data_type}\" is not valid, supported options are 'pacbio' or 'illumina' or 'ont'" 
  }
}
