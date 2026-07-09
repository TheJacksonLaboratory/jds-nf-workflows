#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/wgs_sv_bam.nf"
include {param_log} from "../bin/log/wgs_sv_bam.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {extract_csv_bam} from "../bin/shared/extract_csv_bam.nf"
include {WGS_SV} from "../subworkflows/wgs_sv"

// main workflow
workflow WGS_SV_BAM {

    // help if needed
    if (params.help){
        help()
        exit 0
    }

    // log params
    message = param_log()

    // Save params to a file for record-keeping
    workflow.onComplete {
        final_run_report(message)
    }

    bam_ch = extract_csv_bam(file(params.csv_input, checkIfExists: true))

    bam_ch = bam_ch.map{ sampleID, meta, bam, bai -> [sampleID, bam, bai] }
    WGS_SV(bam_ch)
    // workflow found in: subworkflows/wgs_sv.nf
}
