#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/chip_somatic.nf"
include {param_log} from "../bin/log/chip_somatic.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {extract_csv} from "../bin/shared/extract_csv.nf"
include {CONCATENATE_LOCAL_FILES} from "../subworkflows/concatenate_local_files"

include {MM_CHIP_SOMATIC} from "../subworkflows/mm_chip_somatic"
include {HS_CHIP_SOMATIC} from "../subworkflows/hs_chip_somatic"

// Tumor-only somatic calling for a targeted CHIP gene panel, mouse + human,
// dispatched on `--gen_org` exactly like `--workflow pta` already is
// (pta.nf:41-45). Unlike pta, there is no tumor/normal pairing or
// proxy-normal substitution -- every sample sheet row is an independent
// tumor-only sample, filtered against a population germline-resource +
// panel-of-normals instead. See notes/germline_resource_design.md and the
// chip_somatic implementation plan for the full design reasoning.
//
// Sample sheet uses the SHARED `bin/shared/extract_csv.nf` schema
// (sampleID,fastq_1,fastq_2 + optional sex/lane/replicate/ind) -- not PTA's
// local extract_csv (patient,sampleID,sex,status), since that schema's
// patient/status columns exist only to drive pairing logic this workflow
// doesn't have.

workflow CHIP_SOMATIC {

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

    if (!params.csv_input) {
        exit 1, "No input CSV file was specified with `--csv_input`. A CSV manifest is required. See `--help` for information."
    }

    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
    CONCATENATE_LOCAL_FILES(ch_input_sample)

    if (params.gen_org == "mouse") {

        MM_CHIP_SOMATIC(CONCATENATE_LOCAL_FILES.out.read_meta_ch)

    } else if (params.gen_org == "human") {

        HS_CHIP_SOMATIC(CONCATENATE_LOCAL_FILES.out.read_meta_ch)

    }

}
