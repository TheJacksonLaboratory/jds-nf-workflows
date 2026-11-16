#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/clonal_hema.nf"
include {param_log} from "../bin/log/clonal_hema.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {extract_csv} from "../bin/shared/extract_csv.nf"
include {CONCATENATE_LOCAL_FILES} from "../subworkflows/concatenate_local_files"

include {MM_CLONAL_HEMA} from "../subworkflows/mm_clonal_hema"
include {HS_CLONAL_HEMA} from "../subworkflows/hs_clonal_hema"

// Tumor-only somatic calling for a targeted CHIP gene panel, mouse + human,
// dispatched on `--gen_org` exactly like `--workflow pta` already is
// (pta.nf:41-45). Unlike pta, there is no tumor/normal pairing or
// proxy-normal substitution -- every sample sheet row is an independent
// tumor-only sample, filtered against a population germline-resource +
// panel-of-normals instead. See notes/germline_resource_design.md and the
// clonal_hema implementation plan for the full design reasoning.
//
// Sample sheet uses the SHARED `bin/shared/extract_csv.nf` schema
// (sampleID,fastq_1,fastq_2 + optional sex/lane/replicate/ind) -- not PTA's
// local extract_csv (patient,sampleID,sex,status), since that schema's
// patient/status columns exist only to drive pairing logic this workflow
// doesn't have.

workflow CLONAL_HEMA {

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

    if (!(params.min_alt_molecules.toString() ==~ /[1-9][0-9]*/)) {
        exit 1, "`--min_alt_molecules` must be a whole number of 1 or more (got: ${params.min_alt_molecules}). See `--help` for information."
    }

    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
    CONCATENATE_LOCAL_FILES(ch_input_sample)

    if (params.gen_org == "mouse") {

        MM_CLONAL_HEMA(CONCATENATE_LOCAL_FILES.out.read_meta_ch)

    } else if (params.gen_org == "human") {

        HS_CLONAL_HEMA(CONCATENATE_LOCAL_FILES.out.read_meta_ch)

    }

}
