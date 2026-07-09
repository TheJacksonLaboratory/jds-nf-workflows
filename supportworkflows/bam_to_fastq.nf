#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {help} from "../bin/help/wgs.nf"
include {param_log} from "../bin/log/bam_to_fastq.nf"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {extract_csv_bam_rnaseq} from "../bin/shared/extract_csv_bam.nf"
// note, this is used because it does not require a BAI file. The extract_csv_bam function requires a BAI file, and will throw an error if it is not provided.
include {SAMTOOLS_VIEW} from "../modules/samtools/samtools_view"
include {SAMTOOLS_SORT} from "../modules/samtools/samtools_sort"
include {SAMTOOLS_FASTQ} from "../modules/samtools/samtools_fastq"

// BEGIN main workflow
workflow BAM_TO_FASTQ {

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
        error("Missing required parameter: --csv_input. Please provide a CSV file path.")
    }

    def csv_file = file(params.csv_input)
    if (!csv_file.exists()) {
        error("CSV input file does not exist: ${params.csv_input}")
    }

    if (!csv_file.isFile()) {
        error("CSV input path is not a file: ${params.csv_input}")
    }


    bam_input_ch = extract_csv_bam_rnaseq(csv_file)
    bam_file = bam_input_ch.map{it -> [it[0], file(it[1])]} 

    SAMTOOLS_VIEW(bam_file, '-b -F 2304', 'F_2304')
    SAMTOOLS_SORT(SAMTOOLS_VIEW.out.bam, '-O bam -n', 'bam')
    SAMTOOLS_FASTQ(SAMTOOLS_SORT.out.sorted_file)
}
