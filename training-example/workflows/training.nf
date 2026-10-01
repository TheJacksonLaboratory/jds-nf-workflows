#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { extract_csv } from '../bin/shared/extract_csv'
include { CREATE_MESSAGE } from '../modules/bash/create_message'
include { UPPERCASE_MESSAGE } from '../modules/bash/uppercase_message'
include { SUMMARIZE_MESSAGE } from '../modules/bash/summarize_message'
include { GATHER_SUMMARIES } from '../modules/bash/gather_summaries'

workflow TRAINING_WORKFLOW {
    main:
    input_samples = extract_csv(file(params.csv_input, checkIfExists: true))

    CREATE_MESSAGE(input_samples)
    UPPERCASE_MESSAGE(CREATE_MESSAGE.out.message_file)
    SUMMARIZE_MESSAGE(UPPERCASE_MESSAGE.out.uppercase_file)
    summary_files = SUMMARIZE_MESSAGE.out.summary_file
        .map { sampleID, summaryFile -> summaryFile }
        .collect()
    GATHER_SUMMARIES(summary_files)

    emit:
    summary_file = SUMMARIZE_MESSAGE.out.summary_file
    gathered_file = GATHER_SUMMARIES.out.gathered_file
}