#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {RSEM_PREPAREREFERENCE} from "${projectDir}/modules/rsem/rsem_preparereference"
include {RSEM_SAM_VALIDATOR} from "${projectDir}/modules/rsem/rsem_sam_validator"
include {CONVERT_SAM_FOR_RSEM} from "${projectDir}/modules/rsem/convert_sam_for_rsem"
include {RSEM_EXPRESSION} from "${projectDir}/modules/rsem/rsem_expression"
include {SEX_DETERMINATION} from "${projectDir}/modules/r/sex_determination"
include {MERGE_RSEM_COUNTS} from "${projectDir}/modules/utility_modules/merge_rsem_counts"
include {PICARD_ADDORREPLACEREADGROUPS} from "${projectDir}/modules/picard/picard_addorreplacereadgroups"
include {PICARD_REORDERSAM} from "${projectDir}/modules/picard/picard_reordersam"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_COLLECTRNASEQMETRICS} from "${projectDir}/modules/picard/picard_collectrnaseqmetrics"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

workflow RNA_FROM_BAM {

    take:
        bam_ch
        // [sampleID, bam]
    
    main:
        if (!params.rsem_reference_path || params.rsem_reference_path == '' || params.rsem_reference_path == null) {
            bowtie2_input = Channel.of([params.ref_fa, params.ref_gtf, 'bowtie2', ''])
            RSEM_PREPAREREFERENCE(bowtie2_input)

            rsem_ref_files = RSEM_PREPAREREFERENCE.out.all_files.toList()

            rsem_ref_name = file(params.ref_fa).baseName

        } else {

            rsem_ref_files = Channel
            .fromPath( "${params.rsem_reference_path}/**" )
            .collect()
            .toList()

            rsem_ref_name = params.rsem_reference_name
        
        }

        RSEM_SAM_VALIDATOR(bam_ch)
        
        bam_ch
        .join(RSEM_SAM_VALIDATOR.out.validation)
        .branch { it ->
            failed: it[2] == "FALSE"
                return [it[0], it[1]]
            passed: it[2] == "TRUE"
                return [it[0], it[1]]
        }
        .set { validation_result }

        CONVERT_SAM_FOR_RSEM(validation_result.failed)

        rsem_input = CONVERT_SAM_FOR_RSEM.out.rsem_bam
                    .mix(validation_result.passed)
                    .combine(rsem_ref_files)
                    .map{ it -> [it[0], it[1], it[2], rsem_ref_name] }

        RSEM_EXPRESSION(rsem_input)
}
