#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/joint_gvcf_calling.nf"
include {param_log} from "${projectDir}/bin/log/joint_gvcf_calling.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv_gvcf.nf"
include {SAMTOOLS_FAIDX_TO_BED} from "${projectDir}/modules/samtools/samtools_faidx_to_bed"
include {GATK_SPLITINTERVALS} from "${projectDir}/modules/gatk/gatk_splitintervals"
include {GATK_INDEXFEATUREFILE} from "${projectDir}/modules/gatk/gatk_indexfeaturefile"
include {GATK_GENOMICSDBIMPORT} from "${projectDir}/modules/gatk/gatk_genomicsdbimport"
include {GATK_GENOTYPEGVCF} from "${projectDir}/modules/gatk/gatk_genotypegvcf"
include {GATK_GATHERVCFS} from "${projectDir}/modules/gatk/gatk_gathervcfs"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if (params.csv_input) {
    def csvFile = file(params.csv_input)
    if (!csvFile.exists()) {
        exit 1, "ERROR: CSV input file does not exist. Parameter csv_input was set to: ${params.csv_input}"
    }
    ch_input_sample = extract_csv(csvFile)
    ch_input_sample.map{it -> [it[0], it[2]]}.set{gvcf_ch}
    ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
} else {
    exit 1, "ERROR: CSV input is required. Parameter csv_input was set to: ${params.csv_input}"
}

// main workflow
workflow JOINT_GVCF_CALLING {

    SAMTOOLS_FAIDX_TO_BED()
    GATK_SPLITINTERVALS(SAMTOOLS_FAIDX_TO_BED.out.bed)

    intervals = GATK_SPLITINTERVALS.out.intervals
                .flatten()

    GATK_INDEXFEATUREFILE(gvcf_ch)

    combine_input_ch = gvcf_ch
    .join(GATK_INDEXFEATUREFILE.out.idx) 
    .join(meta_ch)
    .map{ it -> [it[3]['output'], it[1], it[2]] }
    .groupTuple()
    .combine(intervals)

    GATK_GENOMICSDBIMPORT(combine_input_ch)

    GATK_GENOTYPEGVCF(GATK_GENOMICSDBIMPORT.out.genomicsdb)

    GATK_GATHERVCFS(GATK_GENOTYPEGVCF.out.vcf.groupTuple(size: params.interval_count), 'Genotyped_allChroms')
}
