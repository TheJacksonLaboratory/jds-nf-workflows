#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_wgs_simreads"
include {param_log} from "${projectDir}/bin/log/generate_wgs_simreads"
include {final_run_report} from "${projectDir}/bin/shared/final_run_report.nf"
include {SAMTOOLS_FAIDX_CHR_ONLY} from "${projectDir}/modules/samtools/samtools_faidx_chr_only"
include {SPLIT_FILES} from "${projectDir}/modules/python/pyfaidx_split_files"
include {GENERATE_SIMULATED_WGS_DATA} from "${projectDir}/modules/neat/generate_simulated_WGS_data"
include {COMPRESS_INDEX_VCF_SIMREADS} from "${projectDir}/modules/tabix/compress_vcf_simreads"
include {BCFTOOLS_MERGE_VCF} from "${projectDir}/modules/bcftools/bcftools_merge_vcf"
include {VCF_SORT} from "${projectDir}/modules/vcftools/vcf_sort"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"


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

def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}

workflow GENERATE_WGS_SIMREADS {

    SAMTOOLS_FAIDX_CHR_ONLY(params.fasta)
    SPLIT_FILES(SAMTOOLS_FAIDX_CHR_ONLY.out.chr_fa)
    GENERATE_SIMULATED_WGS_DATA(SPLIT_FILES.out.flatten()) // flatten is the KEY, as it converts the channel of lists into a channel of individual files, which is what the process expects as input

    COMPRESS_INDEX_VCF_SIMREADS(GENERATE_SIMULATED_WGS_DATA.out.vcf)

    vcf_ch = COMPRESS_INDEX_VCF_SIMREADS.out.vcf.collect()
    tbi_ch = COMPRESS_INDEX_VCF_SIMREADS.out.tbi.collect()
    
    BCFTOOLS_MERGE_VCF(vcf_ch, tbi_ch) 
    VCF_SORT(BCFTOOLS_MERGE_VCF.out.merged_vcf)

    if (params.read_type == 'PE') {
        CONCATENATE_READS_PE([params.sampleID, GENERATE_SIMULATED_WGS_DATA.out.fq1.collect(), GENERATE_SIMULATED_WGS_DATA.out.fq2.collect()])
    } else {
        CONCATENATE_READS_SE([params.sampleID, GENERATE_SIMULATED_WGS_DATA.out.fq1.collect()])
    }

}

