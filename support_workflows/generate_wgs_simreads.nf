#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_wgs_simreads"
include {param_log} from "${projectDir}/bin/log/generate_wgs_simreads"
include {SAMTOOLS_FAIDX_CHR_ONLY} from "${projectDir}/modules/samtools/samtools_faidx_chr_only"
include {SPLIT_FILES} from "${projectDir}/modules/python/pyfaidx_split_files"
include {GENERATE_SIMULATED_WGS_DATA} from "${projectDir}/modules/neat/generate_simulated_WGS_data"
include {COMPRESS_INDEX_VCF_SIMREADS} from "${projectDir}/modules/tabix/compress_vcf_simreads"
include {BCFTOOLS_MERGE_VCF} from "${projectDir}/modules/bcftools/bcftools_merge_vcf"
include {VCF_SORT} from "${projectDir}/modules/vcftools/vcf_sort"
include {CONCATENATE_SIMREADS} from "${projectDir}/modules/utility_modules/concatenate_simreads"


// help if needed
if (params.help){
    help()
    exit 0
}

// log params
message = param_log()


def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}


workflow GENERATE_WGS_SIMREADS {

    SAMTOOLS_FAIDX_CHR_ONLY(params.fasta)
    SPLIT_FILES(SAMTOOLS_FAIDX_CHR_ONLY.out.chr_fa)
    GENERATE_SIMULATED_WGS_DATA(SPLIT_FILES.out.flatten()) // flatten is the KEY  

    COMPRESS_INDEX_VCF_SIMREADS(GENERATE_SIMULATED_WGS_DATA.out.vcf)

    vcf_ch = COMPRESS_INDEX_VCF_SIMREADS.out.vcf.collect()
    tbi_ch = COMPRESS_INDEX_VCF_SIMREADS.out.tbi.collect()
    
    BCFTOOLS_MERGE_VCF(vcf_ch, tbi_ch) 
    VCF_SORT(BCFTOOLS_MERGE_VCF.out.merged_vcf)


    ch_sampleID = params.sampleID ? Channel.value(params.sampleID) : null
    ch_fastq1 = GENERATE_SIMULATED_WGS_DATA.out.fq1.collect()
    ch_fastq2 = GENERATE_SIMULATED_WGS_DATA.out.fq2.collect()


    reads1_with_id = ch_fastq1.map { file -> 
          def id = params.sampleID // Extract 'sample1'
          return [id, file]
    }

    reads2_with_id = ch_fastq2.map { file ->
          def id = params.sampleID // Extract 'sample1'
          return [id, file]
    }

    fq_reads = reads1_with_id.join(reads2_with_id)
    CONCATENATE_SIMREADS(fq_reads)

}

