#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_wes_simreads"
include {param_log} from "${projectDir}/bin/log/generate_wes_simreads"
include {final_run_report} from "${projectDir}/bin/shared/final_run_report.nf"
include {SAMTOOLS_FAIDX_CHR_ONLY} from "${projectDir}/modules/samtools/samtools_faidx_chr_only"
include {SPLIT_FILES} from "${projectDir}/modules/python/pyfaidx_split_files"
include {GENERATE_SIMULATED_WES_DATA_FIRST_PASS} from "${projectDir}/modules/neat/generate_simulated_WES_data_first_pass"
include {GENERATE_SIMULATED_WES_DATA_SECOND_PASS} from "${projectDir}/modules/neat/generate_simulated_WES_data_second_pass"
include {COMPRESS_INDEX_VCF_SIMREADS as COMPRESS_FIRST_PASS;
         COMPRESS_INDEX_VCF_SIMREADS as COMPRESS_SECOND_PASS} from "${projectDir}/modules/tabix/compress_vcf_simreads"
include {BCFTOOLS_MERGE_VCF;
         BCFTOOLS_MERGE_VCF as BCFTOOLS_MERGE_VCF_SECOND_PASS} from "${projectDir}/modules/bcftools/bcftools_merge_vcf"
include {VCF_SORT} from "${projectDir}/modules/vcftools/vcf_sort"
include {VCFTOOLS_SIMVAR} from "${projectDir}/modules/vcftools/vcftools_simvar"
include {CONCATENATE_SIMREADS} from "${projectDir}/modules/utility_modules/concatenate_simreads"


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


workflow GENERATE_WES_SIMREADS {

    SAMTOOLS_FAIDX_CHR_ONLY(params.fasta)
    SPLIT_FILES(SAMTOOLS_FAIDX_CHR_ONLY.out.chr_fa)
    GENERATE_SIMULATED_WES_DATA_FIRST_PASS(SPLIT_FILES.out.flatten()) // flatten is the KEY 
    
    COMPRESS_FIRST_PASS(GENERATE_SIMULATED_WES_DATA_FIRST_PASS.out.vcf)

    vcf_ch = COMPRESS_FIRST_PASS.out.vcf.collect()
    tbi_ch = COMPRESS_FIRST_PASS.out.tbi.collect()

    BCFTOOLS_MERGE_VCF(vcf_ch, tbi_ch, 'ALL')
    VCF_SORT(BCFTOOLS_MERGE_VCF.out.merged_vcf)
    VCFTOOLS_SIMVAR(VCF_SORT.out.vcf, VCF_SORT.out.tbi)    

    GENERATE_SIMULATED_WES_DATA_SECOND_PASS(SPLIT_FILES.out.flatten(), VCFTOOLS_SIMVAR.out.int_vcf )
    COMPRESS_SECOND_PASS(GENERATE_SIMULATED_WES_DATA_SECOND_PASS.out.vcf)

    vcf2_ch = COMPRESS_SECOND_PASS.out.vcf.collect()
    tbi2_ch = COMPRESS_SECOND_PASS.out.tbi.collect()

    BCFTOOLS_MERGE_VCF_SECOND_PASS(vcf2_ch, tbi2_ch, 'FINAL_ALL')

    ch_sampleID = params.sampleID ? Channel.value(params.sampleID) : null
    ch_fastq1 = GENERATE_SIMULATED_WES_DATA_SECOND_PASS.out.fq1.collect()
    ch_fastq2 = GENERATE_SIMULATED_WES_DATA_SECOND_PASS.out.fq2.collect()


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
