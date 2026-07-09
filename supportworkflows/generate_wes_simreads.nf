#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "../bin/help/generate_wes_simreads"
include {param_log} from "../bin/log/generate_wes_simreads"
include {final_run_report} from "../bin/shared/final_run_report.nf"
include {SAMTOOLS_FAIDX_CHR_ONLY} from "../modules/samtools/samtools_faidx_chr_only"
include {SPLIT_FILES} from "../modules/python/pyfaidx_split_files"
include {GENERATE_SIMULATED_WES_DATA} from "../modules/neat/generate_simulated_WES_data"
include {COMPRESS_INDEX_VCF_SIMREADS} from "../modules/tabix/compress_vcf_simreads"
include {BCFTOOLS_MERGE_VCF} from "../modules/bcftools/bcftools_merge_vcf"
include {VCF_SORT} from "../modules/vcftools/vcf_sort"
include {VCFTOOLS_SIMVAR} from "../modules/vcftools/vcftools_simvar"
include {CONCATENATE_READS_SE} from "../modules/utility_modules/concatenate_reads_SE"
include {CONCATENATE_READS_PE} from "../modules/utility_modules/concatenate_reads_PE"


workflow GENERATE_WES_SIMREADS {
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

    SAMTOOLS_FAIDX_CHR_ONLY(params.fasta)
    SPLIT_FILES(SAMTOOLS_FAIDX_CHR_ONLY.out.chr_fa)
    GENERATE_SIMULATED_WES_DATA(SPLIT_FILES.out.flatten(), params.target_bed) // flatten is the KEY, as it converts the channel of lists into a channel of individual files, which is what the process expects as input
    
    COMPRESS_INDEX_VCF_SIMREADS(GENERATE_SIMULATED_WES_DATA.out.vcf)

    vcf_ch = COMPRESS_INDEX_VCF_SIMREADS.out.vcf.collect()
    tbi_ch = COMPRESS_INDEX_VCF_SIMREADS.out.tbi.collect()

    BCFTOOLS_MERGE_VCF(vcf_ch, tbi_ch)
    VCF_SORT(BCFTOOLS_MERGE_VCF.out.merged_vcf) 
    VCFTOOLS_SIMVAR(VCF_SORT.out.vcf, VCF_SORT.out.tbi, params.target_bed) // restrict callset to just target region (i.e., remove off target calls).

    if (params.read_type == 'PE') {
        ch_fastq1 = GENERATE_SIMULATED_WES_DATA.out.fq1.collect()
        ch_fastq2 = GENERATE_SIMULATED_WES_DATA.out.fq2.collect()

        reads1_with_id = ch_fastq1.map { files -> [params.sampleID, files] }
        reads2_with_id = ch_fastq2.map { files -> [params.sampleID, files] }

        fq_reads_pe = reads1_with_id.join(reads2_with_id)

        CONCATENATE_READS_PE(fq_reads_pe)
    } else {
        ch_fastq1 = GENERATE_SIMULATED_WES_DATA.out.fq1.collect()

        reads1_with_id = ch_fastq1.map { files -> [params.sampleID, files] }

        CONCATENATE_READS_SE(reads1_with_id)
    }

}

def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}
