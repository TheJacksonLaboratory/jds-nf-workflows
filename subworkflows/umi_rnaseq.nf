#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {UMITOOLS_EXTRACT} from "${projectDir}/modules/umitools/umitools_extract"
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {GET_READ_LENGTH} from "${projectDir}/modules/utility_modules/get_read_length"
include {CHECK_STRANDEDNESS} from "${projectDir}/modules/python/python_check_strandedness"
include {STAR_ALIGN} from "${projectDir}/modules/star/star_align_rsem"
include {UMITOOLS_DEDUP as UMITOOLS_DEDUP_GENOME;
         UMITOOLS_DEDUP as UMITOOLS_DEDUP_TRANSCRIPT} from "${projectDir}/modules/umitools/umitools_dedup"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {UMITOOLS_PREPAREFORRSEM} from "${projectDir}/modules/umitools/umitools_prepareforrsem"
include {RSEM_EXPRESSION} from "${projectDir}/modules/rsem/rsem_expression_umi"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {SEX_DETERMINATION} from "${projectDir}/modules/r/sex_determination"
include {MERGE_RSEM_COUNTS} from "${projectDir}/modules/utility_modules/merge_rsem_counts"
include {PICARD_ADDORREPLACEREADGROUPS} from "${projectDir}/modules/picard/picard_addorreplacereadgroups"
include {PICARD_REORDERSAM} from "${projectDir}/modules/picard/picard_reordersam"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_COLLECTRNASEQMETRICS} from "${projectDir}/modules/picard/picard_collectrnaseqmetrics"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"


workflow UMI_RNASEQ {

    take:
        read_ch

    main:

        // UMI Extraction if needed
        ch_UMITOOLS_multiqc = Channel.empty() // optional log, depeding on skip umi extract
        if (!params.skip_umi_extract) {
            UMITOOLS_EXTRACT(read_ch)
            fastp_input = UMITOOLS_EXTRACT.out.umi_fastq
            ch_UMITOOLS_multiqc = UMITOOLS_EXTRACT.out.log
        } else {
            fastp_input = read_ch
        }

        // FASTP Quality Trimming
        ch_FASTP_multiqc = Channel.empty() // optional log, depeding on skip trim
        if (!params.skip_read_trimming) {
            FASTP(fastp_input)
            reads = FASTP.out.trimmed_fastq
            ch_FASTP_multiqc = FASTP.out.quality_json // set log file for multiqc
        } else {
            reads = read_ch
        }
        
        GET_READ_LENGTH(read_ch) // set to full reads, regardless of trim state.
        
        FASTQC(reads)

        // Check strand setting
        CHECK_STRANDEDNESS(reads)

        // STAR Alignment
        STAR_ALIGN(reads.join(GET_READ_LENGTH.out.read_length), params.rsem_ref_files, params.rsem_star_prefix)

        // UMI-Tools Deduplication
        UMITOOLS_DEDUP_GENOME(STAR_ALIGN.out.sorted_genomic_bam_bai)
        UMITOOLS_DEDUP_TRANSCRIPT(STAR_ALIGN.out.sorted_transcript_bam_bai)

        // Restore name sorting - RSEM rquires BAM to be name-sorted (i.e., query-name sorted)
        SAMTOOLS_SORT (
            UMITOOLS_DEDUP_TRANSCRIPT.out.bam,
            '-n',
            'name_sorted.bam'
        )

        // RSEM Quantification
        if (params.read_type == 'PE') {
            UMITOOLS_PREPAREFORRSEM(SAMTOOLS_SORT.out.sorted_file)
            rsem_input = UMITOOLS_PREPAREFORRSEM.out.bam.join(GET_READ_LENGTH.out.read_length).join(CHECK_STRANDEDNESS.out.strand_setting)
        } else {
            rsem_input = SAMTOOLS_SORT.out.sorted_file.join(GET_READ_LENGTH.out.read_length).join(CHECK_STRANDEDNESS.out.strand_setting)
        }

        RSEM_EXPRESSION(rsem_input, params.rsem_ref_files, params.rsem_star_prefix)
        
        if (params.gen_org == "human" || params.gen_org == "mouse") {
            SEX_DETERMINATION(RSEM_EXPRESSION.out.rsem_genes)
        }

        if (params.merge_rna_counts) {
            MERGE_RSEM_COUNTS(RSEM_EXPRESSION.out.rsem_genes.collect{it[1]},
                            RSEM_EXPRESSION.out.rsem_isoforms.collect{it[1]},
                            'allSamples')
        }

        // Get Read Group Information
        READ_GROUPS(reads, "picard")

        // Picard Alignment Metrics
        add_replace_groups = READ_GROUPS.out.read_groups.join(UMITOOLS_DEDUP_GENOME.out.bam)
        PICARD_ADDORREPLACEREADGROUPS(add_replace_groups)

        PICARD_REORDERSAM(PICARD_ADDORREPLACEREADGROUPS.out.bam, params.picard_dict)

        // Picard Alignment Metrics
        PICARD_SORTSAM(PICARD_REORDERSAM.out.bam, 'coordinate')
        
        PICARD_COLLECTRNASEQMETRICS(PICARD_SORTSAM.out.bam.join(CHECK_STRANDEDNESS.out.strand_setting), params.ref_flat, params.ribo_intervals)

        // Summary Stats
        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_UMITOOLS_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_FASTP_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CHECK_STRANDEDNESS.out.strandedness_report.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.star_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_DEDUP_TRANSCRIPT.out.log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(RSEM_EXPRESSION.out.rsem_cnt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS.out.picard_metrics.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        
}