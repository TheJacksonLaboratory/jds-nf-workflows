#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {GATK_PRINTREADS} from "${projectDir}/modules/gatk/gatk_printreads_mt"
include {PICARD_REVERTSAM} from "${projectDir}/modules/picard/picard_revertsam"
include {PICARD_SAMTOFASTQ} from "${projectDir}/modules/picard/picard_samtofastq"
include {BWA_MEM as BWA_MEM_MT;
         BWA_MEM as BWA_MEM_SHIFTED_MT} from "${projectDir}/modules/bwa/bwa_mem_mt"
include{PICARD_MERGEBAMALIGNMENT as PICARD_MERGEBAMALIGNMENT_MT;
        PICARD_MERGEBAMALIGNMENT as PICARD_MERGEBAMALIGNMENT_SHIFTED_MT} from "${projectDir}/modules/picard/picard_mergebamalignment"
include {PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT;
         PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_SHIFTED_MT} from "${projectDir}/modules/picard/picard_markduplicates_mt"
include {PICARD_SORTSAM as PICARD_SORTSAM_MT;
         PICARD_SORTSAM as PICARD_SORTSAM_SHIFTED_MT} from "${projectDir}/modules/picard/picard_sortsam"

include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"

include {GATK_MUTECT2_MT as GATK_MUTECT2_MT;
         GATK_MUTECT2_MT as GATK_MUTECT2_SHIFTEDMT} from "${projectDir}/modules/gatk/gatk_mutect2_mt"

include {PICARD_LIFTOVERVCF_MERGEVCF} from "${projectDir}/modules/picard/picard_liftovervcf_mergevcf"

include {GATK_MERGEMUTECTSTATS} from "${projectDir}/modules/gatk/gatk_mergemutectstats"

include {GATK_FILTERMUECTCALLS;
         GATK_FILTERMUECTCALLS as GATK_FILTERMUECTCALLS_PRIMARY} from "${projectDir}/modules/gatk/gatk_filtermutectcalls_mt"

include {GATK_LEFTALIGNANDTRIMVARIANTS;
         GATK_LEFTALIGNANDTRIMVARIANTS as GATK_LEFTALIGNANDTRIMVARIANTS_PASS} from "${projectDir}/modules/gatk/gatk_leftalignandtrimvariants"

include {HAPLOCHECK} from "${projectDir}/modules/haplocheck/haplocheck"

include {PICARD_MT_COVERAGEATEVERYBASE} from "${projectDir}/modules/picard/picard_mt_coverageateverybase"

workflow MT_VARIANT_CALLING {

    take:
        bam_ch
        // this is a tuple with [sampleID, bam, bai]

    main:
        // Extract reads from the mitochondrial contig
        GATK_PRINTREADS(bam_ch)

        PICARD_REVERTSAM(GATK_PRINTREADS.out.bam)

        PICARD_SAMTOFASTQ(PICARD_REVERTSAM.out.bam)

        BWA_MEM_MT(PICARD_SAMTOFASTQ.out.fastq.map { fastq -> tuple(fastq[0], fastq[1], 'mt') }, params.mt_fasta_index)
        BWA_MEM_SHIFTED_MT(PICARD_SAMTOFASTQ.out.fastq.map { fastq -> tuple(fastq[0], fastq[1], 'shifted_mt') }, params.mt_shifted_fasta_index)

        PICARD_MERGEBAMALIGNMENT_MT(BWA_MEM_MT.out.sam.join(PICARD_REVERTSAM.out.bam), 'mt')
        PICARD_MERGEBAMALIGNMENT_SHIFTED_MT(BWA_MEM_SHIFTED_MT.out.sam.join(PICARD_REVERTSAM.out.bam), 'shifted_mt')

        PICARD_MARKDUPLICATES_MT(PICARD_MERGEBAMALIGNMENT_MT.out.bam)
        PICARD_MARKDUPLICATES_SHIFTED_MT(PICARD_MERGEBAMALIGNMENT_SHIFTED_MT.out.bam)

        PICARD_SORTSAM_MT(PICARD_MARKDUPLICATES_MT.out.dedup_bam, 'coordinate')
        PICARD_SORTSAM_SHIFTED_MT(PICARD_MARKDUPLICATES_SHIFTED_MT.out.dedup_bam, 'coordinate')

        mt_alignment = PICARD_SORTSAM_MT.out.bam.join(PICARD_SORTSAM_MT.out.bai)
        shifted_mt_alignment = PICARD_SORTSAM_SHIFTED_MT.out.bam.join(PICARD_SORTSAM_SHIFTED_MT.out.bai)

        PICARD_COLLECTWGSMETRICS(PICARD_SORTSAM_MT.out.bam, 'mt')
        // This needs to be checked for cases where it is wgs and not mt to make sure it works correctly.  

        if (params.gen_org == 'human') {

            GATK_MUTECT2_MT(mt_alignment, 'mt', 'chrM:576-16024')
            GATK_MUTECT2_SHIFTEDMT(shifted_mt_alignment, 'shifted_mt', 'chrM:8024-9145')
            
            merge_stats_input = GATK_MUTECT2_MT.out.stats.mix(GATK_MUTECT2_SHIFTEDMT.out.stats).groupTuple(size: 2)
            GATK_MERGEMUTECTSTATS(merge_stats_input)

            vcf_for_join = GATK_MUTECT2_MT.out.vcf.join(GATK_MUTECT2_SHIFTEDMT.out.vcf)

            PICARD_LIFTOVERVCF_MERGEVCF(vcf_for_join)

            GATK_FILTERMUECTCALLS_PRIMARY(PICARD_LIFTOVERVCF_MERGEVCF.out.vcf.join(GATK_MERGEMUTECTSTATS.out.stats).map{it -> it[0], it[1], it[2], 'primary'})

            GATK_LEFTALIGNANDTRIMVARIANTS_PASS(GATK_FILTERMUECTCALLS_PRIMARY.out.vcf.join(GATK_FILTERMUECTCALLS_PRIMARY.out.tbi), 'pass-only')

            HAPLOCHECK(GATK_LEFTALIGNANDTRIMVARIANTS_PASS.out.vcf)

            GATK_FILTERMUECTCALLS(GATK_FILTERMUECTCALLS_PRIMARY.out.vcf.join(GATK_MERGEMUTECTSTATS.out.stats).join(HAPLOCHECK.out.contamination))
            
            // MITY?

        } else if (params.gen_org == 'mouse') {

            GATK_MUTECT2_MT(mt_alignment, 'mt', 'MT:1-15423')
            GATK_MUTECT2_SHIFTEDMT(shifted_mt_alignment, 'shifted_mt', 'MT:7423-8299')

            merge_stats_input = GATK_MUTECT2_MT.out.stats.mix(GATK_MUTECT2_SHIFTEDMT.out.stats).groupTuple(size: 2)
            GATK_MERGEMUTECTSTATS(merge_stats_input)

            vcf_for_join = GATK_MUTECT2_MT.out.vcf.join(GATK_MUTECT2_SHIFTEDMT.out.vcf)

            PICARD_LIFTOVERVCF_MERGEVCF(vcf_for_join)

            filter_calls_input = PICARD_LIFTOVERVCF_MERGEVCF.out.vcf.join(GATK_MERGEMUTECTSTATS.out.stats).map{it -> [it[0], it[1], it[2], '/mouse']}
            GATK_FILTERMUECTCALLS(filter_calls_input)
            // Note: 'mouse' is a placeholder for path(contamination), which is not used in mouse
            
            // MITY?

        }
        // for human the original positions are: chrM:576-16024 and the shifted calling is done in chrM:8024-9145.
        // for mouse the original positions are: MT:1-15423, and the shifted calling is done in MT:7423-8299.
        // Note that GRCm38 and GRCm39 use the same MT coordinates, so no need to check the genome version.

        PICARD_MT_COVERAGEATEVERYBASE(mt_alignment.join(shifted_mt_alignment))

        GATK_LEFTALIGNANDTRIMVARIANTS(GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 'all-calls')

    emit:
        output_string = "nothing yet" // This is a placeholder;
}

// Punch list: 

// Add MITY variant caller
// Check Picard CollectWgsMetrics works for WGS
// Add annotation
// Wiki page
// Emit strings and connection to WGS and PTA workflows
// Test human workflow
// Write tests and add test data

// quay.io/biocontainers/mity:2.0.0--pyhdfd78af_0

//   output {
//     File subset_bam = SubsetBamToChrM.output_bam
//     File subset_bai = SubsetBamToChrM.output_bai
//     File mt_aligned_bam = AlignAndCall.mt_aligned_bam
//     File mt_aligned_bai = AlignAndCall.mt_aligned_bai
//     File out_vcf = AlignAndCall.out_vcf
//     File out_vcf_index = AlignAndCall.out_vcf_index
//     File split_vcf = SplitMultiAllelicSites.split_vcf
//     File split_vcf_index = SplitMultiAllelicSites.split_vcf_index
//     File input_vcf_for_haplochecker = AlignAndCall.input_vcf_for_haplochecker
//     File duplicate_metrics = AlignAndCall.duplicate_metrics
//     File coverage_metrics = AlignAndCall.coverage_metrics
//     File theoretical_sensitivity_metrics = AlignAndCall.theoretical_sensitivity_metrics
//     File contamination_metrics = AlignAndCall.contamination_metrics
//     File base_level_coverage_metrics = CoverageAtEveryBase.table
//     Int mean_coverage = AlignAndCall.mean_coverage
//     Float median_coverage = AlignAndCall.median_coverage
//     String major_haplogroup = AlignAndCall.major_haplogroup
//     Float contamination = AlignAndCall.contamination
//   }