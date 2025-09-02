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

include {MITY_RUNALL} from "${projectDir}/modules/mity/mity_runall"
include {MUTSERVE} from "${projectDir}/modules/mutserve/mutserve"

include {PREP_MTDNA_VCF} from "${projectDir}/modules/r/prep_mtdna_vcf"
include {BCFTOOLS_MERGECALLERS} from "${projectDir}/modules/bcftools/bcftools_merge_mt_callers"
include {COMPUTE_AVG_AF} from "${projectDir}/modules/r/compute_avg_af"

include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_COSMIC} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_DBNSFP} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"

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

        MITY_RUNALL(mt_alignment)

        MUTSERVE(mt_alignment)

        PICARD_COLLECTWGSMETRICS(PICARD_SORTSAM_MT.out.bam, 'mt')
        // This needs to be checked for cases where it is wgs and not mt to make sure it works correctly.  

        if (params.gen_org == 'human') {

            GATK_MUTECT2_MT(mt_alignment, 'mt', 'chrM:576-16024')
            GATK_MUTECT2_SHIFTEDMT(shifted_mt_alignment, 'shifted_mt', 'chrM:8024-9145')
            
            merge_stats_input = GATK_MUTECT2_MT.out.stats.mix(GATK_MUTECT2_SHIFTEDMT.out.stats).groupTuple(size: 2)
            GATK_MERGEMUTECTSTATS(merge_stats_input)

            vcf_for_join = GATK_MUTECT2_MT.out.vcf.join(GATK_MUTECT2_SHIFTEDMT.out.vcf)

            PICARD_LIFTOVERVCF_MERGEVCF(vcf_for_join)

            primary_filter_input = PICARD_LIFTOVERVCF_MERGEVCF.out.vcf
                                   .join(GATK_MERGEMUTECTSTATS.out.stats)
                                   .map{ it -> [it[0], it[1], it[2], '/primary'] }

            GATK_FILTERMUECTCALLS_PRIMARY(primary_filter_input, 'primary')

            GATK_LEFTALIGNANDTRIMVARIANTS_PASS(GATK_FILTERMUECTCALLS_PRIMARY.out.vcf_tbi, 'pass-only')

            HAPLOCHECK(GATK_LEFTALIGNANDTRIMVARIANTS_PASS.out.interm_vcf_tbi)

            contam_filter_input = GATK_FILTERMUECTCALLS_PRIMARY.out.vcf_tbi
                                  .join(GATK_MERGEMUTECTSTATS.out.stats)
                                  .join(HAPLOCHECK.out.contam_value)
                                  .map{it -> [it[0], it[1], it[3], it[4]]}

            GATK_FILTERMUECTCALLS(contam_filter_input, 'final')

        } else if (params.gen_org == 'mouse') {

            GATK_MUTECT2_MT(mt_alignment, 'mt', 'MT:1-15423')
            GATK_MUTECT2_SHIFTEDMT(shifted_mt_alignment, 'shifted_mt', 'MT:7423-8299')

            merge_stats_input = GATK_MUTECT2_MT.out.stats.mix(GATK_MUTECT2_SHIFTEDMT.out.stats).groupTuple(size: 2)
            GATK_MERGEMUTECTSTATS(merge_stats_input)

            vcf_for_join = GATK_MUTECT2_MT.out.vcf.join(GATK_MUTECT2_SHIFTEDMT.out.vcf)

            PICARD_LIFTOVERVCF_MERGEVCF(vcf_for_join)

            filter_calls_input = PICARD_LIFTOVERVCF_MERGEVCF.out.vcf
                                .join(GATK_MERGEMUTECTSTATS.out.stats)
                                .map{it -> [it[0], it[1], it[2], '/mouse']}
            
            GATK_FILTERMUECTCALLS(filter_calls_input, 'final')
            // Note: 'mouse' is a placeholder for path(contamination), which is not used in mouse

        }
        // for human the original positions are: chrM:576-16024 and the shifted calling is done in chrM:8024-9145.
        // for mouse the original positions are: MT:1-15423, and the shifted calling is done in MT:7423-8299.
        // Note that GRCm38 and GRCm39 use the same MT coordinates, so no need to check the genome version.

        GATK_LEFTALIGNANDTRIMVARIANTS(GATK_FILTERMUECTCALLS.out.vcf_tbi, 'pass-only-final')

        PICARD_MT_COVERAGEATEVERYBASE(mt_alignment.join(shifted_mt_alignment))

        prep_input = GATK_LEFTALIGNANDTRIMVARIANTS.out.vcf_tbi.map { it -> [it[0], it[1], it[2], 'mutect2'] }
                    .mix(MITY_RUNALL.out.vcf_tbi.map { it -> [it[0], it[1], it[2], 'mity'] })
                    .mix(MUTSERVE.out.vcf_tbi.map { it -> [it[0], it[1], it[2], 'mutserve'] })
        PREP_MTDNA_VCF(prep_input)

        BCFTOOLS_MERGECALLERS(PREP_MTDNA_VCF.out.bcf_tbi.groupTuple(size: 3))
        COMPUTE_AVG_AF(BCFTOOLS_MERGECALLERS.out.vcf)

        // Annotation of merged calls
        if (params.gen_org == 'mouse' || params.gen_org == 'human') {
            SNPSIFT_ANNOTATE_DBSNP(COMPUTE_AVG_AF.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
        }

        // If Human
        if (params.gen_org == 'human') {
            SNPSIFT_ANNOTATE_COSMIC(SNPSIFT_ANNOTATE_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
            SNPEFF(SNPSIFT_ANNOTATE_COSMIC.out.vcf, 'MTDNA', 'vcf')
            SNPSIFT_DBNSFP(SNPEFF.out.vcf, 'MTDNA')
            SNPEFF_ONEPERLINE(SNPSIFT_DBNSFP.out.vcf, 'MTDNA')
            SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf, 'mtdna')
        }

        // If Mouse
        if (params.gen_org == 'mouse') {
            SNPEFF(SNPSIFT_ANNOTATE_DBSNP.out.vcf, 'MTDNA', 'vcf')
            SNPEFF_ONEPERLINE(SNPEFF.out.vcf, 'MTDNA')
            SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf, 'mtdna')
        }

}
