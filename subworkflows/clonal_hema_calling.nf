#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {GATK_GETPILEUPSUMMARIES} from "../modules/gatk/gatk_getpileupsummaries_tumorOnly"
include {GATK_CALCULATECONTAMINATION} from "../modules/gatk/gatk_calculatecontamination_tumorOnly"
include {GATK_MUTECT2} from "../modules/gatk/gatk_mutect2_tumorOnly"
include {GATK_LEARNREADORIENTATIONMODEL} from "../modules/gatk/gatk_learnreadorientationmodel"
include {GATK_FILTERMUECTCALLS as GATK_FILTERMUECTCALLS_CONTAM} from "../modules/gatk/gatk_filtermutectcalls_tumoronly_chip_contam"
include {GATK_FILTERMUECTCALLS as GATK_FILTERMUECTCALLS_NOCONTAM} from "../modules/gatk/gatk_filtermutectcalls_tumoronly_chip"
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "../modules/gatk/gatk_selectvariants_chip"
include {BCFTOOLS_STATS} from "../modules/bcftools/bcftools_stats"
include {CHIP_ANNOTATION} from "./chip_annotation"

// Tumor-only Mutect2 chain shared by both species (mm_clonal_hema.nf,
// hs_clonal_hema.nf). No matched normal, no proxy-normal substitution --
// germline filtering is entirely population-resource-based (`gnomad_ref`,
// `pon_ref`), not comparator-based. See notes/germline_resource_design.md and
// the clonal_hema implementation plan (findings 1, 3, 5, 6).
//
// Contamination estimation (GetPileupSummaries/CalculateContamination) only
// runs when `params.contam_ref` is set (human: always; mouse: not yet --
// finding 3). This also decides which FilterMutectCalls variant runs, both
// chip-specific modules: gatk_filtermutectcalls_tumoronly_chip_contam.nf (takes
// the contamination table + segmentation) or gatk_filtermutectcalls_tumoronly_chip.nf
// (stats + orientation-bias model only). The shared WES module is not used: it
// only applies the orientation model when params.ffpe is set. The existing
// `gatk_filtermutectcalls.nf` "plain" module does NOT fit here -- its input
// signature carries PTA's paired-calling meta/normal_name/tumor_name/tool
// fields, which this tumor-only chain never produces.
//
// Orientation-bias modeling (GATK_LEARNREADORIENTATIONMODEL) always runs,
// unconditionally -- unlike WES's `params.ffpe` gate, clonal_hema drops that
// gate entirely (finding 6): at panel depth, rare strand-specific artifacts
// get statistically "confirmed" by sheer read count in a way WES depth rarely
// reaches, and CHIP samples are fresh blood/leukocyte-derived anyway, so the
// FFPE framing doesn't apply.
//
// Exception (2026-09-25): when the input reads are duplex consensus
// (`duplex_consensus` = true, human duplex_seq mode), the learned model is
// NOT passed to FilterMutectCalls. Duplex consensus reads are all F1R2 (see
// the filter modules), so the model would flag nearly every variant; duplex
// consensus with min-reads "2 1 1" already requires both strands. The model
// is still learned in that mode, as a diagnostic output only.

workflow CLONAL_HEMA_CALLING {

    take:
        tumor_name_ch    // [sampleID, bam, bai, tumor_name]
        duplex_consensus // boolean: reads are fgbio duplex consensus (orientation priors off)

    main:
        bam_bai = tumor_name_ch.map{ it -> [it[0], it[1], it[2]] }

        GATK_MUTECT2(tumor_name_ch)
        GATK_LEARNREADORIENTATIONMODEL(GATK_MUTECT2.out.f1r2)

        if (params.contam_ref) {
            GATK_GETPILEUPSUMMARIES(bam_bai)
            GATK_CALCULATECONTAMINATION(GATK_GETPILEUPSUMMARIES.out.pileup_summary)

            filtermutectcalls_input = GATK_MUTECT2.out.vcf_tbi_stats
                .join(GATK_CALCULATECONTAMINATION.out.contam_segments)
                .join(GATK_LEARNREADORIENTATIONMODEL.out.model_file)
                .map{ it -> it + [!duplex_consensus] }
            GATK_FILTERMUECTCALLS_CONTAM(filtermutectcalls_input)

            filtered_vcf = GATK_FILTERMUECTCALLS_CONTAM.out.mutect2_vcf_tbi
            stats = GATK_FILTERMUECTCALLS_CONTAM.out.stats
        } else {
            filtermutectcalls_input = GATK_MUTECT2.out.vcf_tbi_stats
                .join(GATK_LEARNREADORIENTATIONMODEL.out.model_file)
                .map{ it -> it + [!duplex_consensus] }
            GATK_FILTERMUECTCALLS_NOCONTAM(filtermutectcalls_input)

            filtered_vcf = GATK_FILTERMUECTCALLS_NOCONTAM.out.mutect2_vcf_tbi
            stats = GATK_FILTERMUECTCALLS_NOCONTAM.out.stats
        }

        // SNP/INDEL split. `-L ${params.target_gatk}` and
        // `--af-of-alleles-not-in-resource` (left at GATK's dynamic default --
        // finding 5, see the code comment in gatk_mutect2_tumorOnly.nf for the
        // UM-HET3-specific justification, which does not automatically extend
        // to other backgrounds) both live inside GATK_MUTECT2 itself.
        GATK_SELECTVARIANTS_SNP(filtered_vcf, 'SNP', 'selected_SNP')
        GATK_SELECTVARIANTS_INDEL(filtered_vcf, 'INDEL', 'selected_INDEL')

        // Deliberately run on the full filtered VCF (PASS + germline + any
        // other filter reason), not PASS-only -- germline-flagged variants
        // are tagged, not deleted, by FilterMutectCalls, and should stay
        // visible in reporting rather than silently drop out of the variant
        // counts. See notes/milestones.md open item 5.
        BCFTOOLS_STATS(filtered_vcf)

        // Annotation (dbSNP/COSMIC IDs, SnpEff consequence, dbNSFP scores where
        // applicable -- see subworkflows/chip_annotation.nf for what runs per
        // species). Also on the full filtered VCF, same reasoning as above.
        CHIP_ANNOTATION(GATK_SELECTVARIANTS_SNP.out.vcf, GATK_SELECTVARIANTS_INDEL.out.vcf)

    emit:
        filtered_vcf     = filtered_vcf
        stats            = stats
        snp_vcf          = GATK_SELECTVARIANTS_SNP.out.vcf
        indel_vcf        = GATK_SELECTVARIANTS_INDEL.out.vcf
        bcftools_stats   = BCFTOOLS_STATS.out.stats
        annotated_vcf    = CHIP_ANNOTATION.out.annotated_vcf
        annotation_table = CHIP_ANNOTATION.out.annotation_table
}
