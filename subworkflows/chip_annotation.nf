#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_COSMIC} from "../modules/snpeff_snpsift/snpsift_annotate_chip"
include {BCFTOOLS_ANNOTATE_CLINVAR as BCFTOOLS_ANNOTATE_CLINVAR_SNP;
         BCFTOOLS_ANNOTATE_CLINVAR as BCFTOOLS_ANNOTATE_CLINVAR_INDEL} from "../modules/bcftools/bcftools_annotate_clinvar_chip"
include {SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "../modules/snpeff_snpsift/snpeff_chip"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "../modules/snpeff_snpsift/snpsift_dbnsfp_chip"
include {SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "../modules/snpeff_snpsift/snpeff_oneperline_chip"
include {GATK_MERGEVCF} from "../modules/gatk/gatk_mergevcf_chip"
include {SNPSIFT_EXTRACTFIELDS} from "../modules/snpeff_snpsift/snpsift_extractfields_chip"

// Annotation for clonal_hema, both species, added 2026-09-29 -- learned from
// --workflow somatic_wes's own annotation chain (dbSNP + COSMIC ID annotation,
// SnpEff consequence calls, dbNSFP in-silico scores), reusing the same
// resources somatic_wes already uses (verified present in omics_share) via new
// chip-specific modules. WES's own modules/subworkflow are untouched.
//
// Deliberately NOT copied from WES: the GATK_VARIANTFILTRATION hard-filter step
// (DP < 25, FS > 60) that runs between SelectVariants and annotation in WES.
// At clonal_hema's panel depths (hundreds to thousands of x even at low VAF),
// DP < 25 essentially never fires, and FS > 60 strand-bias filtering duplicates
// what FilterMutectCalls' own orientation/strand_bias filters already do
// (verified against real FILTER values in our synthetic tests). Adding it would
// risk a second, WES-tuned filter layer overwriting FILTER reasons we've
// already characterized. So this subworkflow runs directly on
// CLONAL_HEMA_CALLING's SelectVariants output.
//
// Human: dbSNP -> COSMIC -> ClinVar -> SnpEff -> dbNSFP -> oneperline, both SNP and INDEL.
// Mouse: dbSNP/COSMIC/dbNSFP are human-only resources (see the comment on
// params.dbSNP in config/clonal_hema.config for why the one mouse dbSNP-like
// resource that exists was NOT wired in -- GRCm38 coordinates, wrong build).
// Mouse gets SnpEff consequence/gene annotation only (a real GRCm39 database
// exists), gated on params.dbSNP/params.cosmic/params.dbNSFP being non-null
// rather than on params.gen_org directly, so this subworkflow doesn't need its
// own species branching -- it just follows what's configured.

workflow CHIP_ANNOTATION {

    take:
        snp_vcf   // [sampleID, selected_SNP.vcf.gz]
        indel_vcf // [sampleID, selected_INDEL.vcf.gz]

    main:
        snp_ch = snp_vcf
        indel_ch = indel_vcf

        if (params.dbSNP) {
            SNPSIFT_ANNOTATE_SNP_DBSNP(snp_ch, params.dbSNP, params.dbSNP_index, 'dbsnpID')
            SNPSIFT_ANNOTATE_INDEL_DBSNP(indel_ch, params.dbSNP, params.dbSNP_index, 'dbsnpID')
            snp_ch = SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf
            indel_ch = SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf
        }

        if (params.cosmic) {
            SNPSIFT_ANNOTATE_SNP_COSMIC(snp_ch, params.cosmic, params.cosmic_index, 'cosmicID')
            SNPSIFT_ANNOTATE_INDEL_COSMIC(indel_ch, params.cosmic, params.cosmic_index, 'cosmicID')
            snp_ch = SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf
            indel_ch = SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf
        }

        // ClinVar (human only: params.clinvar is null for mouse). INFO fields CLNSIG etc.
        if (params.clinvar) {
            BCFTOOLS_ANNOTATE_CLINVAR_SNP(snp_ch, params.clinvar, params.clinvar_index, 'clinvar')
            BCFTOOLS_ANNOTATE_CLINVAR_INDEL(indel_ch, params.clinvar, params.clinvar_index, 'clinvar')
            snp_ch = BCFTOOLS_ANNOTATE_CLINVAR_SNP.out.vcf
            indel_ch = BCFTOOLS_ANNOTATE_CLINVAR_INDEL.out.vcf
        }

        SNPEFF_SNP(snp_ch, 'SNP')
        SNPEFF_INDEL(indel_ch, 'INDEL')
        snp_ch = SNPEFF_SNP.out.vcf
        indel_ch = SNPEFF_INDEL.out.vcf

        if (params.dbNSFP) {
            SNPSIFT_DBNSFP_SNP(snp_ch, 'SNP')
            SNPSIFT_DBNSFP_INDEL(indel_ch, 'INDEL')
            snp_ch = SNPSIFT_DBNSFP_SNP.out.vcf
            indel_ch = SNPSIFT_DBNSFP_INDEL.out.vcf
        }

        SNPEFF_ONEPERLINE_SNP(snp_ch, 'SNP')
        SNPEFF_ONEPERLINE_INDEL(indel_ch, 'INDEL')

        merged = SNPEFF_ONEPERLINE_SNP.out.vcf.join(SNPEFF_ONEPERLINE_INDEL.out.vcf)
        GATK_MERGEVCF(merged, 'annotated_final')

        SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF.out.vcf, params.gen_org)

    emit:
        annotated_vcf   = GATK_MERGEVCF.out.vcf
        annotated_idx   = GATK_MERGEVCF.out.idx
        annotation_table = SNPSIFT_EXTRACTFIELDS.out.txt
}
