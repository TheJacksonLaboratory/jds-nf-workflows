#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {SAMTOOLS_FAIDX} from "${projectDir}/modules/samtools/samtools_faidx"

include {SMOOVE_CALL} from "${projectDir}/modules/smoove/smoove_call_germline"
include {MANTA_CALL} from "${projectDir}/modules/illumina/manta_germline"
include {SVABA} from "${projectDir}/modules/svaba/svaba_germline"
include {GATK_UPDATEVCFSEQUENCEDICTIONARY as SVABA_SV_UPDATE_DICTIONARY} from "${projectDir}/modules/gatk/gatk_updatevcfsequencedictionary_germline"
include {DELLY_CALL_GERMLINE} from "${projectDir}/modules/delly/delly_call_germline"

include {DELLY_CNV_GERMLINE} from "${projectDir}/modules/delly/delly_cnv_germline"
include {BCFTOOLS_QUERY_DELLY_CNV} from "${projectDir}/modules/bcftools/bcftools_query_delly_cnv_germline"

include {BCFTOOLS_REHEAD_SORT as REHEAD_SORT_LUMPY;
         BCFTOOLS_REHEAD_SORT as REHEAD_SORT_DELLY;
         BCFTOOLS_REHEAD_SORT as REHEAD_SORT_MANTA;
         BCFTOOLS_REHEAD_SORT as REHEAD_SORT_SVABA} from "${projectDir}/modules/bcftools/bcftools_rehead_sort"

include {DUPHOLD as DUPHOLD_DELLY;
         DUPHOLD as DUPHOLD_LUMPY;
         DUPHOLD as DUPHOLD_MANTA;
         DUPHOLD as DUPHOLD_SVABA} from "${projectDir}/modules/duphold/duphold"

include {BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_DELLY;
         BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_LUMPY;
         BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_MANTA;
         BCFTOOLS_DUPHOLD_FILTER as BCFTOOLS_DUPHOLD_FILTER_SVABA} from "${projectDir}/modules/bcftools/bcftools_duphold_filter"

include {ANNOTATE_DELLY_CNV} from "${projectDir}/modules/r/wgs_sv/annotate_delly_cnv"

include {SV_MERGE} from "${projectDir}/modules/r/wgs_sv/wgs_sv_merge"
include {ANNOTATE_SV;
         ANNOTATE_SV as ANNOTATE_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/wgs_sv/annotate_sv"
include {ANNOTATE_GENES_SV;
         ANNOTATE_GENES_SV as ANNOTATE_GENES_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/wgs_sv/annotate_genes_sv"
include {ANNOTATE_SV_WITH_CNV;
         ANNOTATE_SV_WITH_CNV as ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL} from "${projectDir}/modules/r/wgs_sv/annotate_sv_with_cnv"
include {FILTER_BEDPE;
         FILTER_BEDPE as FILTER_BEDPE_SUPPLEMENTAL} from "${projectDir}/modules/r/wgs_sv/filter_bedpe"

workflow WGS_SV {

    take:
        bam_bai_ch
        // [sampleID, bam, bai]
    
    main:
        // Read a list of chromosome names from a parameter. These are provided to several tools. 
        chroms = Channel
            .fromPath("${params.chrom_contigs}")
            .splitText()
            .map{it -> it.trim()}

        // Get a list of primary chromosomes and exclude chrM (dropRight(1))
        chrom_list = chroms.collect().dropRight(1)

        // Index reference fasta
        faidx_input = ['ref_fasta', params.ref_fa]
        SAMTOOLS_FAIDX(faidx_input)
        fasta_index = SAMTOOLS_FAIDX.out.fai.map{it -> [params.ref_fa, it[1]]}

        // Call SV with SMOOVE
        SMOOVE_CALL(bam_bai_ch)
        REHEAD_SORT_LUMPY(SMOOVE_CALL.out.lumpy_vcf, "lumpy", fasta_index)

        // * Manta

        // Call SV with Manta
        MANTA_CALL(bam_bai_ch, fasta_index)
        REHEAD_SORT_MANTA(MANTA_CALL.out.manta_sv, "manta", fasta_index)

        // * SVABA

        SVABA(bam_bai_ch)
        SVABA_SV_UPDATE_DICTIONARY(SVABA.out.svaba_sv_vcf_tbi, 'svaba') 
        // Note: SVABA jumbles the dict header in the VCF, and must be adjusted prior to use in GATK tools.
        REHEAD_SORT_SVABA(SVABA_SV_UPDATE_DICTIONARY.out.vcf_tbi.map{it -> [it[0], it[1]]}, "svaba", fasta_index)

        // * Delly

        // Call SV with Delly
        DELLY_CALL_GERMLINE(bam_bai_ch, fasta_index)
        REHEAD_SORT_DELLY(DELLY_CALL_GERMLINE.out.delly_bcf, "delly_sv", fasta_index)
    
        // Call CNV with Delly
        DELLY_CNV_GERMLINE(bam_bai_ch, fasta_index)
        // Note: For multi-sample germline calling: DELLY_CLASSIFY should be used on merged BCF. See: https://github.com/dellytools/delly?tab=readme-ov-file#germline-cnv-calling
        //       Due to single sample nature of this workflow, calls are not passed through the classifier, and reported as is.  
        BCFTOOLS_QUERY_DELLY_CNV(DELLY_CNV_GERMLINE.out.delly_bcf.join(DELLY_CNV_GERMLINE.out.delly_csi))
        // PLOT_DELLY_CNV could be included here. For this, the bin/pta/delly_cnv_plot.r script needs to be modified for germline.
        ANNOTATE_DELLY_CNV(BCFTOOLS_QUERY_DELLY_CNV.out.segmentation_file, chrom_list)

        // Duphold
        DUPHOLD_DELLY(bam_bai_ch.join(REHEAD_SORT_DELLY.out.vcf_sort), fasta_index, 'delly_sv') 
        DUPHOLD_MANTA(bam_bai_ch.join(REHEAD_SORT_MANTA.out.vcf_sort), fasta_index, 'manta')
        DUPHOLD_SVABA(bam_bai_ch.join(REHEAD_SORT_SVABA.out.vcf_sort), fasta_index, 'svaba')

        // Duphold Filter
        BCFTOOLS_DUPHOLD_FILTER_DELLY(DUPHOLD_DELLY.out.vcf, 'delly_sv')
        BCFTOOLS_DUPHOLD_FILTER_MANTA(DUPHOLD_MANTA.out.vcf, 'manta')
        BCFTOOLS_DUPHOLD_FILTER_LUMPY(REHEAD_SORT_LUMPY.out.vcf_sort, 'lumpy')
        BCFTOOLS_DUPHOLD_FILTER_SVABA(DUPHOLD_SVABA.out.vcf, 'svaba')


        // * Merge callers and annotate results

        // Join VCFs together by sampleID and merge via NYGC based script
        merge_input = BCFTOOLS_DUPHOLD_FILTER_MANTA.out.vcf.join(BCFTOOLS_DUPHOLD_FILTER_DELLY.out.vcf).join(BCFTOOLS_DUPHOLD_FILTER_LUMPY.out.vcf).join(BCFTOOLS_DUPHOLD_FILTER_SVABA.out.vcf)
                        .map { it -> tuple(it[0], tuple(it[1], it[2], it[3], it[4])) }
        SV_MERGE(merge_input, chrom_list)
        ANNOTATE_SV(SV_MERGE.out.bedpe, "main")
        ANNOTATE_SV_SUPPLEMENTAL(SV_MERGE.out.supp_bedpe, "supplemental")

        ANNOTATE_GENES_SV(ANNOTATE_SV.out.annot_sv_bedpe, "main")
        ANNOTATE_GENES_SV_SUPPLEMENTAL(ANNOTATE_SV_SUPPLEMENTAL.out.annot_sv_bedpe, "supplemental")
        
        annot_sv_cnv_input = ANNOTATE_DELLY_CNV.out.delly_annot.join(ANNOTATE_GENES_SV.out.annot_sv_genes_bedpe)
        ANNOTATE_SV_WITH_CNV(annot_sv_cnv_input, "main")
        FILTER_BEDPE(ANNOTATE_SV_WITH_CNV.out.sv_genes_cnv_bedpe, "main")

        annot_sv_cnv_suppl_input = ANNOTATE_DELLY_CNV.out.delly_annot.join(ANNOTATE_GENES_SV_SUPPLEMENTAL.out.annot_sv_genes_bedpe)
        ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL(annot_sv_cnv_suppl_input, "supplemental")
        FILTER_BEDPE_SUPPLEMENTAL(ANNOTATE_SV_WITH_CNV_SUPPLEMENTAL.out.sv_genes_cnv_bedpe, "supplemental")   
}
