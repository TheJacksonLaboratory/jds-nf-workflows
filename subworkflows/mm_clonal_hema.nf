#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CHIP_UMI_CONSENSUS} from "./chip_umi_consensus"
include {GATK_GETSAMPLENAME} from "../modules/gatk/gatk_getsamplename_noMeta"
include {CLONAL_HEMA_CALLING} from "./clonal_hema_calling"
include {MULTIQC} from "../modules/multiqc/multiqc_clonal_hema"
include {CHIP_GENE_REPORT} from "../modules/clonal_hema/chip_gene_report"
include {REORGANIZE_PUBDIR} from "../modules/clonal_hema/reorganize_pubdir"

// Mouse clonal_hema: same UMI/duplex-consensus chain as human
// (subworkflows/chip_umi_consensus.nf), replacing the BWA_MEM ->
// PICARD_SORTSAM -> PICARD_MARKDUPLICATES chain on 2026-09-28 -- the mouse
// pools use the same UMI library prep. Needs the GRCm39 bwa-mem2 index
// (params.ref_fa_bwamem2_indices, mouse override in clonal_hema.config).
// No contamination estimate for mouse (params.contam_ref is null; see
// clonal_hema_calling.nf). No BQSR (finding 2). No SV/CNV/MT/germline-
// HaplotypeCaller/proxy-normal machinery from mm_pta.nf.
workflow MM_CLONAL_HEMA {

    take:
        concat_ch // [sampleID, meta, [fastq_1, fastq_2]]

    main:
        concat_ch.map{it -> [it[0], it[2]]}.set{read_ch}

        CHIP_UMI_CONSENSUS(read_ch)
        bam_bai = CHIP_UMI_CONSENSUS.out.bam_bai

        GATK_GETSAMPLENAME(bam_bai)
        tumor_name_ch = bam_bai.join(GATK_GETSAMPLENAME.out.sample_name)

        CLONAL_HEMA_CALLING(tumor_name_ch, params.duplex_seq_bool)

        ch_multiqc_files = CHIP_UMI_CONSENSUS.out.multiqc_files
        ch_multiqc_files = ch_multiqc_files.mix(CLONAL_HEMA_CALLING.out.stats.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CLONAL_HEMA_CALLING.out.bcftools_stats.collect{ it -> it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            params.multiqc_config
        )

        // Cross-sample HTML variant/gene report (gene name, coordinate,
        // REF/ALT, VAF, FILTER, dbSNP/COSMIC/SnpEff/dbNSFP annotation --
        // combines every sample's own annotation_table.txt.gz, same DAG
        // position as MULTIQC). See modules/clonal_hema/chip_gene_report.nf.
        CHIP_GENE_REPORT(
            CLONAL_HEMA_CALLING.out.annotation_table.collect{ it -> it[1]},
            file("${projectDir}/modules/clonal_hema/bin/build_gene_report.py"),
            file("${projectDir}/assets/clonal_hema/jax-gt.svg"),
            'mouse',
            // the default mouse (B6) germline resource is an empty placeholder: no germline filtering
            params.gnomad_ref.toString().contains('empty_germline_resource') ? 'no' : 'yes',
            CHIP_UMI_CONSENSUS.out.gene_coverage,
            CHIP_UMI_CONSENSUS.out.on_target
        )

        // Last step: fold two shared modules' (samtools_stats.nf,
        // picard_collecthsmetrics.nf's raw call) non-clonal_hema-convention
        // publish locations back under each sample's own folder -- see
        // bin/clonal_hema/reorganize_pubdir.py. Takes MULTIQC's report as
        // an (unused-content) input purely to run after every per-sample
        // file has already been published.
        sample_ids = read_ch.map{ it -> it[0] }.collect()
        REORGANIZE_PUBDIR(sample_ids, MULTIQC.out.report)
}
