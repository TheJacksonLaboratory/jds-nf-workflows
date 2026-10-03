#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CHIP_UMI_CONSENSUS} from "./chip_umi_consensus"
include {GATK_GETSAMPLENAME} from "../modules/gatk/gatk_getsamplename_noMeta"
include {CHIP_SOMATIC_CALLING} from "./chip_somatic_calling"
include {MULTIQC} from "../modules/multiqc/multiqc"
include {CHIP_GENE_REPORT} from "../modules/chip_somatic/chip_gene_report"
include {REORGANIZE_PUBDIR} from "../modules/chip_somatic/reorganize_pubdir"

// Mouse chip_somatic: same UMI/duplex-consensus chain as human
// (subworkflows/chip_umi_consensus.nf), replacing the BWA_MEM ->
// PICARD_SORTSAM -> PICARD_MARKDUPLICATES chain on 2026-09-28 -- the mouse
// pools use the same UMI library prep. Needs the GRCm39 bwa-mem2 index
// (params.ref_fa_bwamem2_indices, mouse override in chip_somatic.config).
// No contamination estimate for mouse (params.contam_ref is null; see
// chip_somatic_calling.nf). No BQSR (finding 2). No SV/CNV/MT/germline-
// HaplotypeCaller/proxy-normal machinery from mm_pta.nf.
workflow MM_CHIP_SOMATIC {

    take:
        concat_ch // [sampleID, meta, [fastq_1, fastq_2]]

    main:
        concat_ch.map{it -> [it[0], it[2]]}.set{read_ch}

        CHIP_UMI_CONSENSUS(read_ch)
        bam_bai = CHIP_UMI_CONSENSUS.out.bam_bai

        GATK_GETSAMPLENAME(bam_bai)
        tumor_name_ch = bam_bai.join(GATK_GETSAMPLENAME.out.sample_name)

        CHIP_SOMATIC_CALLING(tumor_name_ch, params.duplex_seq_bool)

        ch_multiqc_files = CHIP_UMI_CONSENSUS.out.multiqc_files
        ch_multiqc_files = ch_multiqc_files.mix(CHIP_SOMATIC_CALLING.out.stats.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CHIP_SOMATIC_CALLING.out.bcftools_stats.collect{ it -> it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            params.multiqc_config
        )

        // Cross-sample HTML variant/gene report (gene name, coordinate,
        // REF/ALT, VAF, FILTER, dbSNP/COSMIC/SnpEff/dbNSFP annotation --
        // combines every sample's own annotation_table.txt.gz, same DAG
        // position as MULTIQC). See modules/chip_somatic/chip_gene_report.nf.
        CHIP_GENE_REPORT(
            CHIP_SOMATIC_CALLING.out.annotation_table.collect{ it -> it[1]},
            file("${projectDir}/modules/chip_somatic/bin/build_gene_report.py"),
            file("${projectDir}/assets/chip_somatic/jax-gt.svg")
        )

        // Last step: fold two shared modules' (samtools_stats.nf,
        // picard_collecthsmetrics.nf's raw call) non-chip_somatic-convention
        // publish locations back under each sample's own folder -- see
        // bin/chip_somatic/reorganize_pubdir.py. Takes MULTIQC's report as
        // an (unused-content) input purely to run after every per-sample
        // file has already been published.
        sample_ids = read_ch.map{ it -> it[0] }.collect()
        REORGANIZE_PUBDIR(sample_ids, MULTIQC.out.report)
}
