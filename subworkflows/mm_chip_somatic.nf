#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {FASTP} from "../modules/fastp/fastp"
include {FASTQC} from "../modules/fastqc/fastqc"
include {READ_GROUPS} from "../modules/utility_modules/read_groups"
include {BWA_MEM} from "../modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "../modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "../modules/picard/picard_markduplicates"
include {PICARD_COLLECTHSMETRICS} from "../modules/picard/picard_collecthsmetrics"
include {GATK_GETSAMPLENAME} from "../modules/gatk/gatk_getsamplename_noMeta"
include {CHIP_SOMATIC_CALLING} from "./chip_somatic_calling"
include {MULTIQC} from "../modules/multiqc/multiqc"

// Mouse alignment for chip_somatic. Deliberately plain -- BWA_MEM ->
// PICARD_SORTSAM(coordinate) -> PICARD_MARKDUPLICATES, matching mm_pta.nf's
// existing core chain, no BQSR (finding 2: not a mouse-resource gap, panel
// depth doesn't supply BQSR's needed covariate diversity for either species,
// plus a real risk BQSR suppresses true recurrent low-VAF signal at exactly
// the hotspot positions this panel targets). No short-alignment-marking
// either (finding 4) -- that chain exists in this repo only to protect PTA's
// SV/CNV caller stack, which this workflow doesn't run.
//
// No SV/CNV/MT/germline-HaplotypeCaller/proxy-normal machinery from
// mm_pta.nf -- none of it applies to tumor-only panel calling.
//
// Dedup: PICARD_MARKDUPLICATES assumes hybrid-capture (random-fragment)
// library prep. This is a known open item (see the implementation plan) --
// if the CHIP panel turns out to be amplicon-based or UMI-tagged, this step
// needs to change to primer-clipping or UMI-consensus dedup instead. Not
// resolved by this file.

workflow MM_CHIP_SOMATIC {

    take:
        concat_ch // [sampleID, meta, [fastq_1, fastq_2]]

    main:
        concat_ch.map{it -> [it[0], it[2]]}.set{read_ch}

        FASTP(read_ch)
        FASTQC(FASTP.out.trimmed_fastq)
        READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")

        bwa_mem_mapping = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
                          .map{it -> [it[0], it[1], 'aln', it[2]]}
        BWA_MEM(bwa_mem_mapping)
        PICARD_SORTSAM(BWA_MEM.out.sam, 'coordinate')
        PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

        bam_bai = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
        PICARD_COLLECTHSMETRICS(bam_bai)

        GATK_GETSAMPLENAME(bam_bai)
        tumor_name_ch = bam_bai.join(GATK_GETSAMPLENAME.out.sample_name)

        // No UMI consensus in the mouse path (PICARD_MARKDUPLICATES), so reads
        // keep both pair orientations and the orientation model applies.
        CHIP_SOMATIC_CALLING(tumor_name_ch, false)

        ch_multiqc_files = channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTHSMETRICS.out.hsmetrics.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CHIP_SOMATIC_CALLING.out.stats.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CHIP_SOMATIC_CALLING.out.bcftools_stats.collect{ it -> it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            params.multiqc_config
        )
}
