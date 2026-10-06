#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {FASTQC} from "../modules/fastqc/fastqc"
include {FGBIO_FASTQTOBAM} from "../modules/fgbio/fgbio_fastqtobam"
include {FGBIO_ALIGN_BAM as FGBIO_ALIGN_RAW_BAM;
         FGBIO_ALIGN_BAM as FGBIO_ALIGN_CONS} from "../modules/fgbio/fgbio_align_bam"
include {FGBIO_GROUPREADSBYUMI} from "../modules/fgbio/fgbio_groupreadsbyumi"
include {FGBIO_CALLANDFILTERDUPLEXCONSENSUSREADS} from "../modules/fgbio/fgbio_callandfilterduplexconsensusreads"
include {FGBIO_CALLANDFILTERMOLECULARCONSENSUSREADS} from "../modules/fgbio/fgbio_callandfiltermolecularconsensusreads"
include {FGBIO_COLLECTDUPLEXSEQMETRICS} from "../modules/fgbio/fgbio_collectduplexseqmetrics"
include {DUPLEX_METRICS_MQC} from "../modules/fgbio/duplex_metrics_mqc"
include {SAMTOOLS_STATS as SAMTOOLS_STATS_PREDEDUP;
         SAMTOOLS_STATS as SAMTOOLS_STATS_POSTDEDUP} from "../modules/samtools/samtools_stats"
include {SAMTOOLS_SORT} from "../modules/samtools/samtools_sort"
include {SAMTOOLS_INDEX} from "../modules/samtools/samtools_index"
include {PICARD_COLLECTHSMETRICS;
         PICARD_COLLECTHSMETRICS as PICARD_COLLECTHSMETRICS_RAW} from "../modules/picard/picard_collecthsmetrics_chip"

// Shared UMI/duplex-consensus chain for clonal_hema, both species (Twist
// Bioscience UMI adapter system -- DOC-001337/DOC-001506), replacing the
// plain BWA_MEM -> PICARD_SORTSAM -> PICARD_MARKDUPLICATES chain used until
// 2026-09-14. See notes/milestones.md open item 6 for the full rationale:
// the panel's real library prep uses inline dual UMIs (5bp + 2bp skip, both
// mates), so extraction is a correctness requirement, not an optional
// feature -- feeding UMI+skip bases into a normal aligner as if they were
// genomic sequence would corrupt every alignment. Moved here from
// hs_clonal_hema.nf on 2026-09-28 so mouse uses the same chain (needs a
// bwa-mem2 index per species: params.ref_fa_bwamem2_indices).
//
// FASTP is deliberately NOT run before UMI extraction -- Twist's UMI sits
// at a fixed position at the very start of each read (read structure
// 5M2S+T), and any 5' trimming/adapter-detection could shift or damage
// those exact bases before fgbio ever sees them. FASTQC still runs on the
// raw reads for QC visibility only (matches the reference pipeline this
// chain was ported from, notes/milestones.md open item 6).
//
// No BQSR (finding 2, panel depth doesn't supply the covariate diversity
// BQSR needs, plus real risk of suppressing true recurrent low-VAF signal).
//
// Dropped entirely, not gated: HLA typing/OptiType, MSIsensor2, PDX/Xengsort,
// the ancestry SNP panel, and the COSMIC/dbNSFP/SnpEff annotation cascade --
// none of these fit a targeted CHIP panel. See the implementation plan for
// the full reasoning on each.

workflow CHIP_UMI_CONSENSUS {

    take:
        read_ch // [sampleID, [fastq_1, fastq_2]]

    main:
        FASTQC(read_ch)

        FGBIO_FASTQTOBAM(read_ch)
        FGBIO_ALIGN_RAW_BAM(FGBIO_FASTQTOBAM.out.bam, "template-coordinate")
        FGBIO_GROUPREADSBYUMI(FGBIO_ALIGN_RAW_BAM.out.bam)

        // Duplex (both-strand-confirmed error correction) vs. molecular
        // (single-strand consensus, no duplex requirement -- e.g. min_reads=1
        // to retain singletons) consensus calling, per params.duplex_seq --
        // both branches ported from the same reference pipeline (5baseTAPS's
        // own params.duplex_seq toggle). See config/clonal_hema.config for
        // the call_min_reads format difference between the two modes.
        if (params.duplex_seq_bool) {
            FGBIO_CALLANDFILTERDUPLEXCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam)
            ch_consensus_bam = FGBIO_CALLANDFILTERDUPLEXCONSENSUSREADS.out.bam
        } else {
            FGBIO_CALLANDFILTERMOLECULARCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam)
            ch_consensus_bam = FGBIO_CALLANDFILTERMOLECULARCONSENSUSREADS.out.bam
        }
        FGBIO_ALIGN_CONS(ch_consensus_bam, "coordinate")

        // QC only, not joined back into the main chain: total raw reads
        // (pre-dedup) vs. total reads after UMI/duplex-consensus collapsing
        // (post-dedup), replacing the read-count visibility PICARD_MARKDUPLICATES'
        // dedup_metrics used to give. `_raw` suffix avoids a filename collision
        // with the post-dedup report -- SAMTOOLS_STATS names output files after
        // its sampleID input directly, and both calls would otherwise write
        // "{sampleID}.flagstat" to the same published directory. Restored to
        // the real sampleID immediately after, for joining below.
        SAMTOOLS_STATS_POSTDEDUP(FGBIO_ALIGN_CONS.out.bam)

        // Real raw (pre-consensus) on-target depth, as designed -- FGBIO_ALIGN_RAW_BAM
        // is template-coordinate sorted (required by GroupReadsByUmi), not usable
        // directly by HsMetrics, so it's re-sorted to standard coordinate order here.
        // Suffixed like SAMTOOLS_STATS_PREDEDUP above to avoid a filename collision
        // with the consensus HsMetrics output in the same published directory; kept
        // out of MultiQC's own file feed for the same reason (folded into
        // DUPLEX_METRICS_MQC's "Raw depth" columns instead, see below).
        raw_bam_suffixed = FGBIO_ALIGN_RAW_BAM.out.bam.map{ it -> ["${it[0]}_raw", it[1]] }
        SAMTOOLS_SORT(raw_bam_suffixed, '', 'bam')
        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_file)
        // Pre-dedup stats run on the coordinate-sorted copy, not the
        // template-coordinate FGBIO_ALIGN_RAW_BAM output: `samtools idxstats`
        // without an index requires position order, and template-coordinate
        // order breaks it at any discordant pair (mates on different
        // chromosomes), which real data always contains.
        SAMTOOLS_STATS_PREDEDUP(SAMTOOLS_SORT.out.sorted_file)
        prededup_flagstat = SAMTOOLS_STATS_PREDEDUP.out.flagstat.map{ it -> [it[0].replaceAll('_raw$', ''), it[1]] }
        raw_bam_bai = SAMTOOLS_SORT.out.sorted_file.join(SAMTOOLS_INDEX.out.bai)
        PICARD_COLLECTHSMETRICS_RAW(raw_bam_bai)
        raw_hsmetrics = PICARD_COLLECTHSMETRICS_RAW.out.hsmetrics.map{ it -> [it[0].replaceAll('_raw$', ''), it[1]] }

        // Duplex-yield/fold-reduction MultiQC summary (bin/shared/duplex_metrics.py,
        // ported from the same reference pipeline as the fgbio modules above)
        // -- duplex-specific, gated the same way 5baseTAPS gates its own
        // DUPLEX_QC call: CollectDuplexSeqMetrics expects the A/B duplex MI
        // structure that only the Paired GroupReadsByUmi strategy produces,
        // so it's meaningless (and not run) in molecular mode. Molecular
        // mode currently has no equivalent raw-depth/fold-reduction table --
        // a follow-up, not built here.
        ch_duplex_mqc = channel.empty()
        if (params.duplex_seq_bool) {
            FGBIO_COLLECTDUPLEXSEQMETRICS(FGBIO_GROUPREADSBYUMI.out.bam)
            duplex_mqc_input = prededup_flagstat
                .join(SAMTOOLS_STATS_POSTDEDUP.out.flagstat)
                .join(raw_hsmetrics)
                .join(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics)
            DUPLEX_METRICS_MQC(duplex_mqc_input)
            ch_duplex_mqc = DUPLEX_METRICS_MQC.out.mqc.collect{ it -> it[1]}.ifEmpty([]).flatten()
        }

        bam_bai = FGBIO_ALIGN_CONS.out.bam.join(FGBIO_ALIGN_CONS.out.bai)
        PICARD_COLLECTHSMETRICS(bam_bai)

        ch_multiqc_files = channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{ it -> it[1]}.ifEmpty([]))
        // SAMTOOLS_STATS_PREDEDUP not fed to MultiQC -- its own file is
        // suffixed "_raw" (to avoid a filename collision with the
        // post-dedup report) which MultiQC would otherwise render as a
        // second row for the same sample. DUPLEX_METRICS_MQC's own table
        // already reports raw vs. consensus reads as columns in one row.
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_POSTDEDUP.out.flagstat.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.collect{ it -> it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_duplex_mqc)
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTHSMETRICS.out.hsmetrics.collect{ it -> it[1]}.ifEmpty([]))

    emit:
        bam_bai       = bam_bai          // [sampleID, consensus bam, bai]
        multiqc_files = ch_multiqc_files // QC files for MULTIQC
}
