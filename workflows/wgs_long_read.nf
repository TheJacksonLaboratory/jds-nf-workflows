#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/wgs_long_read.nf"
include {param_log} from "${projectDir}/bin/log/wgs_long_read.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"

include {FASTP_LONG} from "${projectDir}/modules/fastp/fastp_long"
include {PBMM2_CALL} from "${projectDir}/modules/pbmm2/pbmm2_call"

include {SAMTOOLS_MERGE;
         SAMTOOLS_MERGE as SAMTOOLS_MERGE_IND} from "${projectDir}/modules/samtools/samtools_merge"
include {SAMTOOLS_STATS} from "${projectDir}/modules/samtools/samtools_stats"
include {MOSDEPTH} from "${projectDir}/modules/mosdepth/mosdepth"

include {DEEPVARIANT} from "${projectDir}/modules/deepvariant/deepvariant"
include {SAMTOOLS_INDEX;
         SAMTOOLS_INDEX as SAMTOOLS_INDEX_IND;
         SAMTOOLS_INDEX as SAMTOOLS_INDEX_SINGLE;} from "${projectDir}/modules/samtools/samtools_index"

include {BCFTOOLS_MERGEDEEPVAR as BCFTOOLS_MERGEDEEPVAR_VCF;
         BCFTOOLS_MERGEDEEPVAR as BCFTOOLS_MERGEDEEPVAR_GVCF} from "${projectDir}/modules/bcftools/bcftools_merge_deepvar_vcfs"

include {GATK_MERGEVCF;
         GATK_MERGEVCF as GATK_MERGEVCF_UNANNOTATED;
         GATK_MERGEVCF as GATK_MERGEVCF_ANNOTATED} from "${projectDir}/modules/gatk/gatk_mergevcf"

include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_DBNSFP} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"

include {PBSV_DISCOVER} from "${projectDir}/modules/pbsv/pbsv_discover"
include {PBSV_CALL} from "${projectDir}/modules/pbsv/pbsv_call"
include {SNIFFLES} from "${projectDir}/modules/sniffles/sniffles"

include {SV_MERGE} from "${projectDir}/modules/r/wgs_sv/wgs_sv_merge"
include {ANNOTATE_SV;
         ANNOTATE_SV as ANNOTATE_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/wgs_sv/annotate_sv"
include {ANNOTATE_GENES_SV;
         ANNOTATE_GENES_SV as ANNOTATE_GENES_SV_SUPPLEMENTAL} from "${projectDir}/modules/r/wgs_sv/annotate_genes_sv"

include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

//help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

if (!params.csv_input) {
    exit 1, "No input CSV file was specified with `--csv_input`. This workflow requires an input CSV file. See `--help` for information."
}

if (params.concat_lanes) {
    exit 1, "Concatenation of lanes was specified with `--concat_lanes`. However, this is not applicable for PacBio long read data. Please remove `--concat_lanes` from your command."
}

// Prepare reads channel
// Note that for PacBio data there is no concept of PE / SE reads. There is only 1 input fastq for each sample. 
// The concept of 'lanes' is also not applicable, as PacBio reads are not split into lanes like Illumina reads. Therefore, lanes are not concatenated.
if (params.csv_input) {
    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
    ch_input_sample.map{it -> [it[0], it[2]]}.set{read_ch}
    ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
} else {
    exit 1, "ERROR: No input CSV file was specified with `--csv_input`. This workflow requires an input CSV file. See `--help` for information."
}

// BEGIN main workflow
workflow wgs_long_read {

    // ADD A DEMUX STEP FOR DATA NOT IN FASTQ FORMAT??

    // FASTP for quality control
    FASTP_LONG(read_ch)

    // PBMM2 Alignment
    ch_minimap2_index = file("${params.minimap2_index}")

    // Map reads to indexed genome
    PBMM2_CALL(FASTP_LONG.out.trimmed_fastq, ch_minimap2_index)

    bam_file = PBMM2_CALL.out.pbmm2_bam.map { it -> [it[0], it[1]] }

    // Begin Merge on Individuals
    if (params.merge_inds) {
        merge_ch = bam_file.join(meta_ch, by: 0)
            .map { it -> [it[2].ind, it[1]] }
            .groupTuple()
            .map { it -> [it[0], it[1], it[1].size()] }
            .branch {
                merge: it[2] > 1
                pass: it[2] == 1
            }
        // BAM files are joined to the meta_ch (which contains individual IDs in meta.ind)
        // Individual IDs are taken from the meta_ch, and set to the '0' index of a tuple
        // The '0' index defines the sampleID moving forward, which is `ind` in this case.
        // The size of the group is taken from the 2nd index of the tuple, and used to branch
        // the 'merge' and 'pass' channels.
        // The 'merge' channel is passed to the merge step, and the 'pass' channel is used 
        // to pass the single bam files through to the index step.

        merge_input = merge_ch.merge
            .map { it -> [it[0], it[1]] }

        pass_input = merge_ch.pass
            .map { it -> [it[0], it[1][0]] }
        // Adjust the tuples for input to merge and index. [sampleID, bam]
        // This removes the 'size' taken above. For single samples, the '0' index of the BAM
        // array is taken, as it was 'grouped' into an array above 
        // and can't be an array going forward. 

        SAMTOOLS_MERGE_IND(merge_input, 'ind_merged_file')
        SAMTOOLS_INDEX_IND(SAMTOOLS_MERGE_IND.out.bam)

        SAMTOOLS_INDEX_SINGLE(pass_input)

        bam_file = SAMTOOLS_MERGE_IND.out.bam
            .mix(pass_input)

        index_file = SAMTOOLS_INDEX_IND.out.bai.mix(SAMTOOLS_INDEX_SINGLE.out.bai)

    } else {
        // If not merging on individuals, just index the BAM files
        SAMTOOLS_INDEX_SINGLE(bam_file)
        index_file = SAMTOOLS_INDEX_SINGLE.out.bai
    } // END merge on individual

    // BAMs are split into two channels, one for SNP/INDEL calling and one for SV calling

    SAMTOOLS_STATS(bam_file)
    MOSDEPTH(bam_file.join(index_file))

    // Make chromsome channel
    chroms = Channel.fromPath("${params.chrom_contigs}")
        .splitText()
        .map { it -> it.trim() }
    num_chroms = file(params.chrom_contigs).countLines().toInteger()
    chrom_channel = bam_file.join(index_file).combine(chroms)

    // Find X and Y chromosomes in chroms channel
    haploid_chroms = chroms.filter { it ==~ /(?i).*\b(chr)?X\b.*/ }.map { it[0] }
        .combine(chroms.filter { it ==~ /(?i).*\b(chr)?Y\b.*/ }.map { it[0] })
    // Filter the chrom channel to only X and Y. 
    // Because of channel vs. value the filter produces a channel, 
    // which must be manipulated with map to get the value of that channel. 
    if (params.merge_inds) {
        meta_ch = meta_ch
            .map { it -> [it[1].ind, it[1].sex] }
            .unique()
    } else {
        meta_ch = meta_ch
            .map { it -> [it[0], it[1].sex] }
            .unique()
    }
    deepvariant_channel = chrom_channel.combine(meta_ch, by: 0).combine(haploid_chroms)

    // Use chrom channel with sex and haploid chrom information in DeepVariant
    DEEPVARIANT(deepvariant_channel)

    // Merge DeepVariant calls
    BCFTOOLS_MERGEDEEPVAR_VCF(DEEPVARIANT.out.vcf_channel.groupTuple(size: num_chroms), 'vcf')
    
    if (params.run_gvcf) {
        BCFTOOLS_MERGEDEEPVAR_GVCF(DEEPVARIANT.out.gvcf_channel.groupTuple(size: num_chroms), 'gvcf')
    }

    // For other genome, expectation is that dbSNP will not exist.  
    if (params.gen_org == 'mouse' || params.gen_org == 'human') {
        SNPSIFT_ANNOTATE_DBSNP(BCFTOOLS_MERGEDEEPVAR_VCF.out.vcf_idx.map { it -> [it[0], it[1]] }, params.dbSNP, params.dbSNP_index, 'dbsnpID')
    }

    // If Human
    if (params.gen_org == 'human') {
        SNPSIFT_ANNOTATE_COSMIC(SNPSIFT_ANNOTATE_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
        SNPEFF(SNPSIFT_ANNOTATE_COSMIC.out.vcf, 'DEEPVAR', 'vcf')
        SNPSIFT_DBNSFP(SNPEFF.out.vcf, 'BOTH')
        SNPEFF_ONEPERLINE(SNPSIFT_DBNSFP.out.vcf, 'BOTH')
        SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf, 'wgs')
    }

    // If Mouse
    if (params.gen_org == 'mouse') {
        SNPEFF(SNPSIFT_ANNOTATE_DBSNP.out.vcf, 'DEEPVAR', 'vcf')
        SNPEFF_ONEPERLINE(SNPEFF.out.vcf, 'BOTH')
        SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf, 'wgs')
    }

    // SV Calling
    // Call SV with PBSV
    PBSV_DISCOVER(bam_file.join(index_file))
    ch_fasta = params.ref_fa
    PBSV_CALL(PBSV_DISCOVER.out.pbsv_svsig, ch_fasta)

    // Call SV with sniffles
    SNIFFLES(bam_file.join(index_file))

    // Merge SV calls
    sv_merge_input = SNIFFLES.out.sniffles_vcf.join(PBSV_CALL.out.pbsv_vcf).map{
        it -> [it[0], [it[1], it[2]]]
    }

    // Get a list of primary chromosomes and exclude chrM (dropRight(1))
    chrom_list = chroms.collect().dropRight(1)
    SV_MERGE(sv_merge_input, chrom_list)

    // Annotate SV calls
    ANNOTATE_SV(SV_MERGE.out.bedpe, "main")
    ANNOTATE_SV_SUPPLEMENTAL(SV_MERGE.out.supp_bedpe, "supplemental")

    ANNOTATE_GENES_SV(ANNOTATE_SV.out.annot_sv_bedpe, "main")
    ANNOTATE_GENES_SV_SUPPLEMENTAL(ANNOTATE_SV_SUPPLEMENTAL.out.annot_sv_bedpe, "supplemental")

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTP_LONG.out.quality_json.collect{ it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.flagstat.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.idxstat.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.mosdepth.collect { it[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect()
    )
}
