#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/wgs_long_reads.nf"
include {param_log} from "${projectDir}/bin/log/wgs_long_reads.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"

include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"

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

include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration"
include {GATK_MERGEVCF;
         GATK_MERGEVCF as GATK_MERGEVCF_UNANNOTATED;
         GATK_MERGEVCF as GATK_MERGEVCF_ANNOTATED} from "${projectDir}/modules/gatk/gatk_mergevcf"

include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_DBSNP} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"

include {PBSV_DISCOVER} from "${projectDir}/modules/pbsv/pbsv_discover"
include {PBSV_CALL} from "${projectDir}/modules/pbsv/pbsv_call"
include {SNIFFLES} from "${projectDir}/modules/sniffles/sniffles"
include {SURVIVOR_MERGE} from "${projectDir}/modules/survivor/survivor_merge"
include {SURVIVOR_VCF_TO_TABLE} from "${projectDir}/modules/survivor/survivor_vcf_to_table"
include {SURVIVOR_SUMMARY} from "${projectDir}/modules/survivor/survivor_summary"
include {SURVIVOR_TO_BED} from "${projectDir}/modules/survivor/survivor_to_bed"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

//help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

if (params.download_data && !params.csv_input) {
    exit 1, "Data download was specified with `--download_data`. However, no input CSV file was specified with `--csv_input`. This is an invalid parameter combination. `--download_data` requires a CSV manifest. See `--help` for information."
}

// prepare reads channel
if (params.csv_input) {

    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))

    if (params.read_type == 'PE'){
        ch_input_sample.map{it -> [it[0], [it[2], it[3]]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    } else if (params.read_type == 'SE') {
        ch_input_sample.map{it -> [it[0], it[2]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

} else if (params.concat_lanes){
  
  if (params.read_type == 'PE'){
    read_ch = Channel
            .fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true, flat:true )
            .map { file, file1, file2 -> tuple(getLibraryId(file), file1, file2) }
            .groupTuple()
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}", checkExists:true, size:1 )
                .map { file, file1 -> tuple(getLibraryId(file), file1) }
                .groupTuple()
                .map{t-> [t[0], t[1].flatten()]}
  }
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}

} else {
  
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}
}


// BEGIN main workflow
workflow WGS_LONG_READS {
// Step 0: Download data and concat Fastq files if needed. 
  if (params.download_data){
      FILE_DOWNLOAD(ch_input_sample)
      FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
      FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
  }

  // Step 00: Concat local Fastq files from CSV input if required.
  if (!params.download_data && params.csv_input){
      CONCATENATE_LOCAL_FILES(ch_input_sample)
      CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
      CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
  }

  // Step 00: Concat local Fastq files if required.
  if (params.concat_lanes && !params.csv_input){
      if (params.read_type == 'PE'){
          CONCATENATE_READS_PE(read_ch)
          read_ch = CONCATENATE_READS_PE.out.concat_fastq
      } else if (params.read_type == 'SE'){
          CONCATENATE_READS_SE(read_ch)
          read_ch = CONCATENATE_READS_SE.out.concat_fastq
      }
  }


  // Read quality and adapter trimming
  FASTP(read_ch)
  
  // FASTQC  
  FASTQC(FASTP.out.trimmed_fastq)

  // START Split FASTQ
  if (params.split_fastq) {
    split_fastq_files = FASTP.out.trimmed_fastq
                         .map{it -> [it[0], it[1]]}
                         .splitFastq(by: params.split_fastq_bin_size, file: true)
                         .map{it -> [it[0], it[1], it[1].name.split('\\.')[-2]]}

                         // from fastp the naming convention will always be *R*.fastq. 
                         // splitFastq adds an increment between *R* and .fastq. 
                         // This can be used to set an 'index' value to make file names unique.
    split_fastq_count = split_fastq_files
                    .groupTuple()
                    .map{sample, reads, index -> [sample, groupKey(sample, index.size())]}
    
    pbmm2_mapping = split_fastq_count
                .combine(split_fastq_files, by:0)
                .map{it -> [it[1], it[2], it[3]] }


    } else {
    
    pbmm2_mapping = FASTP.out.trimmed_fastq

    }
    
    // PBMM2 Alignment
    ch_minimap2_index = file("${params.minimap2_index}")

    // Map reads to indexed genome
    PBMM2_CALL(pbmm2_mapping, ch_minimap2_index)
    
    // 
    if (params.split_fastq) {
      SAMTOOLS_MERGE(PBMM2_CALL.out.pbmm2_bam.map{it -> [it[0], it[1]] }.groupTuple(), 'merged_file')
      bam_file = SAMTOOLS_MERGE.out.bam
    } else {
      bam_file = PBMM2_CALL.out.pbmm2_bam.map{it -> [it[0], it[1]] }
    }
        
    // Begin Merge on Individuals
    if (params.merge_inds) {
      merge_ch = bam_file.join(meta_ch, by: 0)
                         .map{it -> [it[2].ind, it[1]]}
                         .groupTuple()
                         .map{it -> [it[0], it[1], it[1].size()]}
                         .branch{
                                merge: it[2] > 1
                                pass:  it[2] == 1
                          }
      // BAM files are joined to the meta_ch (which contains individual IDs in meta.ind)
      // Individual IDs are taken from the meta_ch, and set to the '0' index of a tuple
      // The '0' index defines the sampleID moving forward, which is `ind` in this case.
      // The size of the group is taken from the 2nd index of the tuple, and used to branch
      // the 'merge' and 'pass' channels.
      // The 'merge' channel is passed to the merge step, and the 'pass' channel is used 
      // to pass the single bam files through to the index step.
    
      merge_input = merge_ch.merge
                          .map{it -> [it[0], it[1]]}

      pass_input = merge_ch.pass
                          .map{it -> [it[0], it[1][0]]}
      // Adjust the tuples for input to merge and index. [sampleID, bam]
      // This removes the 'size' taken above. For single samples, the '0' index of the BAM
      // array is taken, as it was 'grouped' into an array above 
      // and can't be an array going forward. 

      SAMTOOLS_MERGE_IND(merge_input, 'ind_merged_file')
      SAMTOOLS_INDEX_IND(SAMTOOLS_MERGE_IND.out.bam)

      SAMTOOLS_INDEX_SINGLE(pass_input)

      bam_file = SAMTOOLS_MERGE_IND.out.bam
        .mix(pass_input)
      
      index_file  = SAMTOOLS_INDEX_IND.out.bai.mix(SAMTOOLS_INDEX_SINGLE.out.bai)

    } // END merge on individual

    // BAMs are split into two channels, one for SNP/INDEL calling and one for SV calling

    SAMTOOLS_STATS(bam_file)
    MOSDEPTH(bam_file.join(index_file))

    // Make chromsome channel
    chroms = Channel.fromPath("${params.chrom_contigs}")
                    .splitText()
                    .map{it -> it.trim()}
    num_chroms = file(params.chrom_contigs).countLines().toInteger()
    chrom_channel = bam_file.join(index_file).combine(chroms)

    // Find X and Y chromosomes in chroms channel
    haploid_chroms = chroms.filter { it ==~ /(?i).*\b(chr)?X\b.*/ }.map{ it[0] }
            .combine(chroms.filter { it ==~ /(?i).*\b(chr)?Y\b.*/ }.map{ it[0] })
    // Filter the chrom channel to only X and Y. 
    // Because of channel vs. value the filter produces a channel, 
    // which must be manipulated with map to get the value of that channel. 
    if (params.merge_inds) {
      meta_ch = meta_ch
        .map{it -> [it[1].ind, it[1].sex]}
        .unique()
    } else {
      meta_ch = meta_ch
        .map{it -> [it[0], it[1].sex]}
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
    // create select var channels
    select_var_snp = BCFTOOLS_MERGEDEEPVAR_VCF.out.vcf_idx
    select_var_indel = BCFTOOLS_MERGEDEEPVAR_VCF.out.vcf_idx

    // SNP
    GATK_SELECTVARIANTS_SNP(select_var_snp, 'SNP', 'selected_SNP')
    var_filter_snp = GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.idx)
    GATK_VARIANTFILTRATION_SNP(var_filter_snp, 'SNP')

    // INDEL
    GATK_SELECTVARIANTS_INDEL(select_var_indel, 'INDEL', 'selected_INDEL')
    var_filter_indel = GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.idx)
    GATK_VARIANTFILTRATION_INDEL(var_filter_indel, 'INDEL')

    // For other genome, expectation is that dbSNP will not exist.  
    if (params.gen_org=='mouse' | params.gen_org=='human'){
      SNPSIFT_ANNOTATE_SNP_DBSNP(GATK_VARIANTFILTRATION_SNP.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
      SNPSIFT_ANNOTATE_INDEL_DBSNP(GATK_VARIANTFILTRATION_INDEL.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
    }

    // If Human
    if (params.gen_org=='human'){

      // SNP
      SNPSIFT_ANNOTATE_SNP_COSMIC(SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
      SNPEFF_SNP(SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf, 'SNP', 'vcf')
      SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
      SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
      // INDEL
      SNPSIFT_ANNOTATE_INDEL_COSMIC(SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
      SNPEFF_INDEL(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf, 'INDEL', 'vcf')
      SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
      SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')
      
      // Merge SNP and INDEL and Aggregate Stats
      vcf_files_unannotated = SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf.join(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf)
      GATK_MERGEVCF_UNANNOTATED(vcf_files_unannotated, 'SNP_INDEL_filtered_unannotated_final')

      vcf_files_annotated = SNPEFF_ONEPERLINE_SNP.out.vcf.join(SNPEFF_ONEPERLINE_INDEL.out.vcf)
      GATK_MERGEVCF_ANNOTATED(vcf_files_annotated, 'SNP_INDEL_filtered_annotated_final')
      
      SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF_ANNOTATED.out.vcf)
    }

    // If Mouse
    if (params.gen_org=='mouse'){
      // Merge SNP and INDEL

      vcf_files = SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf.join(SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf)

      GATK_MERGEVCF(vcf_files, 'SNP_INDEL_filtered_unannotated_final')

      SNPEFF(GATK_MERGEVCF.out.vcf, 'BOTH', 'vcf')

      SNPEFF_ONEPERLINE(SNPEFF.out.vcf, 'BOTH')

      SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf)
    }


    // If 'Other'
    if (params.gen_org=='other'){
    // For other genomes, there will likely not be SNP EFF annotations, but merge still needs to happen. 
      vcf_files = GATK_VARIANTFILTRATION_SNP.out.vcf.join(GATK_VARIANTFILTRATION_INDEL.out.vcf)

      GATK_MERGEVCF(vcf_files, 'SNP_INDEL_filtered_unannotated_final')
    }

    // SV Calling
    // Call SV with PBSV
    PBSV_DISCOVER(bam_file.join(index_file))
    ch_fasta = params.ref_fa
    PBSV_CALL(PBSV_DISCOVER.out.pbsv_svsig, ch_fasta)

    // Call SV with sniffles
    SNIFFLES(bam_file.join(index_file))

    // Merge caller results

    // Join VCFs together by sampleID and run SURVIVOR merge
    survivor_input = PBSV_CALL.out.pbsv_vcf.join(SNIFFLES.out.sniffles_vcf)
                     .map { it -> tuple(it[0], tuple(it[1], it[2]))}
    SURVIVOR_MERGE(survivor_input)
    SURVIVOR_VCF_TO_TABLE(SURVIVOR_MERGE.out.vcf)
    SURVIVOR_SUMMARY(SURVIVOR_MERGE.out.vcf)

    bed_prep_input = SURVIVOR_VCF_TO_TABLE.out.annotation.join(SURVIVOR_SUMMARY.out.csv)
    SURVIVOR_TO_BED(bed_prep_input)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.idxstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.mosdepth.collect{it[1]}.ifEmpty([]))

    MULTIQC (
      ch_multiqc_files.collect()
    )

}
