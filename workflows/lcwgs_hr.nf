#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/wgs.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"

include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"

include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {SAMTOOLS_MERGE;
         SAMTOOLS_MERGE as SAMTOOLS_MERGE_IND} from "${projectDir}/modules/samtools/samtools_merge"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {SAMTOOLS_STATS} from "${projectDir}/modules/samtools/samtools_stats"
include {create_bamlist} from "${projectDir}/bin/lcwgs/create_bamlist.nf"
include {MOSDEPTH} from "${projectDir}/modules/mosdepth/mosdepth"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

// help if needed
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


workflow LCWGS_HR{
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
  trimmer_input = read_ch
  FASTP(trimmer_input)
    
  FASTQC(FASTP.out.trimmed_fastq)

  READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")
  bwa_mem_mapping = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
                                           .map{it -> [it[0], it[1], 'aln', it[2]]}
  
  // BWA-MEM Alignment
  BWA_MEM(bwa_mem_mapping)
  PICARD_SORTSAM(BWA_MEM.out.sam, 'coordinate')
  
  // Mark duplicates
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  bam_file = PICARD_MARKDUPLICATES.out.dedup_bam
  index_file = PICARD_MARKDUPLICATES.out.dedup_bai

  // QC stats
  SAMTOOLS_STATS(bam_file)
  MOSDEPTH(bam_file.join(index_file))

  // Prepare for QUILT
  bams = bam_file.map {tuple -> tuple[1]}
                 .collect()
                 .map{bamlist -> [bamlist, "no_downsample"]}
  bamlist_ch = create_bamlist(bams)
  bamlist_ch.view()

    

  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.flagstat.collect { it[1] }.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.idxstat.collect { it[1] }.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { it[1] }.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.mosdepth.collect { it[1] }.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )

}

