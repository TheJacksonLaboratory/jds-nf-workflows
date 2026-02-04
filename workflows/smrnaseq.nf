#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Adapted from: nf-core/smrnaseq 2.2.4 and 2.4.0 workflows

// import modules
include {help} from "${projectDir}/bin/help/smrnaseq.nf"
include {param_log} from "${projectDir}/bin/log/smrnaseq.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {INPUT_CHECK} from "${projectDir}/subworkflows/input_check"
include {FASTQC as FASTQC_RAW;
         FASTQC as FASTQC_TRIM} from "${projectDir}/modules/fastqc/fastqc"
include {FASTP} from "${projectDir}/modules/fastp/fastp_smrna"

include {SEQCLUSTER_SEQUENCES} from "${projectDir}/modules/seqcluster/seqcluster_collapse"

include {MIRTRACE_RUN} from "${projectDir}/modules/mirtrace/mirtrace"

include { BOWTIE_MAP_CONTAMINANTS as MAP_TRNA
          BOWTIE_MAP_CONTAMINANTS as MAP_CDNA
          BOWTIE_MAP_CONTAMINANTS as MAP_NCRNA
          BOWTIE_MAP_CONTAMINANTS as MAP_OTHER } from "${projectDir}/modules/bowtie2/bowtie2_map_contaminants"

include {BOWTIE_MAP_SEQ  as BOWTIE_MAP_MATURE;
         BOWTIE_MAP_SEQ  as BOWTIE_MAP_HAIRPIN;
         BOWTIE_MAP_SEQ  as BOWTIE_MAP_GENOME;
         BOWTIE_MAP_SEQ  as BOWTIE_MAP_SEQCLUSTER } from "${projectDir}/modules/bowtie/bowtie_map_mirna"

include {SAMTOOLS_SORT as SAMTOOLS_SORT_MATURE;
         SAMTOOLS_SORT as SAMTOOLS_SORT_HAIRPIN;
         SAMTOOLS_SORT as SAMTOOLS_SORT_GENOME;
         SAMTOOLS_SORT as SAMTOOLS_SORT_SEQCLUSTER } from "${projectDir}/modules/samtools/samtools_sort"

include {SAMTOOLS_INDEX as SAMTOOLS_INDEX_MATURE;
         SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAIRPIN;
         SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENOME;
         SAMTOOLS_INDEX as SAMTOOLS_INDEX_SEQCLUSTER } from "${projectDir}/modules/samtools/samtools_index"

include {SAMTOOLS_STATS as SAMTOOLS_STATS_MATURE;
         SAMTOOLS_STATS as SAMTOOLS_STATS_HAIRPIN;
         SAMTOOLS_STATS as SAMTOOLS_STATS_GENOME;
         SAMTOOLS_STATS as SAMTOOLS_STATS_SEQCLUSTER } from "${projectDir}/modules/samtools/samtools_stats"

include { MIRTOP_QUANT         } from "${projectDir}/modules/mirtop/mirtop_quant"
include { TABLE_MERGE          } from "${projectDir}/modules/r/datatable_merge"
include { EDGER_QC             } from "${projectDir}/modules/r/edger_qc"

include { MIRDEEP2_PIGZ        } from "${projectDir}/modules/mirdeep2/mirdeep2_prepare"
include { MIRDEEP2_MAPPER      } from "${projectDir}/modules/mirdeep2/mirdeep2_mapper"
include { MIRDEEP2_RUN         } from "${projectDir}/modules/mirdeep2/mirdeep2_run"

include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"



///////////////////////

// Function that parses fastp json output file to get total number of reads after trimming

import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    return new JsonSlurper().parseText(json_file.text)
    ?.get('summary')
    ?.get('after_filtering')
    ?.get('total_reads')
    ?.toInteger()
}

String getFastpAdapterSequence(json_file){
    return new JsonSlurper().parseText(json_file.text)
    ?.get('adapter_cutting')
    ?.get('read1_adapter_sequence')
}

def add_suffix(row, suffix) {
    sampleID = "${row[0]}_${suffix}"
    def array = []
    array = [ sampleID, row[1] ]
    return array
}

///////////////////////


// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()



// main workflow
workflow SMRNASEQ {

  if (params.input)     { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }


  // SUBWORKFLOW: Read in samplesheet, validate and stage input files
  INPUT_CHECK(file(params.input)
  )
  .reads
  .dump(tag: 'group')
  .branch {
      meta, fastq ->
          single  : fastq.size() == 1
              return [ meta, fastq.flatten() ]
          multiple: fastq.size() > 1
              return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }


  // Channel setup
  ch_read = ch_fastq.single
  ch_read.map{it -> [it[0].id, it[1]]}.set{read_ch}


//////////////////////////////
  
  // Processes :

  // FastQC Raw reads
  FASTQC_RAW(read_ch)

  // FastP Raw reads
  FASTP(read_ch, params.fastp_known_mirna_adapters, params.save_trimmed_fail, params.save_merged)

  // FastQC Trim reads
  FASTQC_TRIM(FASTP.out.reads)


//////////////////////////////

  // Get adapter sequence and Channel setup for mirtrace inputs

  // This will return with all jsons
  jsons = FASTP.out.json.collect{it[1]}.ifEmpty([])

  // Just pick the first jason file
  json_1 = jsons.map { it[0] }

  // Set adapter sequence
  json_1
    .map { json -> [getFastpAdapterSequence(json)] }
    .set { adapterseq }

  // Combine adapter sequence with reads  and rearrange the columns order
  // and set the channel
  FASTP.out.reads
    .combine(adapterseq)
    .map { id, reads, adapterseq -> [adapterseq, id, reads] }
    .groupTuple()
    .set { ch_mirtrace_inputs }

//////////////////////////////

  // Processes :


  // Mitrace
  MIRTRACE_RUN(ch_mirtrace_inputs)


  // tRNA, cDNA, ncRNA
  // mature, hairpin, genome

  // bowtie indexes

  // Set up trna / cdna / ncrna bowtie2 index channel
  trna_index = Channel.fromPath(params.bowtie2_index_trna+'*')
  cdna_index = Channel.fromPath(params.bowtie2_index_cdna+'*')
  ncrna_index = Channel.fromPath(params.bowtie2_index_ncrna+'*')

  // Set up mature / hairpin / genome bowtie1 index channel
  mature_index = Channel.fromPath(params.bowtie_index_mature+'*')
  hairpin_index = Channel.fromPath(params.bowtie_index_hairpin+'*')
  genome_index = Channel.fromPath(params.ref_fa_indices+'*')


  // Bowtie2
  // Set up input reads channel : tRNA
  FASTP.out.reads
       .dump (tag:'tsux')
       .set { reads_trna }


  // Map tRNA : bowtie2
  MAP_TRNA ( reads_trna, trna_index.collect(), 'tRNA')
  MAP_TRNA.out.unmapped.set { cdna_reads }


  // Map cDNA : bowtie2
  MAP_CDNA ( cdna_reads, cdna_index.collect() , 'cDNA' )
  MAP_CDNA.out.unmapped.set { ncrna_reads }

  
  // Map ncRNA : bowtie2 
  MAP_NCRNA ( ncrna_reads, ncrna_index.collect() , 'ncRNA' )
  MAP_NCRNA.out.unmapped.set { mature_reads }


  // Bowtie1
  // Set up input reads channel : Mature
  mature_reads
       .map { add_suffix(it, "mature") }
       .dump (tag:'msux')
       .set { reads_mirna }

  // Map mature : bowtie1
  BOWTIE_MAP_MATURE ( reads_mirna, mature_index.collect() )

  // Set up input reads channel : Hairpin
  BOWTIE_MAP_MATURE.out.unmapped
        .map { add_suffix(it, "hairpin") }
        .dump (tag:'hsux')
        .set { reads_hairpin }


  // Map hairpin : bowtie1
  BOWTIE_MAP_HAIRPIN ( reads_hairpin, hairpin_index.collect() )


  // EdgeR

  // Samtools sort
  SAMTOOLS_SORT_MATURE(BOWTIE_MAP_MATURE.out.bam, '-O bam', 'bam')
  SAMTOOLS_SORT_HAIRPIN(BOWTIE_MAP_HAIRPIN.out.bam, '-O bam', 'bam')

  // Samtools index ( mature, hairpin, seqcluster )
  SAMTOOLS_INDEX_MATURE(SAMTOOLS_SORT_MATURE.out.sorted_file)
  SAMTOOLS_INDEX_HAIRPIN(SAMTOOLS_SORT_HAIRPIN.out.sorted_file)

  // Samtools stats ( mature, hairpin, seqcluster )
  SAMTOOLS_STATS_MATURE(SAMTOOLS_SORT_MATURE.out.sorted_file)
  SAMTOOLS_STATS_HAIRPIN(SAMTOOLS_SORT_HAIRPIN.out.sorted_file)


  // Set up inputs for EdgeR
  SAMTOOLS_STATS_MATURE.out.idxstat.collect{it[1]}
        .mix(SAMTOOLS_STATS_HAIRPIN.out.idxstat.collect{it[1]})
        .dump(tag:'edger')
        .flatten()
        .collect()
        .set { edger_input }

  // EdgeR
  EDGER_QC ( edger_input )
 
  
  // Seqcluster
  // Set up input reads channel : Seqcluster
  mature_reads
       .map { add_suffix(it, "seqcluster") }
       .dump (tag:'ssux')
       .set { reads_seqcluster }


  // Seqcluster Test here :
  SEQCLUSTER_SEQUENCES ( reads_seqcluster ).collapsed.set { reads_collapsed }

  // Bowtie map seqcluster
  BOWTIE_MAP_SEQCLUSTER ( reads_collapsed, hairpin_index.collect() )

   
  // Mirtop & Table merge
  ch_mirtop_logs = Channel.empty()

  if (params.mirtrace_species){

        //Mirtop quant
        MIRTOP_QUANT(BOWTIE_MAP_SEQCLUSTER.out.bam.collect{it[1]}, params.formatted_hairpin , params.gtf)
        ch_mirtop_logs = MIRTOP_QUANT.out.logs

        // Table merge
        TABLE_MERGE ( MIRTOP_QUANT.out.mirtop_table )
  }


  // Genome
  // Set up input reads : genome
  BOWTIE_MAP_HAIRPIN.out.unmapped
        .map { add_suffix(it, "genome") }
        .dump (tag:'gsux')
        .set { reads_genome }

  // Bowtie map genome
  BOWTIE_MAP_GENOME ( reads_genome, genome_index.collect() )
  SAMTOOLS_SORT_GENOME(BOWTIE_MAP_GENOME.out.bam, '-O bam', 'bam')
  SAMTOOLS_INDEX_GENOME(SAMTOOLS_SORT_GENOME.out.sorted_file)
  SAMTOOLS_STATS_GENOME(SAMTOOLS_SORT_GENOME.out.sorted_file)


  // Mirdeep2
  MIRDEEP2_MAPPER ( mature_reads, genome_index.collect() )
  MIRDEEP2_RUN (params.ref_fa, MIRDEEP2_MAPPER.out.mirdeep2_inputs, params.formatted_hairpin, params.formatted_mature )

  
  // Create channels for multi input files
  ch_multiqc_files = Channel.empty()

  ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_GENOME.out.flagstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_GENOME.out.idxstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_GENOME.out.stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_MATURE.out.flagstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_MATURE.out.idxstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_MATURE.out.stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_HAIRPIN.out.flagstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_HAIRPIN.out.idxstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_HAIRPIN.out.stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(ch_mirtop_logs.collect().ifEmpty([]))
  //ch_multiqc_files = ch_multiqc_files.mix(MIRTOP_QUANT.out.mirtop_logs.collect().ifEmpty([]))

  ch_multiqc_files = ch_multiqc_files.mix(MIRTRACE_RUN.out.mirtrace.collect{it[1]}.ifEmpty([]))


  // Step 41 : MultiQC
  MULTIQC (
      ch_multiqc_files.collect()
  )

}
