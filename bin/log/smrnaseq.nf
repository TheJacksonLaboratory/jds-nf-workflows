import Logos

logo = new Logo()
println '\n'
println logo.show()

if (!(params.gen_org in ['human', 'mouse'])) {
  log.error "This workflow supports only --gen_org mouse or --gen_org human. Please retry with correct option."
  System.exit(1)
}

if (!params.csv_input) {
  log.error "The --csv_input parameter is required. Please provide a CSV input file."
  System.exit(1)
}


def param_log(){
if (params.gen_org == 'human') {
  if (params.mirtrace_species != 'hsa') {
    log.error "--mirtrace_species should be set to 'hsa' when '--gen_org human'"
    System.exit(1)
  }


  log.info """
  SMRNASEQ PARAMETER LOG

  --comment: ${params.comment}

  Results Published to: ${params.pubdir}
  ______________________________________________________
  --workflow                      ${params.workflow}
  --gen_org                       ${params.gen_org}
  --genome_build                  ${params.genome_build}
  --read_type                     ${params.read_type}
  --csv_input                     ${params.csv_input}
  -w                              ${workDir}
  -c                              ${params.config}
  --pubdir                        ${params.pubdir}
  --gtf                           ${params.gtf}
  --ref_fa                        ${params.ref_fa}
  --ref_fa_indices                ${params.ref_fa_indices}
  --bowtie_index_mature           ${params.bowtie_index_mature}
  --formatted_mature              ${params.formatted_mature}
  --bowtie_index_hairpin          ${params.bowtie_index_hairpin}
  --formatted_hairpin             ${params.formatted_hairpin}
  --bowtie2_index_trna            ${params.bowtie2_index_trna}
  --bowtie2_index_cdna            ${params.bowtie2_index_cdna}
  --bowtie2_index_ncrna           ${params.bowtie2_index_ncrna}
  --mirtrace_species              ${params.mirtrace_species}
  --clip_r1                       ${params.clip_r1} 
  --three_prime_clip_r1           ${params.three_prime_clip_r1} 
  --three_prime_adapter           ${params.three_prime_adapter}
  --trim_fastq                    ${params.trim_fastq}
  --fastp_min_length              ${params.fastp_min_length}
  --fastp_max_length              ${params.fastp_max_length}
  --adapter_fasta                 ${params.adapter_fasta}
  --tmpdir                        ${params.tmpdir}

  Project Directory: ${projectDir}

  Command line call: 
  ${workflow.commandLine}
  ______________________________________________________
"""
} else {
  if (params.mirtrace_species != 'mmu') {
    log.error "--mirtrace_species should be set to 'mmu' when '--gen_org mouse'"
    System.exit(1)
  }

  log.info """
  SMRNASEQ PARAMETER LOG

  --comment: ${params.comment}

  Results Published to: ${params.pubdir}
  ______________________________________________________
  --workflow                      ${params.workflow}
  --gen_org                       ${params.gen_org}
  --genome_build                  ${params.genome_build}
  --read_type                     ${params.read_type}
  --csv_input                     ${params.csv_input}
  -w                              ${workDir}
  -c                              ${params.config}
  --pubdir                        ${params.pubdir}
  --gtf                           ${params.gtf}
  --ref_fa                        ${params.ref_fa}
  --ref_fa_indices                ${params.ref_fa_indices}
  --bowtie_index_mature           ${params.bowtie_index_mature}
  --formatted_mature              ${params.formatted_mature}
  --bowtie_index_hairpin          ${params.bowtie_index_hairpin}
  --formatted_hairpin             ${params.formatted_hairpin}
  --bowtie2_index_trna            ${params.bowtie2_index_trna}
  --bowtie2_index_cdna            ${params.bowtie2_index_cdna}
  --bowtie2_index_ncrna           ${params.bowtie2_index_ncrna}
  --mirtrace_species              ${params.mirtrace_species}
  --clip_r1                       ${params.clip_r1}
  --three_prime_clip_r1           ${params.three_prime_clip_r1}
  --three_prime_adapter           ${params.three_prime_adapter}
  --trim_fastq                    ${params.trim_fastq}
  --fastp_min_length              ${params.fastp_min_length}
  --fastp_max_length              ${params.fastp_max_length}
  --adapter_fasta                 ${params.adapter_fasta}
  --tmpdir                        ${params.tmpdir}

  Project Directory: ${projectDir}

  Command line call: 
  ${workflow.commandLine}
  ______________________________________________________
  """

  }
}
