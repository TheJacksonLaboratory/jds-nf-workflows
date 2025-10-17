import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
if (params.gen_org=='human')
  log.info """
SMRNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--input                         ${params.input}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--gtf                           ${params.gtf}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--bowtie2_index_mature          ${params.bowtie2_index_mature}
--formatted_mature              ${params.formatted_mature}
--bowtie2_index_hairpin         ${params.bowtie2_index_hairpin}
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
--multiqc_config              	${params.multiqc_config}
--tmpdir                        ${params.tmpdir}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
log.info """
SMRNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--input                         ${params.input}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--gtf                           ${params.gtf}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--bowtie2_index_mature          ${params.bowtie2_index_mature}
--formatted_mature              ${params.formatted_mature}
--bowtie2_index_hairpin         ${params.bowtie2_index_hairpin}
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
--multiqc_config                ${params.multiqc_config}
--tmpdir                        ${params.tmpdir}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

}
