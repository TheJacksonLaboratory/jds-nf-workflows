import Logos

logo = new Logo()
println '\n'
println logo.show()

if (params.gen_org != "human" && params.gen_org != "mouse") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', or 'human'" 
}

def param_log(){
if (params.gen_org=='mouse')
log.info """
mtDNA VARIANT CALLING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                              ${params.workflow}
-w                                      ${workDir}
-c                                      ${params.config}
--csv_input                             ${params.csv_input}
--gen_org                               ${params.gen_org}
--genome_build                          ${params.genome_build}
--mt_contig_name                        ${params.mt_contig_name}
--mt_fasta                              ${params.mt_fasta}
--mt_genome                             ${params.mt_genome}
--mt_shifted_fasta                      ${params.mt_shifted_fasta}
--shift_back_chain                      ${params.shift_back_chain}
--mt_fasta_index                        ${params.mt_fasta_index}
--mt_shifted_fasta_index                ${params.mt_shifted_fasta_index}
--max_allele_count                      ${params.max_allele_count}
--exclusion_sites                       ${params.exclusion_sites}
--non_control_region_interval_list      ${params.non_control_region_interval_list}
--control_region_shifted_interval_list  ${params.control_region_shifted_interval_list}
--detection_limit                       ${params.detection_limit}
--mapQ                                  ${params.mapQ}
--baseQ                                 ${params.baseQ}
--gen_ver                               ${params.gen_ver}
--dbSNP                                 ${params.dbSNP}
--dbSNP_index                           ${params.dbSNP_index}
--snpEff_config                         ${params.snpEff_config}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
log.info """
mtDNA VARIANT CALLING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                              ${params.workflow}
-w                                      ${workDir}
-c                                      ${params.config}
--csv_input                             ${params.csv_input}
--gen_org                               ${params.gen_org}
--genome_build                          ${params.genome_build}
--mt_contig_name                        ${params.mt_contig_name}
--mt_fasta                              ${params.mt_fasta}
--mt_genome                             ${params.mt_genome}
--mt_shifted_fasta                      ${params.mt_shifted_fasta}
--shift_back_chain                      ${params.shift_back_chain}
--mt_fasta_index                        ${params.mt_fasta_index}
--mt_shifted_fasta_index                ${params.mt_shifted_fasta_index}
--max_allele_count                      ${params.max_allele_count}
--exclusion_sites                       ${params.exclusion_sites}
--non_control_region_interval_list      ${params.non_control_region_interval_list}
--control_region_shifted_interval_list  ${params.control_region_shifted_interval_list}
--detection_limit                       ${params.detection_limit}
--mapQ                                  ${params.mapQ}
--baseQ                                 ${params.baseQ}
--gen_ver                               ${params.gen_ver}
--dbSNP                                 ${params.dbSNP}
--dbSNP_index                           ${params.dbSNP_index}
--snpEff_config                         ${params.snpEff_config}
--cosmic                                ${params.cosmic}
--cosmic_index                          ${params.cosmic_index}
--dbNSFP                                ${params.dbNSFP}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}



