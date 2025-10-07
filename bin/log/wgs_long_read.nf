import Logos

logo = new Logo()
println '\n'
println logo.show()


def param_log(){
if (params.gen_org != "human" && params.gen_org != "mouse" && params.gen_org != "other") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', 'human', or 'other'" 
}


if (params.pbmode != "CCS" && params.pbmode != "CLR") {
    error "'--pbmode': \"${params.pbmode}\" is not valid, supported options are 'CCS' or 'CLR'" 
}

if (params.deepvariant != true) {
    error "'--deepvariant': \"${params.deepvariant}\" must be true. Workflow is currently only supports deepvariant." 
}


if (params.gen_org=='human')
log.info """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--gen_ver                       ${params.gen_ver}
--csv_input                     ${params.csv_input}
--merge_inds                    ${params.merge_inds}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--data_type                     ${params.data_type}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--ref_fa                        ${params.ref_fa}
--chrom_contigs                 ${params.chrom_contigs}
--minimap2_index                ${params.minimap2_index}
--pbmode                        ${params.pbmode}
--deepvariant                   ${params.deepvariant}
--deepvariant_model_type        ${params.deepvariant_model_type}
--run_gvcf                      ${params.run_gvcf}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--snpEff_config                 ${params.snpEff_config}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--cosmic_index                  ${params.cosmic_index}
--pbsv_tandem                   ${params.pbsv_tandem}
--tandem_repeats                ${params.tandem_repeats}
--vep_cache_directory           ${params.vep_cache_directory}
--vep_fasta                     ${params.vep_fasta}
--min_sv_length                 ${params.min_sv_length}
--sv_slop                       ${params.sv_slop}
--sizemargin                    ${params.sizemargin}
--multiqc_config                ${params.multiqc_config}


Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""
else
log.info """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--gen_ver                       ${params.gen_ver}
--csv_input                     ${params.csv_input}
--merge_inds                    ${params.merge_inds}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--data_type                     ${params.data_type}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--ref_fa                        ${params.ref_fa}
--chrom_contigs                 ${params.chrom_contigs}
--minimap2_index                ${params.minimap2_index}
--pbmode                        ${params.pbmode}
--deepvariant                   ${params.deepvariant}
--deepvariant_model_type        ${params.deepvariant_model_type}
--run_gvcf                      ${params.run_gvcf}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--snpEff_config                 ${params.snpEff_config}
--pbsv_tandem                   ${params.pbsv_tandem}
--tandem_repeats                ${params.tandem_repeats}
--vep_cache_directory           ${params.vep_cache_directory}
--vep_fasta                     ${params.vep_fasta}
--min_sv_length                 ${params.min_sv_length}
--sv_slop                       ${params.sv_slop}
--sizemargin                    ${params.sizemargin}
--multiqc_config                ${params.multiqc_config}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""

}
