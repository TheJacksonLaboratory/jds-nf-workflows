import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log() {

if (params.gen_org != "mouse" && params.gen_org != "human") {
    error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'human' or 'mouse'" 
}

def baseParams = """
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--pubdir                        ${params.pubdir}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--deduplicate_reads             ${params.deduplicate_reads}
--split_fastq                   ${params.split_fastq}
--split_fastq_bin_size          ${params.split_fastq_bin_size}
--coverage_cap                  ${params.coverage_cap}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--trim_poly_g                   ${params.trim_poly_g}
--trim_poly_x                   ${params.trim_poly_x}
${params.trim_poly_x ? "--poly_x_min_len                ${params.poly_x_min_len}" : ""}
--xengsort_host_fasta           ${params.xengsort_host_fasta}
--xengsort_idx_path             ${params.xengsort_idx_path}
--xengsort_idx_name             ${params.xengsort_idx_name}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--ref_fa_dict                   ${params.ref_fa_dict}
--combined_reference_set        ${params.combined_reference_set}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
--chrom_contigs                 ${params.chrom_contigs}
--chrom_intervals               ${params.chrom_intervals}
--excludeIntervalList           ${params.excludeIntervalList}
--intervalListBed               ${params.intervalListBed}
--lancet_beds_directory         ${params.lancet_beds_directory}
--strelka_config                ${params.strelka_config}
--vep_cache_directory           ${params.vep_cache_directory}
--vep_fasta                     ${params.vep_fasta}
--cytoband                      ${params.cytoband}
--callRegions                   ${params.callRegions}
"""

def humanParams = """
--pdx                           ${params.pdx}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--hapmap                        ${params.hapmap}
--omni                          ${params.omni}
--pon_bed                       ${params.pon_bed}
--germline_filtering_vcf        ${params.germline_filtering_vcf}
--gripss_pon                    ${params.gripss_pon}
--msisensor_model               ${params.msisensor_model}
--cosmic_cgc                    ${params.cosmic_cgc}
--cosmic_cancer_resistance_muts ${params.cosmic_cancer_resistance_muts}
--ensembl_entrez                ${params.ensembl_entrez}
--dgv                           ${params.dgv}
--thousandG                     ${params.thousandG}
--cosmicUniqueBed               ${params.cosmicUniqueBed}
--cancerCensusBed               ${params.cancerCensusBed}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
--gap                           ${params.gap}
--dgvBedpe                      ${params.dgvBedpe}
--thousandGVcf                  ${params.thousandGVcf}
--svPon                         ${params.svPon}
--cosmicBedPe                   ${params.cosmicBedPe}
--na12878_bam                   ${params.na12878_bam}
--na12878_bai                   ${params.na12878_bai}
--na12878_sampleName            ${params.na12878_sampleName}
"""

def mouseParams = """
--delly_exclusion               ${params.delly_exclusion}
--delly_mappability             ${params.delly_mappability}
--cnv_window                    ${params.cnv_window}
--cnv_min_size                  ${params.cnv_min_size}
--cnv_germline_prob             ${params.cnv_germline_prob}
--known_del                     ${params.known_del}
--known_ins                     ${params.known_ins}
--known_inv                     ${params.known_inv}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
--gap                           ${params.gap}
--exclude_list                  ${params.exclude_list}
--proxy_normal_bam              ${params.proxy_normal_bam}
--proxy_normal_bai              ${params.proxy_normal_bai}
--proxy_normal_sampleName       ${params.proxy_normal_sampleName}
"""


def mtParams = """
--mt_contig_name                       ${params.mt_contig_name}
--mt_fasta                             ${params.mt_fasta}
--mt_genome                            ${params.mt_genome}
--mt_shifted_fasta                     ${params.mt_shifted_fasta}
--shift_back_chain                     ${params.shift_back_chain}
--mt_fasta_index                       ${params.mt_fasta_index}
--mt_shifted_fasta_index               ${params.mt_shifted_fasta_index}
--exclusion_sites                      ${params.exclusion_sites}
--non_control_region_interval_list     ${params.non_control_region_interval_list}
--control_region_shifted_interval_list ${params.control_region_shifted_interval_list}
--gen_ver                              ${params.gen_ver}
--dbNSFP                               ${params.dbNSFP}
--snpEff_config                        ${params.snpEff_config}
--cosmic                               ${params.cosmic}
--cosmic_index                         ${params.cosmic_index}
"""


def logHeader = """
PTA PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
"""

def projectInfo = """
Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

def mtHeader = "____________________ MT PARAMS _______________________"

def msg = logHeader + baseParams

if (params.gen_org == 'human') {
  msg += humanParams
} else if (params.gen_org == 'mouse') {
  msg += mouseParams
}

if (params.run_mt_calling) {
  msg += "\n" + mtHeader + "\n"
  msg += mtParams
}

msg += projectInfo

log.info msg

return(msg)

}