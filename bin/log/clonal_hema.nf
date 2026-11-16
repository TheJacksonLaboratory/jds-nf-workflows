def param_log() {

if (params.gen_org != "mouse" && params.gen_org != "human") {
    error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'human' or 'mouse'"
}

// Panel intervals have no default (see config/clonal_hema.config). Check them
// here so a missing file stops the run at startup, not at HsMetrics/Mutect2
// after UMI processing and alignment.
if (!params.ref_fa_bwamem2_indices) {
    error "'--ref_fa_bwamem2_indices' is not set for --gen_org ${params.gen_org}"
}
if (!file(params.ref_fa_bwamem2_indices + '.0123').exists()) {
    error "'--ref_fa_bwamem2_indices': no bwa-mem2 index found at prefix ${params.ref_fa_bwamem2_indices} (expected ${params.ref_fa_bwamem2_indices}.0123)"
}

['target_gatk', 'target_picard', 'bait_picard'].each { p ->
    if (!params[p]) {
        error "'--${p}' is required for --workflow clonal_hema (no default until the panel's final probe files exist)"
    }
    if (!file(params[p]).exists()) {
        error "'--${p}': file not found: ${params[p]}"
    }
}

def msg = """
CLONAL_HEMA PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________

--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--pubdir                        ${params.pubdir}
--project_name                  ${params.project_name}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}

--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--trim_poly_g                   ${params.trim_poly_g}
--trim_poly_x                   ${params.trim_poly_x}
${params.trim_poly_x ? "--poly_x_min_len                ${params.poly_x_min_len}" : ""}

--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--ref_fa_bwamem2_indices        ${params.ref_fa_bwamem2_indices}
--duplex_seq                    ${params.duplex_seq_bool}
--call_min_reads                ${params.call_min_reads}
--min_alt_molecules             ${params.min_alt_molecules}

--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--read_on_target_padding        ${params.read_on_target_padding}

--germline_resource_vcf         ${params.germline_resource_vcf}
--pon_ref                       ${params.pon_ref}
--genotype_pon                  ${params.genotype_pon}
--genotype_germline             ${params.genotype_germline}
--contam_ref                    ${params.contam_ref}
--dbSNP                          ${params.dbSNP}
--cosmic                         ${params.cosmic}
--clinvar                        ${params.clinvar}
--dbNSFP                         ${params.dbNSFP}
--gen_ver                        ${params.gen_ver}
______________________________________________________
"""

def projectInfo = """
Project Directory: ${projectDir}

Command line call:
${workflow.commandLine}
______________________________________________________
"""

msg += projectInfo

log.info msg

return(msg)

}
