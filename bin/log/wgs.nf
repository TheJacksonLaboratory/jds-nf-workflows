import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
<<<<<<< HEAD
    if (params.gen_org != "human" && params.gen_org != "mouse" && params.gen_org != "other") {
        error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', 'human', or 'other'" 
    }

    if (params.gen_org == 'other' && params.run_sv) {
        error "'--gen_org': When '--gen_org other' is specified, '--run_sv' must be false. SV calling is only supported for mouse and human samples."
    }

    if (params.gen_org == 'other' && params.run_mt_calling) {
        error "'--gen_org': When '--gen_org other' is specified, '--run_mt_calling' must be false. MT calling is only supported for mouse and human samples."
    }
    
    if (params.bam_input && !params.csv_input) {
        error "When `--bam_input` is specified, input must be provided with `--csv_input`." 
    }

    def baseParams = """
    --workflow                      ${params.workflow}
    --gen_org                       ${params.gen_org}
    --genome_build                  ${params.genome_build}
    --gen_ver                       ${params.gen_ver}
    --read_type                     ${params.read_type}
    --sample_folder                 ${params.sample_folder}
    --pattern                       ${params.pattern}
    --extension                     ${params.extension}
    --concat_lanes                  ${params.concat_lanes}
    --csv_input                     ${params.csv_input}
    --download_data                 ${params.download_data}
    --merge_inds                    ${params.merge_inds}
    -w                              ${workDir}
    -c                              ${params.config}
    --pubdir                        ${params.pubdir}
    --ref_fa                        ${params.ref_fa}
    --ref_fa_indices                ${params.ref_fa_indices}
    --deduplicate_reads             ${params.deduplicate_reads}
    --split_fastq                   ${params.split_fastq}
    --split_fastq_bin_size          ${params.split_fastq_bin_size}
    --coverage_cap                  ${params.coverage_cap}
    --quality_phred                 ${params.quality_phred}
    --unqualified_perc              ${params.unqualified_perc}
    --detect_adapter_for_pe         ${params.detect_adapter_for_pe}
    --deepvariant                   ${params.deepvariant}
    --run_gvcf                      ${params.run_gvcf}
    --dbSNP                         ${params.dbSNP}
    --snpEff_config                 ${params.snpEff_config}
    --mismatch_penalty              ${params.mismatch_penalty}
    --call_val                      ${params.call_val}
    --ploidy_val                    ${params.ploidy_val}
    """

    def humanParams = """
    --gold_std_indels               ${params.gold_std_indels}
    --phase1_1000G                  ${params.phase1_1000G}
    --dbNSFP                        ${params.dbNSFP}
    --cosmic                        ${params.cosmic}
    --snpEff_config                 ${params.snpEff_config}
    """

    def svParams = """
    --run_sv                        ${params.run_sv}
    --ref_fa_dict                   ${params.ref_fa_dict}
    --smoove_support                ${params.smoove_support}
    --exclude_regions               ${params.exclude_regions}
    --delly_exclusion               ${params.delly_exclusion}
    --delly_mappability             ${params.delly_mappability}
    --cnv_window                    ${params.cnv_window}
    --cnv_min_size                  ${params.cnv_min_size}
    --callRegions                   ${params.callRegions}
    --combined_reference_set        ${params.combined_reference_set}
    --cytoband                      ${params.cytoband}
    --sizemargin                    ${params.sizemargin}
    --sv_slop                       ${params.sv_slop}
    --min_sv_length                 ${params.min_sv_length}
    --cnv_distance_limit            ${params.cnv_distance_limit}
    --ensemblUniqueBed              ${params.ensemblUniqueBed}
    """

    def svHumanParams = """
    --dgv                           ${params.dgv}
    --thousandG                     ${params.thousandG}
    --cosmicUniqueBed               ${params.cosmicUniqueBed}
    --gap                           ${params.gap}
    --dgvBedpe                      ${params.dgvBedpe}
    --thousandGVcf                  ${params.thousandGVcf}
    --svPon                         ${params.svPon}
    --cosmicBedPe                   ${params.cosmicBedPe}
    """

    def svMouseParams = """
    --known_del                     ${params.known_del}
    --known_ins                     ${params.known_ins}
    --known_inv                     ${params.known_inv}
    """

    def mtParams = """
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
    """

    def logHeader = """
    WGS PARAMETER LOG

    --comment: ${params.comment}

    Results Published to: ${params.pubdir}
    ________________________________________________________________________________________
    """

    def svHeader = "________________________________     SV PARAMS     _____________________________________"
    def mtHeader = "________________________________     MT PARAMS     _____________________________________"

    def projectInfo = """
    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ________________________________________________________________________________________
    """

    def msg = logHeader + baseParams

    if (params.gen_org == 'human') {
        msg += humanParams
    }

    if (params.run_sv) {
        msg += "\n" + svHeader + "\n"
        msg += svParams
        if (params.gen_org == 'human') {
            msg += svHumanParams
        } else if (params.gen_org == 'mouse') {
            msg += svMouseParams
        }
    }

    if (params.run_mt_calling) {
        msg += "\n" + mtHeader + "\n"
        msg += mtParams
    }
    
    msg += projectInfo

    log.info msg
=======
if (params.gen_org != "human" && params.gen_org != "mouse" && params.gen_org != "other") {
    error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', 'human', or 'other'" 
}

if (params.gen_org == 'other' && params.run_sv) {
    error "'--gen_org': When '--gen_org other' is specified, '--run_sv' must be false. SV calling is only supported for mouse and human samples."
}

if (params.gen_org == 'other' && params.run_mt_calling) {
    error "'--gen_org': When '--gen_org other' is specified, '--run_mt_calling' must be false. MT calling is only supported for mouse and human samples."
}

def baseParams = """
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--gen_ver                       ${params.gen_ver}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
--merge_inds                    ${params.merge_inds}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
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
--deepvariant                   ${params.deepvariant}
--run_gvcf                      ${params.run_gvcf}
--dbSNP                         ${params.dbSNP}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
"""

def humanParams = """
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--snpEff_config                 ${params.snpEff_config}
"""

def svParams = """
--run_sv                        ${params.run_sv}
--ref_fa_dict                   ${params.ref_fa_dict}
--smoove_support                ${params.smoove_support}
--exclude_regions               ${params.exclude_regions}
--delly_exclusion               ${params.delly_exclusion}
--delly_mappability             ${params.delly_mappability}
--cnv_window                    ${params.cnv_window}
--cnv_min_size                  ${params.cnv_min_size}
--callRegions                   ${params.callRegions}
--combined_reference_set        ${params.combined_reference_set}
--cytoband                      ${params.cytoband}
--sizemargin                    ${params.sizemargin}
--sv_slop                       ${params.sv_slop}
--min_sv_length                 ${params.min_sv_length}
--cnv_distance_limit            ${params.cnv_distance_limit}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
"""

def svHumanParams = """
--dgv                           ${params.dgv}
--thousandG                     ${params.thousandG}
--cosmicUniqueBed               ${params.cosmicUniqueBed}
--gap                           ${params.gap}
--dgvBedpe                      ${params.dgvBedpe}
--thousandGVcf                  ${params.thousandGVcf}
--svPon                         ${params.svPon}
--cosmicBedPe                   ${params.cosmicBedPe}
"""

def svMouseParams = """
--known_del                     ${params.known_del}
--known_ins                     ${params.known_ins}
--known_inv                     ${params.known_inv}
"""

def mtParams = """
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
"""

def logHeader = """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
"""

def svHeader = "________________________________     SV PARAMS     _____________________________________"
def mtHeader = "________________________________     MT PARAMS     _____________________________________"

def projectInfo = """
Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""

def msg = logHeader + baseParams

if (params.gen_org == 'human') {
    msg += humanParams
}

if (params.run_sv) {
    msg += "\n" + svHeader + "\n"
    msg += svParams
    if (params.gen_org == 'human') {
        msg += svHumanParams
    } else if (params.gen_org == 'mouse') {
        msg += svMouseParams
    }
}

if (params.run_mt_calling) {
    msg += "\n" + mtHeader + "\n"
    msg += mtParams
}

msg += projectInfo

log.info msg

return(msg)
>>>>>>> remotes/origin/release/0.9.3
}
