import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

if (params.rsem_aligner != "bowtie2" && params.rsem_aligner != "star") {
  error "'--rsem_aligner': \"${params.rsem_aligner}\" is not valid, supported options are 'star' or 'bowtie2'" 
}

if (!params.bam_input && params.gen_org != "mouse" && params.gen_org != "human" && params.gen_org != "other") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', 'human', or 'other'" 
}

if (params.bam_input && !params.csv_input) {
  error "When `--bam_input` is specified, input must be provided with `--csv_input`." 
}

if (params.umi && params.bam_input) {
  error "`--umi` and `--bam_input` are not compatible." 
}

if (params.umi && params.pdx) {
  error "`--umi` and `--pdx` are not compatible at this time." 
}

if (params.umi && params.rsem_aligner != "star") {
  error "Use of `--umi` data requires alignment with STAR. Specify `--rsem_aligner star` in your command." 
}

if (params.umi && !params.skip_umi_extract && !params.umitools_bc_pattern && !params.umitools_bc_pattern2) {
  error "When `--umi` is specified and `--skip_umi_extract` is not set, at least one of `--umitools_bc_pattern` or `--umitools_bc_pattern2` must also be provided. See workflow Wiki for details."
}

if (params.bam_input && params.rsem_reference_path && params.rsem_reference_name == null) {
  error "When `--rsem_reference_path` is specified, the RSEM reference name must also be provided with `--rsem_reference_name`." 
}

if (params.bam_input && params.rsem_reference_path == null && (params.ref_fa == null || params.ref_gtf == null)) {
  error "When `--rsem_reference_path` is not specified, `--ref_fa` and `--ref_gtf` must be provided to allow an RSEM reference to be built." 
}

// This version of strandedness is used in RSEM_EXPRESSION when bam input is provided, and strandedness must be provided by the user. 
if (params.bam_input && !['forward', 'reverse', 'none'].contains(params.bam_strandedness)) {
    error "When `--bam_input` is specified, `--bam_strandedness` must also be specified. Options are 'none', 'forward', 'reverse'. See parameter details here: https://deweylab.github.io/RSEM/rsem-calculate-expression.html#BASIC-OPTIONS for details" 
}

// This version of strandedness is used in CHECK_STRANDEDNESS to override cases where the tool can't determine stranding. 
if (!params.bam_input && params.strandedness != null && !['reverse_stranded', 'forward_stranded', 'non_stranded'].contains(params.strandedness)) {
  error "'--strandedness': \"${params.strandedness}\" is not valid, supported options are 'reverse_stranded' or 'forward_stranded' or 'non_stranded'" 
}

// Parameter blocks
def logHeader = """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________________________
"""

def baseParams = """
--workflow                              ${params.workflow}
--pubdir                                ${params.pubdir}
-w                                      ${workDir}
-c                                      ${params.config}
--profile                               ${params.profile}

--gen_org                               ${params.gen_org}
--genome_build                          ${params.genome_build}

--read_type                             ${params.read_type}
--csv_input                             ${params.csv_input}
--download_data                         ${params.download_data}
"""

def inputFileParams = """
--sample_folder                         ${params.sample_folder}
--extension                             ${params.extension}
--pattern                               ${params.pattern}
--concat_lanes                          ${params.concat_lanes}
"""

def commonParams = """
--quality_phred                         ${params.quality_phred}
--unqualified_perc                      ${params.unqualified_perc}
--detect_adapter_for_pe                 ${params.detect_adapter_for_pe}
--seed_length                           ${params.seed_length}
--strandedness_ref                      ${params.strandedness_ref}
--strandedness_gtf                      ${params.strandedness_gtf}
--strandedness                          ${params.strandedness}
--rsem_aligner                          ${params.rsem_aligner}
--merge_rna_counts                      ${params.merge_rna_counts}
--skip_read_trimming                    ${params.skip_read_trimming}
"""

def seParams = """
--fragment_length_mean                  ${params.fragment_length_mean}
--fragment_length_sd                    ${params.fragment_length_sd}
"""

def umiParas = """
--umi                                   ${params.umi}
--skip_umi_extract                      ${params.skip_umi_extract}
--umitools_extract_method               ${params.umitools_extract_method}
--umitools_bc_pattern                   ${params.umitools_bc_pattern}
--umitools_bc_pattern2                  ${params.umitools_bc_pattern2}
--umitools_grouping_method              ${params.umitools_grouping_method}
"""

def bowtie2params = """
--rsem_ref_prefix                       ${params.rsem_ref_prefix}
--rsem_ref_files                        ${params.rsem_ref_files}
--picard_dict                           ${params.picard_dict}
--ref_flat                              ${params.ref_flat}
--ribo_intervals                        ${params.ribo_intervals}
"""

def starParams = """
--rsem_ref_prefix                       ${params.rsem_ref_prefix}
--rsem_ref_files                        ${params.rsem_ref_files}
--rsem_star_prefix                      ${params.rsem_star_prefix}
--picard_dict                           ${params.picard_dict}
--ref_flat                              ${params.ref_flat}
--ribo_intervals                        ${params.ribo_intervals}
"""

def pdxParams = """
--pdx                                   ${params.pdx}

Xengsort:             
--ref_fa                                ${params.ref_fa}
--xengsort_host_fasta                   ${params.xengsort_host_fasta}
--xengsort_idx_path                     ${params.xengsort_idx_path}
--xengsort_idx_name                     ${params.xengsort_idx_name}

EBV classifier:               
--classifier_table                      ${params.classifier_table}
"""

def pdxSTARparams = """
Human Specific:
--rsem_ref_prefix_human                 ${params.rsem_ref_prefix_human}
--rsem_star_prefix_human                ${params.rsem_star_prefix_human}
--rsem_ref_files_human                  ${params.rsem_ref_files_human}
--picard_dict_human                     ${params.picard_dict_human}
--ref_flat_human                        ${params.ref_flat_human}
--ribo_intervals_human                  ${params.ribo_intervals_human}

Mouse Specific:
--rsem_ref_prefix_mouse                 ${params.rsem_ref_prefix_mouse}
--rsem_star_prefix_mouse                ${params.rsem_star_prefix_mouse}
--rsem_ref_files_mouse                  ${params.rsem_ref_files_mouse}
--picard_dict_mouse                     ${params.picard_dict_mouse}
--ref_flat_mouse                        ${params.ref_flat_mouse}
--ribo_intervals_mouse                  ${params.ribo_intervals_mouse}
"""

def pdxBowtie2params = """
Human Specific:
--rsem_ref_prefix_human                 ${params.rsem_ref_prefix_human}
--rsem_ref_files_human                  ${params.rsem_ref_files_human}
--picard_dict_human                     ${params.picard_dict_human}
--ref_flat_human                        ${params.ref_flat_human}
--ribo_intervals_human                  ${params.ribo_intervals_human}

Mouse Specific:
--rsem_ref_prefix_mouse                 ${params.rsem_ref_prefix_mouse}
--rsem_ref_files_mouse                  ${params.rsem_ref_files_mouse}
--picard_dict_mouse                     ${params.picard_dict_mouse}
--ref_flat_mouse                        ${params.ref_flat_mouse}
--ribo_intervals_mouse                  ${params.ribo_intervals_mouse}
"""

def bamParams = """
--workflow                              ${params.workflow}
--bam_input                             ${params.bam_input}
--rsem_reference_path                   ${params.rsem_reference_path}
--rsem_reference_name                   ${params.rsem_reference_name}
--ref_fa                                ${params.ref_fa}
--ref_gtf                               ${params.ref_gtf}
--strandedness                          ${params.bam_strandedness}
--read_type                             ${params.read_type}
"""

def projectInfo = """
Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________________________
"""


if (params.bam_input) {

  def msg = logHeader + bamParams

  if (params.read_type == 'SE') {
    msg += seParams
  }

  msg += projectInfo

  log.info msg

} else {

  def msg = logHeader + baseParams

  if (!params.csv_input) {
    msg += inputFileParams
  }

  if (params.read_type == 'SE') {
    msg += seParams
  }

  if (params.umi) {
    msg += umiParas
  }

  msg += commonParams

  if (params.pdx && params.rsem_aligner=='bowtie2') {
    msg += pdxParams + pdxBowtie2params
  }

  if (params.pdx && params.rsem_aligner=='star') {
    msg += pdxParams + pdxSTARparams
  }

  if (!params.pdx && params.rsem_aligner=='bowtie2') {
    msg += bowtie2params
  }

  if (!params.pdx && params.rsem_aligner=='star') {
    msg += starParams
  }

  msg += projectInfo

  log.info msg

  return(msg)

}

}
