import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

if (params.rsem_aligner != "bowtie2" && params.rsem_aligner != "star") {
  error "'--rsem_aligner': \"${params.rsem_aligner}\" is not valid, supported options are 'bowtie2' or 'star'" 
}

if (!params.bam_input && params.gen_org != "mouse" && params.gen_org != "human") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse' or 'human'" 
}

if (params.bam_input && !params.csv_input) {
  error "When `--bam_input` is specified, input must be provided with `--csv_input`." 
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


if (params.bam_input)
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                             ${params.workflow}
--bam_input                            ${params.bam_input}
--rsem_reference_path                  ${params.rsem_reference_path}
--rsem_reference_name                  ${params.rsem_reference_name}
--ref_fa                               ${params.ref_fa}
--ref_gtf                              ${params.ref_gtf}
--strandedness                         ${params.bam_strandedness}
--read_type                            ${params.read_type}
--fragment_length_mean (SE only)       ${params.fragment_length_mean}
--fragment_length_sd (SE only)         ${params.fragment_length_sd}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________

"""



else if (params.pdx && params.rsem_aligner=='bowtie2')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                   ${params.workflow}
--gen_org                    ${params.gen_org}
--genome_build               ${params.genome_build}
--read_type                  ${params.read_type}
--sample_folder              ${params.sample_folder}
--extension                  ${params.extension}
--pattern                    ${params.pattern}
--concat_lanes               ${params.concat_lanes}
--csv_input                  ${params.csv_input}
--download_data              ${params.download_data}
--pubdir                     ${params.pubdir}
-w                           ${workDir}
--keep_intermediate          ${params.keep_intermediate}
-c                           ${params.config}
--seed_length                ${params.seed_length}
--quality_phred              ${params.quality_phred}
--unqualified_perc           ${params.unqualified_perc}
--detect_adapter_for_pe      ${params.detect_adapter_for_pe}

--pdx                        ${params.pdx}
--ref_fa                     ${params.ref_fa}
--xengsort_host_fasta        ${params.xengsort_host_fasta}
--xengsort_idx_path          ${params.xengsort_idx_path}
--xengsort_idx_name          ${params.xengsort_idx_name}

--strandedness_ref           ${params.strandedness_ref}
--strandedness_gtf           ${params.strandedness_gtf}
--strandedness               ${params.strandedness}

--rsem_aligner               ${params.rsem_aligner}
--merge_rna_counts           ${params.merge_rna_counts}
--skip_read_trimming         ${params.skip_read_trimming}

Human specific files: 
--rsem_ref_prefix_human      ${params.rsem_ref_prefix_human}
--rsem_ref_files_human       ${params.rsem_ref_files_human}
--picard_dict_human          ${params.picard_dict_human}
--ref_flat_human             ${params.ref_flat_human}
--ribo_intervals_human       ${params.ribo_intervals_human}
--classifier_table           ${params.classifier_table}

Mouse specific files: 
--rsem_ref_prefix_mouse      ${params.rsem_ref_prefix_mouse}
--rsem_ref_files_mouse       ${params.rsem_ref_files_mouse}
--picard_dict_mouse          ${params.picard_dict_mouse}
--ref_flat_mouse             ${params.ref_flat_mouse}
--ribo_intervals_mouse       ${params.ribo_intervals_mouse}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________

"""
else if (params.pdx && params.rsem_aligner=='star')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                   ${params.workflow}
--gen_org                    ${params.gen_org}
--genome_build               ${params.genome_build}
--read_type                  ${params.read_type}
--sample_folder              ${params.sample_folder}
--extension                  ${params.extension}
--pattern                    ${params.pattern}
--concat_lanes               ${params.concat_lanes}
--csv_input                  ${params.csv_input}
--download_data              ${params.download_data}
--pubdir                     ${params.pubdir}
-w                           ${workDir}
--keep_intermediate          ${params.keep_intermediate}
-c                           ${params.config}
--quality_phred              ${params.quality_phred}
--unqualified_perc           ${params.unqualified_perc}
--detect_adapter_for_pe      ${params.detect_adapter_for_pe}
--seed_length                ${params.seed_length}

--pdx                        ${params.pdx}
--xengsort_host_fasta        ${params.xengsort_host_fasta}
--xengsort_idx_path          ${params.xengsort_idx_path}
--xengsort_idx_name          ${params.xengsort_idx_name}

--strandedness_ref           ${params.strandedness_ref}
--strandedness_gtf           ${params.strandedness_gtf}
--strandedness               ${params.strandedness}

--rsem_aligner               ${params.rsem_aligner}
--merge_rna_counts           ${params.merge_rna_counts}
--skip_read_trimming         ${params.skip_read_trimming}

Human specific files: 
--rsem_ref_prefix_human      ${params.rsem_ref_prefix_human}
--rsem_ref_files_human       ${params.rsem_ref_files_human}
--rsem_star_prefix_human     ${params.rsem_star_prefix_human}
--picard_dict_human          ${params.picard_dict_human}
--ref_flat_human             ${params.ref_flat_human}
--ribo_intervals_human       ${params.ribo_intervals_human}
--classifier_table           ${params.classifier_table}

Mouse specific files: 
--rsem_ref_prefix_mouse      ${params.rsem_ref_prefix_mouse}
--rsem_ref_files_mouse       ${params.rsem_ref_files_mouse}
--rsem_star_prefix_mouse     ${params.rsem_star_prefix_mouse}
--picard_dict_mouse          ${params.picard_dict_mouse}
--ref_flat_mouse             ${params.ref_flat_mouse}
--ribo_intervals_mouse       ${params.ribo_intervals_mouse}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________

"""
else if (params.gen_org=='human' && params.rsem_aligner=='bowtie2')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                ${params.workflow}
--gen_org                 ${params.gen_org}
--genome_build            ${params.genome_build}
--read_type               ${params.read_type}
--sample_folder           ${params.sample_folder}
--extension               ${params.extension}
--pattern                 ${params.pattern}
--concat_lanes            ${params.concat_lanes}
--csv_input               ${params.csv_input}
--download_data           ${params.download_data}
--pubdir                  ${params.pubdir}
-w                        ${workDir}
--keep_intermediate       ${params.keep_intermediate}
-c                        ${params.config}
--quality_phred           ${params.quality_phred}
--unqualified_perc        ${params.unqualified_perc}
--detect_adapter_for_pe   ${params.detect_adapter_for_pe}
--strandedness_ref        ${params.strandedness_ref}
--strandedness_gtf        ${params.strandedness_gtf}
--strandedness            ${params.strandedness}
--seed_length             ${params.seed_length}
--rsem_ref_prefix         ${params.rsem_ref_prefix}
--rsem_ref_files          ${params.rsem_ref_files}
--rsem_aligner            ${params.rsem_aligner}
--merge_rna_counts        ${params.merge_rna_counts}
--skip_read_trimming      ${params.skip_read_trimming}
--picard_dict             ${params.picard_dict}
--ref_flat                ${params.ref_flat}
--ribo_intervals          ${params.ribo_intervals}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.gen_org=='human' && params.rsem_aligner=='star')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                ${params.workflow}
--gen_org                 ${params.gen_org}
--genome_build            ${params.genome_build}
--read_type               ${params.read_type}
--sample_folder           ${params.sample_folder}
--extension               ${params.extension}
--pattern                 ${params.pattern}
--concat_lanes            ${params.concat_lanes}
--csv_input               ${params.csv_input}
--download_data           ${params.download_data}
--pubdir                  ${params.pubdir}
-w                        ${workDir}
--keep_intermediate       ${params.keep_intermediate}
-c                        ${params.config}
--quality_phred           ${params.quality_phred}
--unqualified_perc        ${params.unqualified_perc}
--detect_adapter_for_pe   ${params.detect_adapter_for_pe}
--strandedness_ref        ${params.strandedness_ref}
--strandedness_gtf        ${params.strandedness_gtf}
--strandedness            ${params.strandedness}
--seed_length             ${params.seed_length}
--rsem_ref_prefix         ${params.rsem_ref_prefix}
--rsem_ref_files          ${params.rsem_ref_files}
--rsem_aligner            ${params.rsem_aligner}
--merge_rna_counts        ${params.merge_rna_counts}
--skip_read_trimming      ${params.skip_read_trimming}
--rsem_star_prefix        ${params.rsem_star_prefix}
--picard_dict             ${params.picard_dict}
--ref_flat                ${params.ref_flat}
--ribo_intervals          ${params.ribo_intervals}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.gen_org=='mouse' && params.rsem_aligner=='bowtie2')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--extension                     ${params.extension}
--pattern                       ${params.pattern}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
--pubdir                        ${params.pubdir}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--strandedness_ref              ${params.strandedness_ref}
--strandedness_gtf              ${params.strandedness_gtf}
--strandedness                  ${params.strandedness}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--merge_rna_counts              ${params.merge_rna_counts}
--skip_read_trimming            ${params.skip_read_trimming}
--picard_dict                   ${params.picard_dict}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else if (params.gen_org=='mouse' && params.rsem_aligner=='star')
log.info """
RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--extension                     ${params.extension}
--pattern                       ${params.pattern}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
--pubdir                        ${params.pubdir}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--strandedness_ref              ${params.strandedness_ref}
--strandedness_gtf              ${params.strandedness_gtf}
--strandedness                  ${params.strandedness}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--merge_rna_counts              ${params.merge_rna_counts}
--skip_read_trimming            ${params.skip_read_trimming}
--rsem_star_prefix              ${params.rsem_star_prefix}
--picard_dict                   ${params.picard_dict}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

else error "Invalid parameters in '--gen_org': ${params.gen_org} and/or in '--rsem_aligner': ${params.rsem_aligner}. Supported options are 'mouse' or 'human' and 'bowtie2' or 'star'."

}
