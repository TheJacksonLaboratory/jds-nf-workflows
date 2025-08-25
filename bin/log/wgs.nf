import Logos

logo = new Logo()
println '\n'
println logo.show()


def param_log(){
if (params.gen_org != "human" && params.gen_org != "mouse" && params.gen_org != "other") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', 'human', or 'other'" 
}

if (params.gen_org == 'other' && params.run_sv) {
  error "'--gen_org': When '--gen_org other' is specified, '--run_sv' must be false. SV calling is only supported for mouse and human samples."
}


if (params.gen_org=='human' && !params.run_sv)
log.info """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
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
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--snpEff_config                 ${params.snpEff_config}


Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""
else if ((params.gen_org=='mouse' || params.gen_org == 'other') && !params.run_sv)
log.info """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
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

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""
else if (params.gen_org=='human' && params.run_sv)
log.info """
WGS WITH SV CALLING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________     WGS PARAMS     ____________________________________
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
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--snpEff_config                 ${params.snpEff_config}

________________________________     SV PARAMS     _____________________________________
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
--dgv                           ${params.dgv}
--thousandG                     ${params.thousandG}
--cosmicUniqueBed               ${params.cosmicUniqueBed}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
--gap                           ${params.gap}
--dgvBedpe                      ${params.dgvBedpe}
--thousandGVcf                  ${params.thousandGVcf}
--svPon                         ${params.svPon}
--cosmicBedPe                   ${params.cosmicBedPe}
--min_sv_length                 ${params.min_sv_length}
--cnv_distance_limit            ${params.cnv_distance_limit}
--sv_slop                       ${params.sv_slop}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""
else if (params.gen_org=='mouse' && params.run_sv)
log.info """
WGS WITH SV CALLING PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________     WGS PARAMS     ____________________________________
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

________________________________     SV PARAMS     _____________________________________
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
--known_del                     ${params.known_del}
--known_ins                     ${params.known_ins}
--known_inv                     ${params.known_inv}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
--min_sv_length                 ${params.min_sv_length}
--cnv_distance_limit            ${params.cnv_distance_limit}
--sv_slop                       ${params.sv_slop}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""
}
