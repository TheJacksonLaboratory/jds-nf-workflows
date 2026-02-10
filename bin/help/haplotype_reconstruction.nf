def help(){
println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.
--csv_input | /<FILE_PATH> | CSV delimited sample sheet that controls how samples are processed. The required input headers for the first iteration of haplotype reconstruction are: project_id, finalreport_file, covar_file, cross_type. See the repository wiki (https://github.com/TheJacksonLaboratory/cs-nf-pipelines/wiki) for additional information. 
--rerun | <true> | Options: false and true. Default: true. If this boolean is false, the pipeline accepts a different set of input files. The input fields are: project_id (val), covar_file, cross_file, genoprobs_file, alleleprobs_file, viterbi_file, kinship_file, marker_file, cross_type (val), all originating as results from a prior workflow execution.
--correct_ids | <false> | Options: false and true. Default: false. If this boolean is selected and rerun = false, the pipeline will attempt to correct sample IDs in the input files. This is useful if there sample mixups, or if the phenotype file to be used in the QTL mapping workflow has different sample IDs from what is in the FinalReport files.
--remove_markers | <false> | Options: false and true. Default: false. If this boolean is selected and rerun = false, the pipeline with remove the markers from the marker_file field from the cross object and recalculate genotype and allele probabilities.

'''
}