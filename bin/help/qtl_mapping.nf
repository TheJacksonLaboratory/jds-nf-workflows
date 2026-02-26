def help(){
println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.
--csv_input | /<FILE_PATH> | CSV delimited sample sheet that controls how samples are processed. The required input header is: project_id, cross_file, genoprobs_file, alleleprobs_file, kinship_file, covar_file, pheno_file. See the repository wiki (https://github.com/TheJacksonLaboratory/cs-nf-pipelines/wiki) for additional information. 
--n_perms | <INT> | The number of permutations run to determinine the signficance threshold for QTL mapped for each trait. Default is 1000.
--primary_chrom_bed | <FILE_PATH> | BED file that contains the coordinates of the primary chromosome. This is used in generating plots that summarize QTL effects.  
'''
}