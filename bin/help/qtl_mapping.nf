def help(){
println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.
--csv_input | /<FILE_PATH> | CSV delimited sample sheet that controls how samples are processed. The required input header is: sampleID,outputID,gvcf. See the repository wiki (https://github.com/TheJacksonLaboratory/cs-nf-pipelines/wiki) for additional information. 
--gen_org | mouse | Options: mouse, human, other.
--genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCh38. If other, this is ignored.

!!!ADJUST AS NEEDED!!!

'''
}