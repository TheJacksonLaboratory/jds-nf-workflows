def help(){
  println '''
Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. 
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--read_type | <string> | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,bam". See the repository wiki for an example file. BAM files can either be absolute paths to local files, or URLs to remote files.

'''
}
