def help(){
  println '''

Note: This workflow generates simulated RNA reads from a set of transcripts at specified abundances. Each transcript has a predefined abundance in the output.

Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--gen_org | mouse | Options: mouse and human.

--fa_cds | Mouse: 'http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz' 
         | Human: 'http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz'
         | All transcript coding sequences resulting from Ensembl genes.

'''
}
