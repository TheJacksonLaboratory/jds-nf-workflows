def help(){
  println '''

Note: This workflow generates simulated WGS reads using NEAT package to generate fastq, 'gold standard' (truth) VCF, and 'gold standard' (truth) BAM files. .

Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--gen_org | mouse | Options: mouse and human.

--fasta | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
        | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'

'''
}
