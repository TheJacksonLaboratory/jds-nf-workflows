def help(){
  println '''

Note: This workflow generates simulated WES reads using NEAT package to generate fastq, 'gold standard' (truth) VCF, and 'gold standard' (truth) BAM files. .

Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--gen_org | mouse | Options: mouse and human.

--fasta | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
        | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'

--target_bed | Mouse: '/projects/omics_share/mouse/GRCm38/supporting_files/benchmarking_data/WES/target_beds/mm10Exome_v4_12-19.1.mm10.baits_merged.bed'
             | Human: '/projects/omics_share/human/GRCh38/supporting_files/benchmarking_data/WES/target_beds/hg38_liftover_agilent_SureSelect_V4_pChrM_PrimaryOnly_probes.bed'

--target_sort_bed | Human: '/projects/omics_share/human/GRCh38/supporting_files/benchmarking_data/WES/target_beds/hg38_liftover_agilent_SureSelect_V4_pChrM_PrimaryOnly_probes_sorted.bed'

'''
}
