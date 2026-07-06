def help(){
  println '''

Note: This workflow generates simulated WES reads using NEAT package to generate fastq, 'gold standard' (truth) VCF, and 'gold standard' (truth) BAM files. .

Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--gen_org | mouse | If a species other than 'mouse' or 'human' is used, all params must be specified manually.
--sampleID | <sampleID> | A sample ID to be used in the read names and output files.

--fasta | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
        | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'
--chrom_list | '1 2 3 4 5 ...' | A list of valid chromosome names to extract from the supplied fasta
--target_bed | <bed_file> | A BED file specifying the target regions for WES simulation. These should correspond to the chromosome names and fasta
--coverage | 60 | The desired coverage for the simulated WES regions.
--off_target_coverage | 0.02 | The simulated coverage for reads not within WES regions.  
--read_length | 150 | The length of each read in the simulated WES data.
--error_rate | 0.001 | The desired sequencing error rate for the simulated WES reads.
--read_type | PE | Options: PE/SE. Paired or Single end reads
--mutation_rate | null | Average mutation rate. The mutation rate model is rescaled to make this the average value. Must be between 0 and 0.3.
--insert_size | 350 | Used when read_type is PE. Paired-end fragment length mean
--insert_size_sd | 30 | Used when read_type is PE. Paired-end fragment length standard deviation

'''
}
