def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--input | /<PATH> | The path to the samplesheet file that contains all the samples to be run by the pipeline. - For samplesheet file format, please see : 'https://nf-co.re/smrnaseq/2.2.4/docs/usage'   
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).

--gen_org | mouse | Options: mouse and human.
--genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCh38.

--gtf      | GFF/GTF file with coordinates positions of precursor and miRNAs. 
           | miRBase .gff3 file, typically downloaded from https://mirbase.org/ftp/CURRENT/genomes/                                                        

--ref_fa         | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bowtie/genome.fa' 
                 | Human: '/projects/omics_share/human/GRCh38/genome/indices/ensembl/v109/bowtie/genome.fa'
--ref_fa_indices | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bowtie/' | The default value for mm10. 
                 | Human: '/projects/omics_share/human/GRCh38/genome/indices/ensembl/v109/bowtie/'
                 | Pre-compiled BOWTIE index files, points to human reference when --gen_org human.

--bowtie_index_mature  | <PATH> | bowtie index for mature miRNAs
--formatted_mature     | <PATH> | The path to formatted mature miRNAs file
--bowtie_index_hairpin | <PATH> | bowtie index for miRNAs precursors
--formatted_hairpin    | <PATH> | The path to formatted miRNAs precursors file

--bowtie2_index_trna  | bowtie2 index for tRNA contamination database
--bowtie2_index_cdna  | bowtie2 index for cDNA contamination database
--bowtie2_index_ncrna | bowtie2 index for ncRNA contamination database

--mirtrace_species | Species for miRTrace | Example values: hsa for human, mmu for mouse                           

--clip_r1 | integer | The number of basepairs to remove from the 5’ end of read 1
--three_prime_clip_r1 | integer | The number of basepairs to remove from the 3’ end of read 1 AFTER adapter/quality trimming has been performed
--three_prime_adapter | 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'  | Sequencing adapter sequence to use for trimming
--trim_fastq | true | Trim FastQ files
--fastp_min_length | 17 | Minimum filter length for raw reads
--fastp_max_length | 100 | Maximum filter length for raw reads 
--adapter_fasta | known_adapters.fa | Fasta with known miRNA adapter sequences for adapter trimming
--tmpdir  | /<PATH> | Temporary directory to store temp files.  
'''
}

