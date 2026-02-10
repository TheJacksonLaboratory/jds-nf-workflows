def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. 
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specified, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 
--gen_org | mouse | Options: mouse
--genome_build | 'GRCm39' | Mouse specific. Options: GRCm39
--library_type | 'seqwell' | Type of DNA sequencing library prep used. Options: 'seqwell', 'ddRADseq'. This parameter is used to select the appropriate sequence QC process.
--ref_fa | Mouse: '/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa' 
         | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. 
--ref_fa_indices | Mouse: '/projects/omics_share/mouse/GRCm39/genome/indices/ensembl/v105/bwa/Mus_musculus.GRCm39.dna.primary_assembly.fa'
                 | Pre-compiled BWA index files, points to mouse reference files.
--ref_haps_dir | '/projects/omics_share/mouse/GRCm39/supporting_files/lcwgs_hr'
               | Directory containing reference haplotypes for genotype imputation with QUILT.
--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.

--mismatch_penalty | 8 | The BWA penalty for a mismatch.
--gridfile | '/projects/omics_share/mouse/GRCm39/supporting_files/lcwgs_hr/interp_1M_physical_grid.csv'
           | Uniformly spaced marker grid used to filter high-quality SNPs from QUILT genotype imputation
--covar_file | /<PATH> | The path to the .csv file used specify covariates in an R/qtl2 QTL mapping experiment
             | See https://kbroman.org/qtl2/assets/vignettes/input_files.html#Detailed_specifications_for_each_cross_type for formatting covariates for different cross types.
--cross_type | 'do' | Options: do, cc, het3, bxd, genail4, genail8. Parameter specifying the cross type for R/qtl2 QTL mapping.
--downsample | 'false' | Options: true, false. If true, downsample the aligned reads to coverage specified.
--downsampling_coverage_csv | /<PATH> | Path to a CSV file specifying the target coverage for downsampling each sample. The CSV file is just one unheadered column with numeric values (X coverage) in each row.
--smooth_window | 200 | Number of markers to smooth over when calculating genotype probabilities in R/qtl2.
--interp_250k_gridfile | '/projects/omics_share/mouse/GRCm39/supporting_files/lcwgs_hr/interp_0.25M_physical_grid.csv'
                             | Uniformly spaced marker grid used to interpolate genotype probabilities to a standard 250k marker grid for R/qtl2 QTL mapping.
'''  
}

