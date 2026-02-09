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

--gen_org | mouse | Options: mouse and human.
--genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCm38.

--bam_input | false | Options: false, true. If specified, use BAM file input. `--csv_input` is required. See the wiki page for details. 

--umi | false | Options: false, true. If specified, UMI-tools will be used to extract and deduplicate UMIs from reads prior to transcript/gene quantification.
--skip_umi_extract | false | Options: false, true. If specified, UMI extraction step will be skipped. Used in cases when UMIs have already been extracted from reads and placed in read headers.
--umi_separator | ':' | The character that separates the UMI from the read name in FASTQ files. Default is ':' (e.g., @READ_NAME:UMI_SEQUENCE).
--umitools_extract_method | 'regex' | Options: regex, string. Method used by UMI-tools to extract UMIs from reads. 'regex' allows for flexible UMI extraction using regular expressions. 'string' requires a fixed position for UMI extraction. See UMI-tools documentation for details.
--umitools_bc_pattern | null | Pattern used for read 1 UMI extraction. See UMI-tools documentation, and the jds-workflows RNAseq Wiki page for details on how to set this parameter.
--umitools_bc_pattern2 | null | Pattern used for read 2 UMI extraction. See UMI-tools documentation, and the jds-workflows RNAseq Wiki page for details on how to set this parameter.
--umitools_grouping_method | 'unique' | Method used by UMI-tools to group UMIs. Options include 'unique', 'adjacency', 'directional', etc. See UMI-tools documentation for details.

--skip_read_trimming | false | Options: false, true. If specified, skip FastP read trimming.

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.
--trim_poly_g | false | If enabled, polyG trimming is done. For Illumina NextSeq/NovaSeq data, polyG can happen in read tails since G means no signal in the Illumina two-color systems. fastp can detect the polyG in read tails and trim them. 
--trim_poly_x | false | If enabled, polyX trimming is done. If specified with polyG trimming, that is done first then polyX trimming is done. A minimum length can be set with --poly_x_min_len for fastp to detect polyX
--poly_x_min_len | 10 | Minimum length of polyX to be trimmed. Default is 10.

--strandedness_ref | Mouse: '/projects/compsci/omics_share/human/GRCh38/transcriptome/indices/ensembl/v104/kallisto/kallisto_index'
                   | Human: '/projects/compsci/omics_share/human/GRCh38/transcriptome/indices/ensembl/v104/kallisto/kallisto_index' 
                   | Modfied kallisto index file used in strandedness determination. 
--strandedness_gtf | Mouse: '/projects/compsci/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.gtf'
                   | Human: '/projects/compsci/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.gtf' 
                   | GTF file used with kallisto index file used in strandedness determination. 
--strandedness     | null | Library strandedness override. Supported options are 'reverse_stranded' or 'forward_stranded' or 'non_stranded'. This override parameter is only used when the tool `check_strandedness` fails to classify the strandedness of a sample. If the tool provides a strand direction, that determination is used." 

--rsem_ref_files | /projects/omics_share/mouse/GRCm38/transcriptome/indices/ensembl/v102 | Pre-compiled index files. Refers to human indices when --gen_org human.
--rsem_ref_prefix | 'Mus_musculus.GRCm38.dna.toplevel' | Prefix for index files. Refers to human indices when --gen_org human.
--seed_length | 25 | 'Bowtie2 seed length used by the read aligner. Providing the correct value is important for RSEM. If RSEM runs Bowtie2, it uses this value for Bowtie2's seed length parameter.'
--rsem_aligner | 'star' | Options: star or bowtie2. The aligner algorithm used by RSEM.

--fragment_length_mean | 280 | Used when --read_type == 'SE'. "The mean of the fragment length distribution, which is assumed to be a Gaussian."
--fragment_length_sd | 50 | Used when --read_type == 'SE'. "The standard deviation of the fragment length distribution, which is assumed to be a Gaussian. "

--merge_rna_counts | false | Options false, true. If specified, gene and transcript counts are merged across all samples. Typically used in multi-sample cases. 

--picard_dict | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.dict' 
              | Human: '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.toplevel.dict'
              | The coverage metric calculation step requires this file. Refers to human assembly when --gen_org human. 

--ref_flat | Mouse: '/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.chr_patch_hapl_scaff.refFlat.txt' 
           | Human: '/projects/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.refFlat.txt'
           | The coverage metric calculation step requires this file. Refers to human assembly when --gen_org human. 

--ribo_intervals | Mouse: '/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.chr_patch_hapl_scaff.rRNA.interval_list' 
                 | Human: '/projects/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.rRNA.interval_list'
                 | The coverage metric calculation step requires this file. Refers to human assembly when --gen_org human. 

--pdx | false | Options: false, true. If specified, 'Xengsort' is run on reads to deconvolute human and mouse reads. Human only reads are used in analysis. 
--classifier_table | '/projects/compsci/omics_share/human/GRCh38/supporting_files/rna_ebv_classifier/EBVlym_classifier_table_48.txt' | EBV expected gene signatures used in EBV classifier. Only used when '--pdx' is run. 
--ref_fa | '/projects/compsci/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'| Xengsort graft fasta file. Used by Xengsort Index when `--pdx` is run, and xengsort_idx_path is `null` or false.  
--xengsort_host_fasta | '/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/imputed/rel_2112_v8/NOD_ShiLtJ.39.fa' | Xengsort host fasta file. Used by Xengsort Index when `--pdx` is run, and xengsort_idx_path is `null` or false.  
--xengsort_idx_path | '/projects/compsci/omics_share/human/GRCh38/supporting_files/xengsort' | Xengsort index for deconvolution of human and mouse reads. Used when `--pdx` is run. If `null`, Xengsort Index is run using ref_fa and host_fa.  
--xengsort_idx_name | 'hg38_GRCm39-NOD_ShiLtJ' | Xengsort index name associated with files located in `xengsort_idx_path` or name given to outputs produced by Xengsort Index

There are two additional parameters that are human specific. They are: 

Parameter| Default| Description

--ref_fa | '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.toplevel.fa'| Reference fasta to be used in alignment calculation as well as any downstream analysis. 
--ref_fai | '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.toplevel.fa.fai' | Reference fasta index file.  

--rsem_reference_path | null | RSEM reference directory. Used when --bam_input is specified. If not provided --ref_fa / --ref_gtf will be used to generate one. 
--rsem_reference_name | null | User provided path to an RSEM reference name. If not provided, --ref_fa / --ref_gtf will be used to generate one. 
--ref_fa | null | Genomic reference file used to build and RSEM reference when --bam_input is specified, and --rsem_reference_path / --rsem_reference_name are not provided. 
--ref_gtf | null | Transcript GTF reference file used to build and RSEM reference when --bam_input is specified, and --rsem_reference_path / --rsem_reference_name are not provided. 
--bam_strandedness | null | Used when --bam_input is specified. Must be set to 'forward', 'reverse' or 'none'

'''
}




