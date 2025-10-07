def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 
--merge_inds | false | In some use cases, samples are structured by a higher organizational level. If specified, `merge_ind` merges of BAMs to the ind level prior to calling (e.g., Ind_42 <-- sampleA, sampleB, sampleC).

--data_type | 'pacbio' | Sequencing data type. Future options may include 'nanopore' or 'illumina'. Default: 'pacbio'. This parameter is used to determine the appropriate mapping and variant calling tools parameter values.

--run_gvcf | false | Options: false and true. Default: false. If this boolean is specified, GCVF output will be generated.

--gen_org | mouse | Options: mouse, human, other.
--genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCh38. If other, this parameter is not used.

--ref_fa | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
         | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'
         | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. JAX users should not change this parameter.

--chrom_contigs | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.primaryChr.contig_list' 
                | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.primaryChr.contig_list'
                | A list of all chromosomes, unplaced, and unlocalized contigs present in the reference file, points to human reference when --gen_org human. Used to scatter variant calling by chromosome. JAX users should not change this parameter.

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.

--deepvariant | true | Workflow is currently only supports deepvariant. This must be 'true' due to logic in the extract_csv function. 
--deepvariant_model_type        ${params.deepvariant_model_type}

--minimap2_index | 	/<PATH> | Minimap index for reference genome, used in mapping. 

--pbmode | CCS | Options: CCS or CLR. Specify whether input data are from PacBio CCS or CLR data.

--dbSNP | Mouse: '/projects/omics_share/mouse/GRCm38/genome/annotation/snps_indels/GCA_000001635.6_current_ids.vcf.gz' 
        | Human: '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz'
        | The dbSNP database contains known single nucleotide polymorphisms, and is used in the annotation of known variants. Points to human dbSNP when --gen_org human.

--dbSNP_index | </PATH> | The dbSNP index file associated with the dbSNP VCF file. 

--gen_ver | Mouse: 'GRCm38.99'
          | Human: 'hg38'
          | snpEff genome version. Sets to 'hg38' when --gen_org human

--snpEff_config | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/snpEff_5_1/snpEff.config' 
                | Human: '/projects/omics_share/human/GRCh38/genome/indices/snpEff_5_1/snpEff.config'
                | The configuration file used while running snpEff, points to human snpEff file when --gen_org human. JAX users should not change this parameter.

--gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--dbNSFP | '/projects/omics_share/human/GRCh38/genome/annotation/function/dbNSFP4.2a.gatk_formatted.txt.gz' | Human Only - Used in variant annotation.
--cosmic | '/projects/omics_share/human/GRCh38/genome/annotation/function/COSMICv95_Coding_Noncoding.gatk_formatted.vcf' | Human Only - Used in variant annotation.

--tandem_repeats | </PATH> | BED file that lists the coordinates of centromeres and telomeres to exclude as alignment targets. 
--pbsv_tandem | null | Optional tandem repeat annotation .bed file of your reference, used by PBSV in structural variant calling.

--min_sv_length | <INT> | Minimum length of SVs to report.
--sv_slop | <INT> | Number of bases to extend SV breakpoints to merge.
--sizemargin | 0.8 | Error margin in allowable size to prevent matching of SVs of different sizes.






--known_del | /<PATH> | BED file of known deletions for annotation.

--known_ins | /<PATH> | BED file of known insertions for annotation.

--known_inv | /<PATH> | BED file of known inversions for annotation.

--ensemblUniqueBed | /<PATH> | BED file of unique Ensembl genes for annotation.

--gap | /<PATH> | BED file with gap genomic regions for annotation.  
'''
}
