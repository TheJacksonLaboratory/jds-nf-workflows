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
--merge_inds | false | In some use cases, samples are structured by a higher organizational level. If specified, `merge_ind` merges of BAMs to the ind level prior to calling (e.g., Ind_42 <-- sampleA, sampleB, sampleC).

--run_sv | false | Options: false and true. Default: false. If this boolean is specified, structural variant calling will be performed.
--run_mt_calling | false | Options: false and true. Default: false. If this boolean is specified, mitochondrial variant calling will be performed.

--deduplicate_reads | false | Options: false, true. If specified, run bbmap clumpify on input reads. Clumpify will deduplicate reads prior to trimming. This can help with mapping and downstream steps when analyzing high coverage WGS data.

--split_fastq | false | Options false, true. If specified, FASTQ files will be split into chunks sized based on split_fastq_bin_size prior to mapping. This option is recommended for high coverage data. 
--split_fastq_bin_size | 10000000 | If split_fastq is specified, FASTQ files will splint into chunks of this size prior to mapping. 

--coverage_cap | null | If an integer value is specified, jvarkit 'Biostar154220' is used to cap coverage at the that value. See: http://lindenb.github.io/jvarkit/Biostar154220.html
--primary_chrom_bed | '/projects/compsci/omics_share/mouse/GRCm38/genome/annotation/intervals/Mus_musculus.GRCm38.dna.primary_assembly.bed' | A bed file containing the primary chromsomes with positions. Used in limiting jvarkit 'Biostar154220' to those regions with expected coverage.

--run_gvcf | false | Options: false and true. Default: false. If this boolean is specified, GCVF output will be generated.

--gen_org | mouse | Options: mouse, human, other.
--genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCh38. If other, this parameter is not used.

--ref_fa | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
         | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'
         | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. 

--ref_fa_indices | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bwa/Mus_musculus.GRCm38.dna.toplevel.fa'
                 | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta'
                 | Pre-compiled BWA index files, points to human reference when --gen_org human. 

--chrom_contigs | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.primaryChr.contig_list' 
                | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.primaryChr.contig_list'
                | A list of all chromosomes, unplaced, and unlocalized contigs present in the reference file, points to human reference when --gen_org human. Used to scatter variant calling by chromosome. 

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.
--trim_poly_g | false | If enabled, polyG trimming is done. For Illumina NextSeq/NovaSeq data, polyG can happen in read tails since G means no signal in the Illumina two-color systems. fastp can detect the polyG in read tails and trim them. 
--trim_poly_x | false | If enabled, polyX trimming is done. If specified with polyG trimming, that is done first then polyX trimming is done. A minimum length can be set with --poly_x_min_len for fastp to detect polyX
--poly_x_min_len | 10 | Minimum length of polyX to be trimmed. Default is 10.

--deepvariant | false | Options: false and true. Default: false. If this boolean is specified, Google DeepVariant will be used for variant calling rather than GATK HaplotypeCaller. This option requires csv_input with `sex` as a provided column.

--mismatch_penalty | 8 | The BWA penalty for a mismatch.
--call_val | 50 | The minimum phred-scaled confidence threshold at which variants should be called.
--ploidy_val | '2' | Sample ploidy

--dbSNP | Mouse: '/projects/omics_share/mouse/GRCm38/genome/annotation/snps_indels/GCA_000001635.6_current_ids.vcf.gz' 
        | Human: '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz'
        | The dbSNP database contains known single nucleotide polymorphisms, and is used in the annotation of known variants. Points to human dbSNP when --gen_org human.

--gen_ver | Mouse: 'GRCm38.99'
          | Human: 'hg38'
          | snpEff genome version. Sets to 'hg38' when --gen_org human

--snpEff_config | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/snpEff_5_1/snpEff.config' 
                | Human: '/projects/omics_share/human/GRCh38/genome/indices/snpEff_5_1/snpEff.config'
                | The configuration file used while running snpEff, points to human snpEff file when --gen_org human. 

--gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Human Only - Used in GATK BaseRecalibrator. 
--phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Human Only - Used in GATK BaseRecalibrator. 
--dbNSFP | '/projects/omics_share/human/GRCh38/genome/annotation/function/dbNSFP4.2a.gatk_formatted.txt.gz' | Human Only - Used in variant annotation.
--cosmic | '/projects/omics_share/human/GRCh38/genome/annotation/function/COSMICv95_Coding_Noncoding.gatk_formatted.vcf' | Human Only - Used in variant annotation.
'''

if (params.gen_org == 'mouse' && params.run_sv)
  println '''
  --min_sv_length | <INT> | Minimum length of SVs to report.
  --cnv_distance_limit | <INT> | Maximum distance allowed between SV and nearest CNV.
  --sv_slop | <INT> | Number of bases to extend SV breakpoints to merge.
  --sizemargin | 0.8 | Error margin in allowable size to prevent matching of SVs of different sizes.
  --smoove_support | 3 | Minimum number of supporting reads for SV calling with smoove.
  --exclude_regions | /<PATH> | BED file of regions to exclude from SV calling.
  --callRegions | /<PATH> | BED file of regions to call SVs, used with MANTA.
  --ref_fa_dict | /<PATH> | Reference fasta dictionary file for SV calling.
  --delly_exclusion | /<PATH> | TSV file of regions to exclude for Delly SV calling.
  --delly_mappability | /<PATH> | Mappability file for Delly SV calling.
  --cnv_window | 10000 | Window size for CNV calling.
  --cnv_min_size | 10000 | Minimum CNV size to report.
  --combined_reference_set | /<PATH> | Reference fasta for SVABA (must be in same directory as BWA index).
  --cytoband | /<PATH> | Cytoband annotation file for CNV annotation.
  --known_del | /<PATH> | BED file of known deletions for annotation.
  --known_ins | /<PATH> | BED file of known insertions for annotation.
  --known_inv | /<PATH> | BED file of known inversions for annotation.
  --ensemblUniqueBed | /<PATH> | BED file of unique Ensembl genes for annotation.
  --gap | /<PATH> | BED file with gap genomic regions for annotation.  
  '''

if (params.gen_org == 'human' && params.run_sv)
  println '''
  --min_sv_length | <INT> | Minimum length of SVs to report.
  --cnv_distance_limit | <INT> | Maximum distance allowed between SV and nearest CNV.
  --sv_slop | <INT> | Number of bases to extend SV breakpoints to merge.
  --sizemargin | 0.8 | Error margin in allowable size to prevent matching of SVs of different sizes.
  --smoove_support | 3 | Minimum number of supporting reads for SV calling with smoove.
  --exclude_regions | /<PATH> | BED file of regions to exclude from SV calling.
  --ref_fa_dict | /<PATH> | Reference fasta dictionary file for SV calling.
  --delly_exclusion | /<PATH> | TSV file of regions to exclude for Delly SV calling.
  --delly_mappability | /<PATH> | Mappability file for Delly SV calling.
  --cnv_window | 10000 | Window size for CNV calling.
  --cnv_min_size | 10000 | Minimum CNV size to report.
  --callRegions | /<PATH> | BED file of regions to call SVs, used with MANTA.
  --combined_reference_set | /<PATH> | Reference fasta for SVABA (must be in same directory as BWA index).
  --cytoband | /<PATH> | Cytoband annotation file for CNV annotation.
  --dgv | /<PATH> | DGV annotation file for CNV annotation.
  --thousandG | /<PATH> | 1000 Genomes CNV annotation file.
  --cosmicUniqueBed | /<PATH> | COSMIC unique intervals for CNV annotation.
  --ensemblUniqueBed | /<PATH> | Ensembl unique genes for CNV and SV annotation.
  --gap | /<PATH> | BED file with gap genomic regions for annotation.  
  --dgvBedpe | /<PATH> | Database of Genomic Variants BEDPE file for SV annotation.
  --thousandGVcf | /<PATH> | 1000 Genomes SV VCF for SV annotation.
  --svPon | /<PATH> | SV Panel of Normals BEDPE for SV annotation.
  --cosmicBedPe | /<PATH> | COSMIC SV BEDPE for SV annotation.
  '''

  if (params.run_mt_calling)
    println '''
  --mt_contig_name | 'MT' | Name of the mitochondrial contig.
  --mt_fasta | <PATH> | Path to the mitochondrial fasta file.
  --mt_genome | <PATH> | Path to the mitochondrial genome file.
  --mt_shifted_fasta | <PATH> | Path to the shifted mitochondrial fasta file.
  --shift_back_chain | <PATH> | Path to the shift back chain file.
  --mt_fasta_index | <PATH> | Path to the mitochondrial fasta index.
  --mt_shifted_fasta_index | <PATH> | Path to the shifted mitochondrial fasta index.
  --max_allele_count | 4 | Maximum allele count for mitochondrial variant calling.
  --exclusion_sites | <PATH> | BED file of exclusion sites for mitochondrial calling.
  --non_control_region_interval_list | <PATH> | Interval list for non-control region of chrMT.
  --control_region_shifted_interval_list | <PATH> | Interval list for shifted control region of chrMT.
  --detection_limit | 0.01 | Mutserve detection limit.
  --mapQ | 20 | Minimum mapping quality for Mutserve.
  --baseQ | 20 | Minimum base quality for Mutserve.
  '''
  
}

