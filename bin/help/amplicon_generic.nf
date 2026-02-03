def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | Directory where saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | Directory containing cached Singularity containers.
-w | /<PATH> | Working directory for Nextflow (intermediate files and work directories).
--sample_folder | /<PATH> | Path to folder containing sample files (can include symlinks).
--extension | .fastq.gz | Expected read file extension.
--pattern | '*_R{1,2}*' | Pattern to match R1/R2 filenames (e.g. READ_R1.fastq.gz).
--read_type | PE | Read type: `PE` (paired-end) or `SE` (single-end).
--concat_lanes | false | If `true`, concatenate FASTQs across lanes per sample.
--csv_input | null | CSV manifest with header: `sampleID,lane,fastq_1,fastq_2` (fastq_2 optional for PE).
--download_data | null | When using URLs in `--csv_input`, download remote files.
--multiqc_config | /<PATH> | Path to the MultiQC configuration (amplicon.yaml).

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.
--trim_poly_g | false | If enabled, polyG trimming is done. For Illumina NextSeq/NovaSeq data, polyG can happen in read tails since G means no signal in the Illumina two-color systems. fastp can detect the polyG in read tails and trim them. 
--trim_poly_x | false | If enabled, polyX trimming is done. If specified with polyG trimming, that is done first then polyX trimming is done. A minimum length can be set with --poly_x_min_len for fastp to detect polyX
--poly_x_min_len | 10 | Minimum length of polyX to be trimmed. Default is 10.

--ref_fa | /projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta | Reference fasta used for alignment and downstream analysis.
--ref_fa_indices | /projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta | Pre-built BWA index prefix.
--mismatch_penalty | 8 | BWA mismatch penalty value.

--markduplicates | false | This option, if specified, will mark duplicates reads following mapping, while this is possible, it is NOT recommended for amplicon data.

--amplicon_primer_intervals | <path> | GATK interval_list of primer intervals for coverage metrics.
--amplicon_target_intervals | <path> | GATK interval_list of target intervals for coverage metrics.

--gold_std_indels | /projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | (Human only) Indels VCF for GATK BaseRecalibrator.
--phase1_1000G | /projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz | (Human only) 1000G SNP set for BaseRecalibrator.
--dbSNP | /projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz | dbSNP VCF for known variant annotation.
--dbSNP_index | /projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz.tbi | Index for `--dbSNP` VCF.
--call_val | 50 | Minimum phred-scaled confidence threshold for variant calling.
--ploidy_val | '2' | Ploidy value passed to HaplotypeCaller (default 2).
--target_gatk | <path> | BED file of amplicon targets used by GATK (must match assay design).
--bwa_min_score | (unset) | Minimum BWA alignment score filter (if used by downstream steps).

'''
}
