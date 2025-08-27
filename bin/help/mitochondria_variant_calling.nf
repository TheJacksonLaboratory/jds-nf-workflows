def help(){
  println '''
Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.

-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--gen_org | 'mouse' | The organism for the analysis.
--genome_build | 'GRCm38' | The genome build version.
--mt_contig_name | 'MT' | Name of the mitochondrial contig.
--mt_fasta | </fasta> | Path to the mitochondrial FASTA file.
--mt_genome | </.genome> | Path to the full genome file.
--mt_shifted_fasta | </shifted_fasta> | Path to the shifted mitochondrial FASTA file.
--shift_back_chain | </.shiftback> | Path to the shift-back chain file.
--mt_fasta_index | </bwa_index> | Path to the BWA index for the mitochondrial FASTA.
--mt_shifted_fasta_index | </shifted_bwa_index> | Path to the BWA index for the shifted mitochondrial FASTA.
--max_allele_count | 4 | Maximum number of alleles to consider.
--exclusion_sites | </exclusion_sites> | BED file of exclusion sites.
--non_control_region_interval_list | </non_control_region.interval_list> | Interval list for non-control region.
--control_region_shifted_interval_list | </control_region_shifted.interval_list> | Interval list for shifted control region.
--detection_limit | 0.01 | Detection limit for Mutserve.
--mapQ | 20 | Minimum mapping quality.
--baseQ | 20 | Minimum base quality.
--gen_ver | 'hg38' | Genome version for the analysis.
--dbSNP | </dbsnp.vcf.gz> | Path to dbSNP annotation VCF file.
--dbSNP_index | </dbsnp.vcf.gz.tbi> | Path to dbSNP annotation VCF index file.
--snpEff_config | </snpEff.config> | Path to snpEff configuration file.
--cosmic | </cosmic> | Path to COSMIC annotation VCF file. Used when `--gen_org == human`
--cosmic_index | </cosmic_index> | Path to COSMIC annotation VCF index file. Used when `--gen_org == human`
--dbNSFP | </dbNSFP> | Path to dbNSFP annotation file. Used when `--gen_org == human`
'''
}
