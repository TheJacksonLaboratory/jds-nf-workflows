def help(){
  println '''
Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.

-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--gen_org | 'mouse' | The organism for the analysis.
--genome_build | 'GRCm38' | The genome build version.
--mt_contig_name | 'MT' | Name of the mitochondrial contig.
--mt_fasta | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/Mus_musculus.GRCm38.dna.MT.fa' | Path to the mitochondrial FASTA file.
--mt_genome | params.reference_cache+'/mouse/GRCm38/genome/sequence/ensembl/GRCm38.p6/Mus_musculus.GRCm38.dna.primary_assembly.genome' | Path to the full genome file.
--mt_shifted_fasta | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/Mus_musculus.GRCm38.dna.MT.shifted8k.fa' | Path to the shifted mitochondrial FASTA file.
--shift_back_chain | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/Mus_musculus.GRCm38.dna.MT.shiftback' | Path to the shift-back chain file.
--mt_fasta_index | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/bwa_index/Mus_musculus.GRCm38.dna.MT' | Path to the BWA index for the mitochondrial FASTA.
--mt_shifted_fasta_index | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/bwa_index/Mus_musculus.GRCm38.dna.MT.shifted8k' | Path to the BWA index for the shifted mitochondrial FASTA.
--max_allele_count | 4 | Maximum number of alleles to consider.
--blacklisted_sites | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/Mus_musculus.GRCm38.dna.MT.exclusion.bed' | BED file of blacklisted sites.
--non_control_region_interval_list | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/non_control_region.chrMT.GRCm38.interval_list' | Interval list for non-control region.
--control_region_shifted_interval_list | '/projects/compsci/omics_share/meta/benchmarking/mt_dna/shifted_ref/GRCm38/control_region_shifted.chrMT.GRCm38.interval_list' | Interval list for shifted control region.
--detection_limit | 0.01 | Detection limit for Mutserve.
--mapQ | 20 | Minimum mapping quality.
--baseQ | 20 | Minimum base quality.
--ref_fa | params.reference_cache+'/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.primary_assembly.fa' | Reference FASTA file (not used, but set for warning).
'''
}
