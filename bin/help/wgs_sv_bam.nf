def help(){
    
  println '''
  Parameter | Default | Description

  --pubdir | /<PATH> | The directory that the saved outputs will be stored.
  --cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
  -w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

  --csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.

  --gen_org | mouse | Options: mouse, human, other.
  --genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCh38. If other, this parameter is not used.

  --ref_fa | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
          | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'
          | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. JAX users should not change this parameter.

  --ref_fa_indices | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bwa/Mus_musculus.GRCm38.dna.toplevel.fa'
                  | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta'
                  | Pre-compiled BWA index files, points to human reference when --gen_org human. JAX users should not change this parameter.

  --dbSNP | Mouse: '/projects/omics_share/mouse/GRCm38/genome/annotation/snps_indels/GCA_000001635.6_current_ids.vcf.gz' 
          | Human: '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz'
          | The dbSNP database contains known single nucleotide polymorphisms, and is used in the annotation of known variants. Points to human dbSNP when --gen_org human.

  --gen_ver | Mouse: 'GRCm38.99'
            | Human: 'hg38'
            | snpEff genome version. Sets to 'hg38' when --gen_org human

  --snpEff_config | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/snpEff_5_1/snpEff.config' 
                  | Human: '/projects/omics_share/human/GRCh38/genome/indices/snpEff_5_1/snpEff.config'
                  | The configuration file used while running snpEff, points to human snpEff file when --gen_org human. JAX users should not change this parameter.
  '''

  if (params.gen_org == 'mouse')
    println '''
    --min_sv_length | <INT> | Minimum length of SVs to report.
    --cnv_distance_limit | <INT> | Maximum distance allowed between SV and nearest CNV.
    --sv_slop | <INT> | Number of bases to extend SV breakpoints to merge.
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

  if (params.gen_org == 'human')
    println '''
    --min_sv_length | <INT> | Minimum length of SVs to report.
    --cnv_distance_limit | <INT> | Maximum distance allowed between SV and nearest CNV.
    --sv_slop | <INT> | Number of bases to extend SV breakpoints to merge.

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

    --gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Used in GATK BaseRecalibrator. JAX users should not change this parameter.
    --phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Used in GATK BaseRecalibrator. JAX users should not change this parameter.
    --dbNSFP | '/projects/omics_share/human/GRCh38/genome/annotation/function/dbNSFP4.2a.gatk_formatted.txt.gz' | Used in variant annotation.
    --cosmic | '/projects/omics_share/human/GRCh38/genome/annotation/function/COSMICv95_Coding_Noncoding.gatk_formatted.vcf' | Used in variant annotation.
    '''
}
