def help(){
  println '''

Note: This workflow relies on 'biotype' annotations within the GTF to define 'rRNA' regions. These regions are used in quality control steps to check percentage of reads mapping to rRNA. `--mgi` is provided to convert MGI based GFF3 files to GTF with the appropriate biotype annotations for rRNA. If using a custom GTF, please ensure that rRNA regions are annotated with biotype 'rRNA' for accurate QC metrics.

Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--ref_fa | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' or GRCm39: '/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa' or MGI: '/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v114/Mus_musculus.GRCm39.dna.primary_assembly.fa'
         | Human: '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
         | The reference fasta to be used in index generation.

--ref_gtf | Mouse: '/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.gtf' or GRCm39: '/projects/omics_share/mouse/GRCm39/transcriptome/annotation/ensembl/v105/Mus_musculus.GRCm39.105.gtf' or MGI: '/projects/omics_share/mouse/GRCm39/transcriptome/annotation/mgi/v114/MGI.gtf'
          | Human: '/projects/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.gtf'
          | The reference gtf to be used in index generation.

--ref_gff | null | Used in cases where no GTF is available. When specifiying GFF, set `--ref_gtf FALSE` in the Nextflow command when using this param. The workflow first converts to GTF using AGAT_GFFTOGTF before continuing.

--ref_gff3 | null | Used in cases where no GTF is available. When specifiying GFF3, set `--ref_gtf FALSE` in the Nextflow command when using this param. The workflow first converts to GTF using GFFREAD_GFF3TOGTF and MODIFY_MGI_GTF before continuing.

--mgi | false | If specified, the MGI GTF (either provided or converted) is appended to include 'biotype' annotations.

--annotation_source | ensembl | Source of transcriptome annotation. Used for internal run tracking or setting default parameters to "MGI"

--custom_gene_fasta | null | The path to a fasta file with additonal transcript sequences to add to the index. Will be annotated based on the name provided in the sequnece name field. For example: ">New_Gene_42", where New_Gene_42 will be the name of the gene, transcript, and exon. 

'''
}
