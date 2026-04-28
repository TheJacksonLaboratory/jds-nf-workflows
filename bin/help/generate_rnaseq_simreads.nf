def help(){
  println '''

Note: This workflow generates simulated RNA reads from a set of transcripts at specified abundances. Each transcript has a predefined abundance in the output.

Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--gen_org | mouse | Options: mouse and human.

--fa_cds | Mouse: 'http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz' 
         | Human: 'http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz'
         | All transcript coding sequences resulting from Ensembl genes.

--library_size | '1000000, 3000000' | Number of simulated reads. Comma delimited values will generate data for each size.
--library_strategy | 'PE' | Library strategy: PE (paired-end), SE (single-end), BOTH (both paired-end and single-end)
--read_length | 100 | A anumeric value giving the length of each read in the output. Maximum length is 250bp. Transcripts that are shorter than the read length will not be used for generating simulated reads.
--fragment_length_min | 125 | A numeric value giving the minimum fragment length. The minimum fragment length must be equal to or greater than the output read length.
--fragment_length_max | 500 | A numeric value giving the maximum fragment length.
--fragment_length_mean | 250 | A numeric value giving the mean fragment length.
--fragment_length_sd | 40 | A numeric value giving the standard deviation of the fragment length.
--simulate_sequencing_error | false | Logical indicating if sequencing errors should be simulated in the output reads. If TRUE, the quality.reference parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, the quality.reference parameter can be optionally omitted, and the function will use its inbuilt quality strings.
--quality_reference | '' | A character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred+33 format. The sequencing quality strings must have the same length as read.length.
'''
}
