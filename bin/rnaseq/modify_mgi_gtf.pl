#!/usr/bin/perl

use strict;

# Sai.Lek@jax.org
# 01/21/2026,     


my $usage = "
Usage : 
$0  <gtf_in> <ref_file> 
";


my $gtf_in     = $ARGV[0] or die $usage;
my $ref_file   = $ARGV[1] or die $usage;


# Chrom Source  Type	        Start   End     Score   Strand  Phase   Transcript id					  gene_id                biotype 
# 1	ENSEMBL	transcript	3143476	3144545	.	+	.	transcript_id "MGI_C57BL6J_1918292_transcript_1"; gene_id "MGI:1918292"; gene_biotype "processed_pseudogene"; 

my $gtf_out  = 'MGI.gtf';

my %refInfo;

# Read ref file
open (INref, "< $ref_file");
while (<INref>) {
  next if /^(\s+$|ENSG_ID|\#)/;
  chomp($_);
  my $line = $_;
  my @t = split(/\t+/, $_);

  my $id      = $t[0]; chomp $id;  
  my $name    = $t[1]; chomp $name; 
  my $biotype = $t[2]; chomp $biotype;   
  my $mgi_id  = $t[3]; chomp $mgi_id;  

  $refInfo{$mgi_id}[0] = $id;
  $refInfo{$mgi_id}[1] = $name;
  $refInfo{$mgi_id}[2] = $biotype;
}

# To write to modified MGI.gtf
open (OUTgtf, ">$gtf_out");

# Read MGI.gtf.tmp ( output from gff3read utility )
open (INgtf, "< $gtf_in");
while (<INgtf>) {
  next if /^(\s+$|ENSG_ID|\#)/;
  chomp($_);
  my $line = $_;
  my @t = split(/\t+/, $_);
  my $gene_id = $t[8]; chomp $gene_id;
  $gene_id =~ s/\S+\s+\S+\s+\S+\s+\S+MGI:/MGI:/g;
  $gene_id =~ s/";//g;
  chomp $gene_id;

  $gene_id =~ s/gene_id //;
  $gene_id =~ s/;//;
 
  my $biotype = $refInfo{$gene_id}[2];

  $biotype =~ s/biotype=//g;
  $biotype = "gene_biotype " . "\"$biotype\"" . "\;";

  $line .= "\t" . $biotype;
  print OUTgtf "$line\n";

}

close OUTgtf;

