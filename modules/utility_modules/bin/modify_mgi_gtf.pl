#!/usr/bin/perl
use strict;
use warnings;

# Sai.Lek@jax.org
# 01/21/2026

# Usage message displayed if incorrect arguments provided
my $usage = "Usage: perl script.pl <gtf_file> <reference_file>\n";

# Command-line argument validation
# $ARGV[0] = input GTF file path
# $ARGV[1] = reference file path
# 'or die' exits script if arguments missing
my $gtf_in     = $ARGV[0] or die $usage;
my $ref_file   = $ARGV[1] or die $usage;
my $gtf_out    = $ARGV[2] or die $usage;

# Print confirmation of which files are being processed
print "GTF input: $gtf_in\n";
print "Reference file: $ref_file\n";
print "GTF output: $gtf_out\n";

# Hash data structure to store reference information
# Key: MGI_ID (Mouse Genome Informatics ID)
# Value: Array containing [gene_id, gene_name, gene_biotype]
my %refInfo;

# ===== PHASE 1: PARSE REFERENCE FILE =====

# Open reference file for reading
# '<' indicates read mode
open (INref, "< $ref_file");

# Loop through each line of reference file
# $_ is Perl's default variable (current line)
while (<INref>) {
  # Skip lines that are blank, contain header "ENSG_ID", or start with comments (#)
  next if /^(\s+$|ENSG_ID|\#)/;
  
  # Remove trailing newline from $_ 
  chomp;
  
  # Split line on one or more tab characters (\t+)
  # Store resulting array in @t
  my @t = split(/\t+/, $_);
  
  # Extract specific columns from reference file
  # Column indices: gene_id, gene_name, biotype, mgi_id
  my $id      = $t[0];      # Ensembl gene ID
  my $name    = $t[1];      # Gene name
  my $biotype = $t[2];      # Gene biotype (protein_coding, lncRNA, etc.)
  my $mgi_id  = $t[3];      # MGI ID (used as lookup key)
  
  # Store reference data in hash using MGI ID as key
  # Array notation [0], [1], [2] stores multiple values per key
  $refInfo{$mgi_id}[0] = $id;       # Store gene ID
  $refInfo{$mgi_id}[1] = $name;     # Store gene name
  $refInfo{$mgi_id}[2] = $biotype;  # Store gene biotype
}

# Close reference file
close(INref);

# ===== PHASE 2: PROCESS GTF FILE AND APPEND BIOTYPE =====

# Open GTF input file for reading
open (INgtf, "< $gtf_in");

# Open output file for writing
# '>' indicates write mode (overwrites if file exists)
open (OUTgtf, "> $gtf_out");

# Loop through each line of GTF file
while (<INgtf>) {
  # Skip blank lines, header lines with "ENSG_ID", and comment lines (starting with #)
  next if /^(\s+$|ENSG_ID|\#)/;
  
  # Remove trailing newline from current line
  chomp($_);
  
  # Store the complete original line (we'll append to this later)
  my $line = $_;
  
  # Split GTF line on one or more tab characters
  # GTF format has 9 columns, with attributes in column 9 ($t[8])
  my @t = split(/\t+/, $_);
  
  # Extract the attributes column (9th column, index 8)
  my $gene_id = $t[8]; 
  chomp $gene_id;
  
  # Parse out the MGI ID from the attributes string
  # Example: gene_id "ENSMUSG00000000001"; gene_name "Gnai3"; mgi_id "MGI:95773";
  # This regex removes everything before "MGI:" leaving just "MGI:12345"
  # \S+ matches non-whitespace, repeated 4 times to skip tokens before MGI:
  $gene_id =~ s/\S+\s+\S+\s+\S+\s+\S+MGI:/MGI:/g;
  
  # Remove the trailing ";" and any quotes
  $gene_id =~ s/";//g;
  chomp $gene_id;

  # Clean up the gene_id string by removing "gene_id " prefix and trailing semicolon
  # After these operations, $gene_id should be just "MGI:12345"
  $gene_id =~ s/gene_id //;
  $gene_id =~ s/;//;
 
  # Look up the biotype for this MGI ID from the reference hash
  # $refInfo{$gene_id}[2] retrieves the biotype value stored earlier
  my $biotype = $refInfo{$gene_id}[2];

  # Remove the "biotype=" prefix from the biotype value if present
  $biotype =~ s/biotype=//g;
  
  # Format the biotype as a proper GTF attribute
  # Result: gene_biotype "protein_coding";
  $biotype = "gene_biotype " . "\"$biotype\"" . "\;";

  # Append the formatted biotype to the original GTF line with a tab separator
  $line .= "\t" . $biotype;
  
  # Write the enhanced line to output file
  print OUTgtf "$line\n";

}

# Close output file
close OUTgtf;

# Print completion message
print "GTF file enriched and saved to: $gtf_out\n";

