#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
linkContigToGene.pl

PURPOSE:

INPUT:
--infile_gff <string>    : gff file (i.e. output of prodigal)
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_gff);
my $verbose = 0;

GetOptions(
   'infile_gff=s' => \$infile_gff,
   'verbose' 	  => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile_gff is missing\n" unless($infile_gff);

## MAIN
# First, rename contigs (remove spurious characters in headers). For instance convert this: ">contig-0 34353 nucleotides" to this: ">contig-0"
my %hash_gene_to_contig;

# Parse Prodigal's gff output.
open(IN, "<".$infile_gff) or die "Can't open $infile_gff\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   #$contig_id = $1 if $contig_id =~ m/^contig-(\d+)/;

   #1 to 7
   my $source     = $row[1];
   my $type       = $row[2];
   my $start      = $row[3];
   my $end        = $row[4];
   my $score      = $row[5];
   my $strand     = $row[6];
   my $phase      = $row[7];
   my $attributes = $row[8];

   my @attributes = split(/;/, $attributes);
   my $gene_id = $attributes[0];
   #$gene_id =~ s/gene_id=//;
   #$gene_id =~ s/gene_id_//;
 
   #$hash_gene_to_contig{$gene_id} = $contig_id;
   print STDOUT $gene_id."\t".$contig_id."\n";
   
}
close(IN);

#print STDOUT "#gene_id\tcontig_id\n";
#for my $key (keys %hash_gene_to_contig){
#   print STDOUT $key."\t".$hash_gene_to_contig{$key}."\n";
#}
exit;
