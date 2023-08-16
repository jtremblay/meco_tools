#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
getGenesFaaFromGeneList.pl

PURPOSE:

INPUT:
--infile_faa <string>       : Sequence file.
--infile_gene_list <string> : File with list of genes.
            
OUTPUT:
STDOUT                      : Fasta (faa) file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_faa, $infile_gene_list);
my $verbose = 0;

GetOptions(
   'infile_faa=s' 	   => \$infile_faa,
   'infile_gene_list=s' => \$infile_gene_list,
   'verbose' 	         => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

my %hash;
open(IN, "<".$infile_gene_list) or die "Can't open $infile_gene_list\n";
while(<IN>){
   chomp;
   #my @row = split(/\t/, $_);
   my $gene_id = $_;
   $hash{$gene_id} = "";
}
close(IN);

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile_faa) or die("Unable to open Fasta file, $infile_faa\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   
   if(exists $hash{$header}){
     print STDOUT $curr->output; 
   }
}


