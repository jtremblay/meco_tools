#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
splitNcbiFaaByGenes.pl

PURPOSE:

INPUT:
--infile_faa <string>    : Sequence file (all faas)
--infile_blast <string>  : Blast tsv file of a given gene.
				
OUTPUT:
STDOUT                   : Sequence file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_faa, $infile_blast);
my $verbose = 0;

GetOptions(
   'infile_faa=s'   => \$infile_faa,
   'infile_blast=s' => \$infile_blast,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($infile_faa) or die("Unable to open Fasta file, $infile_faa\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my($header) = $curr->header =~ />(gi\|\d+\|)/;
   $hash{$header}{SEQ} = $curr->seq;
   $hash{$header}{FULLNAME} = $curr->header;
}

open(IN, "<".$infile_blast) or die "Can't open $infile_blast\n";
while(<IN>){
   chomp;
   my($header) = $_ =~ /(gi\|\d+\|)/;
   if(exists $hash{$header}){
      print STDOUT $hash{$header}{FULLNAME}."\n".$hash{$header}{SEQ}."\n";
   }

}
