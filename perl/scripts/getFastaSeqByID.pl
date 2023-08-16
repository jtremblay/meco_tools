#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
--term <string>   : Term to match.

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $term);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'term=s' 	=> \$term,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fastq file, $infile\n");
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   my($parsed_header) = $header =~ m/(^\S+)/;
   if($parsed_header eq $term){
      print STDOUT $curr->output;
   }

}

