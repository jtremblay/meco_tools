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
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %seen;
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   if(exists $seen{$curr->header}){

   }else{
      print STDOUT $curr->header."\n".$curr->seq."\n";
      $seen{$curr->header}=1
   }
}

