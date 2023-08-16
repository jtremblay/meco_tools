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
--infile <string>  : Sequence file
--length <int>     : length

OUTPUT:
STDOUT             : Met the min length req
--failed           : Did not met the length req

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $failed, $length);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'length=i'  => \$length,
   'failed=s'  => \$failed,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
      print STDOUT $curr->header."\n".substr($curr->seq, $length)."\n+\n".substr($curr->qual, $length)."\n";
}

