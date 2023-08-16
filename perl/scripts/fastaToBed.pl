#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fastaToBed.pl

PURPOSE:

INPUT:
--fasta <string> : Sequence file
				
OUTPUT:
STDOUT            : Bed file.

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
   'fasta=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    #my $header = $curr->header;
    #$header = substr $header, 1;
    my ($header) = $curr->header =~ /^>(\S+)/;
    print STDOUT $header."\t0\t".length($curr->seq)."\n";
}

