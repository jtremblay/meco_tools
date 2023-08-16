#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fastqToChrSize.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@lbl.gov

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
    'infile=s' 	=> \$infile,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header();
	$header =~ s/>//;
	print STDOUT $header."\t".length($curr->seq())."\n";
}
exit;
