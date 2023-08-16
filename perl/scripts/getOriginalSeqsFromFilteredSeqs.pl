#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
getOriginalSeqsFromFilteredSeqs.pl

PURPOSE:
From a filtered/trimmed fastq file, store all headers in 
a hash, loop through original fastq file and write
original output sequences in a new file.

INPUT:
--infile <string>         : Sequence file fastq
--originalInfile <string> : Sequence file fastq
--outfile <string         : Sequence file fastq
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $originalInfile);
my $verbose = 0;

GetOptions(
    'infile=s' 	       => \$infile,
	'originalInfile=s' => \$originalInfile,
	'outfile=s'        => \$outfile,
    'verbose' 	       => \$verbose,
    'help' 		       => \$help
);
if ($help) { print $usage; exit; }

die "--infile missing\n" unless($infile);
die "--originalInfile missing\n" unless($originalInfile);
die "--outfile missing\n" unless($outfile);

## MAIN

open(OUT, '>'.$outfile) or die "Can't open file $outfile\n";

my %hash;

my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fastq file, $infile\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
	$hash{$curr->header}++;
}

$ref_fastq_db = Iterator::FastqDb->new($originalInfile) or die("Unable to open Fastq file, $originalInfile\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
	if(exists $hash{$curr->header}){
		print OUT $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
	}
}
close(OUT);
exit;


