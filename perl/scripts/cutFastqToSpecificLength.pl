#!/usr/bin/env perl

eval 'exec /jgi/tools/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
    itagger_cut_fastq_seq.pl
PURPOSE:
    To remove a specific number of nucleotides from either
    the beginning or end of reads (or both). Useful if we 
    have a library that consistently shows bad quality after 
    nt 125..., or a library in which the last nucleotides 
    are consistently of bad quality.
    Also useful to trimm reads before assembly with FLASH.
INPUT:
    --infile <string>    : fastq infile
    --final length <int> : Will trim sequences to that length.
    --keep_shorter <int> : If set, will keep reads shorter than <final_length>
OUTPUT:
    --outfile <string>   : trimmed fastq file.
NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
     Julien Tremblay - jtremblay514@gmail.comjtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $final_length, $keep_shorter);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'infile=s' 			=> \$infile,
	'outfile=s'			=> \$outfile,
	'final_length=i' 	=> \$final_length,
	'keep_shorter' 		=> \$keep_shorter,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

#VALIDATE
die("--fastq fastq file required\n") 			unless $infile;
die("--begin or --end int value required\n") 	unless($final_length);
die("--outfile outfile required\n") 			unless $outfile;

#MAIN=====================================================================================================================================================================
open(OUT, ">".$outfile) or die "Can't open file ".$!."\n";

my $in = new Iterator::FastqDb($infile, {trimN=>0}) or die("Unable to open Fastq file, $infile\n");

while( my $curr = $in->next_seq() ){
	my $length = length($curr->seq());
	if($length >= $final_length){
		print OUT $curr->header."\n".substr($curr->seq, 0, $final_length)."\n+\n".substr($curr->qual, 0, $final_length)."\n";
	}elsif($keep_shorter){
		print OUT $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";	
	}
}
close(OUT);

exit;
