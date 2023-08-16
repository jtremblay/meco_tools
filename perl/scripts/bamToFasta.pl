#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::Pipe;
use List::Util qw(min max);
use Statistics::Descriptive;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file. fastq, fasta or bam
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
  'infile=s' 	=> \$infile,
  'verbose' 	=> \$verbose,
  'help' 		  => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $pipe = IO::Pipe->new();
$pipe->reader("samtools view " . $infile);
while(<$pipe>){
  my @row  = split(/\t/, $_);
  my $seq  = $row[9];
  #my $qual = $row[4];
  print STDOUT ">".$row[0]."#ACGTACGT\n".$seq."\n";
}

#my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
#while( my $curr = $ref_fasta_db->next_seq() ) {
  	

#$counter++;
#}


exit;
