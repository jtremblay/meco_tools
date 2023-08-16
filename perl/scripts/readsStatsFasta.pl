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
readsStatsFasta.pl

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
    'infile=s' => \$infile,
    'verbose'  => \$verbose,
    'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $stat = new Statistics::Descriptive::Full->new();
my $pipe = IO::Pipe->new();
my $numberOfReads = 0;
my $numberOfBases = 0;
my @readsLength   = ();

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  my $seq = $curr->seq;
  next if(length($seq) < 500);
  push(@readsLength, length($seq));
  $stat->add_data(length($seq));
  $numberOfReads++;
  $numberOfBases = length($seq) + $numberOfBases;
}
my $min = min(@readsLength);
my $max = max(@readsLength);
my $stdDev = $stat->standard_deviation();
my $average = $stat->mean();

print STDOUT "********Stats for fasta file********:\n";
print STDOUT "Number of reads      : $numberOfReads\n"; 
print STDOUT "Number of bases      : $numberOfBases\n"; 
print STDOUT "Std Dev              : $stdDev\n";
print STDOUT "average reads length : $average\n";
print STDOUT "max reads length     : $max\n";
print STDOUT "min reads length     : $min\n";


exit;
