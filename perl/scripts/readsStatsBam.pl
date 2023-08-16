#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::Pipe;
use List::Util qw(min max);
use Statistics::Descriptive;

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
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $stat = new Statistics::Descriptive::Full->new();
my $stat1 = new Statistics::Descriptive::Full->new();
my $stat2 = new Statistics::Descriptive::Full->new();
my $stat3 = new Statistics::Descriptive::Full->new();
my $stat4 = new Statistics::Descriptive::Full->new();
my $stat5 = new Statistics::Descriptive::Full->new();
my $stat10 = new Statistics::Descriptive::Full->new();
my $stat20 = new Statistics::Descriptive::Full->new();
my $stat30 = new Statistics::Descriptive::Full->new();
my $stat40 = new Statistics::Descriptive::Full->new();
my $stat50 = new Statistics::Descriptive::Full->new();
my $stat60 = new Statistics::Descriptive::Full->new();
my $pipe = IO::Pipe->new();
my $numberOfReads = 0;
my @readsLength;
my @readsLength1;
my @readsLength2;
my @readsLength3;
my @readsLength4;
my @readsLength5;
#my $qual0 = 0;
my $qual1 = 0;
my $qual2 = 0;
my $qual3 = 0;
my $qual4 = 0;
my $qual5 = 0;
my $qual10 = 0;
my $qual20 = 0;
my $qual30 = 0;
my $qual40 = 0;
my $qual50 = 0;
my $qual60 = 0;

#if($infile =~ m/[\.fasta\.gz|\.fastq\.gz]/){
  $pipe->reader("samtools view " . $infile);
  
    while(<$pipe>){
      my @row = split(/\t/, $_);
      my $seq = $row[9];
      my $qual = $row[4];

      #$qual0++ if($qual >= 0);
      if($qual == 1)  {$qual1++;  $stat1->add_data(length($seq))}
      if($qual == 2)  {$qual2++;  $stat2->add_data(length($seq))}
      if($qual == 3)  {$qual3++;  $stat3->add_data(length($seq))}
      if($qual == 4)  {$qual4++;  $stat4->add_data(length($seq))}
      if($qual == 5)  {$qual5++;  $stat5->add_data(length($seq))}
      if($qual == 10) {$qual10++; $stat10->add_data(length($seq))}
      if($qual >= 20) {$qual20++; $stat20->add_data(length($seq))}
      if($qual >= 30) {$qual30++; $stat30->add_data(length($seq))}
      if($qual >= 40) {$qual40++; $stat40->add_data(length($seq))}
      if($qual >= 50) {$qual50++; $stat50->add_data(length($seq))}
      if($qual >= 60) {$qual60++; $stat60->add_data(length($seq))}

      push(@readsLength, length($seq));
      $stat->add_data(length($seq));
      $numberOfReads++;
    }
  
#}

my $min = min(@readsLength);
my $max = max(@readsLength);
my $stdDev = $stat->standard_deviation();
my $average = $stat->mean();

print STDOUT "********Stats for bamfile********:\n";
print STDOUT "********$infile********:\n";
print STDOUT "Number of reads      : $numberOfReads\n"; 
print STDOUT "Std Dev              : $stdDev\n";
print STDOUT "average reads length : $average\n";
print STDOUT "max reads length     : $max\n";
print STDOUT "min reads length     : $min\n";
print STDOUT "MAPQ == 1            : $qual1\t".$stat1->mean()." +/- ".$stat1->standard_deviation()."\n";
print STDOUT "MAPQ == 2            : $qual2\t".$stat2->mean()." +/- ".$stat2->standard_deviation()."\n";
print STDOUT "MAPQ == 3            : $qual3\t".$stat3->mean()." +/- ".$stat3->standard_deviation()."\n";
print STDOUT "MAPQ == 4            : $qual4\t".$stat4->mean()." +/- ".$stat4->standard_deviation()."\n";
print STDOUT "MAPQ == 5            : $qual5\t".$stat5->mean()." +/- ".$stat5->standard_deviation()."\n";
print STDOUT "MAPQ == 10           : $qual10\t".$stat10->mean()." +/- ".$stat10->standard_deviation()."\n";
print STDOUT "MAPQ > 20            : $qual20\t".$stat20->mean()." +/- ".$stat20->standard_deviation()."\n";
print STDOUT "MAPQ > 30            : $qual30\t".$stat30->mean()." +/- ".$stat30->standard_deviation()."\n";
print STDOUT "MAPQ > 40            : $qual40\t".$stat40->mean()." +/- ".$stat40->standard_deviation()."\n";
print STDOUT "MAPQ > 50            : $qual50\t".$stat50->mean()." +/- ".$stat50->standard_deviation()."\n";
print STDOUT "MAPQ > 60            : $qual60\t".$stat60->mean()." +/- ".$stat60->standard_deviation()."\n";


exit;
