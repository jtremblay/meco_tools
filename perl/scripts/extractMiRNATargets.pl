#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
extractMiRNATargets.pl

PURPOSE:
From a blast output file. Ideally from blasting miRNA sequences against a reference sequence db.
extract the target sequence

INPUT:
--infile <string> : Sequence file
--gff    <string> : 

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

# First loop through blast file and only keep the perfect matches

my %hash;

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $qlen = $row[12];
    my $qstart = $row[6];
    my $qend = $row[7];
    my $sstart = $row[8];
    my $send = $row[9];
    my $alen = $row[3];
    my $qid = $row[0];
    my $sid = $row[1];

    if($alen < $qlen){
        next;
    }else{
        print STDOUT $_."\n";
    }
}



