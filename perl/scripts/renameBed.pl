#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
renameBed.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
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
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $name = $row[0];
    my $start = $row[1];
    my $end = $row[2];
    my $abundance = $row[3];
    print STDOUT $name.":".$start."-".$end."\t".$start."\t".$end."\t".$abundance."\n";
}
close(IN);
exit;
