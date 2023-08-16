#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
fetchNRAnnotation.pl

PURPOSE:
From nr blast spreadsheet results, fetch product name
from nr.tsv.

INPUT:
--blast <string>       : blast output (outfmt 6)
--annotations <string> : nr annotation spreadsheet. 2 columns.
				
OUTPUT:
STDOUT                 : blast output annotated.

NOTES:
nr.tsv was generated from nr.fasta. For each accession number, 
it contains many productname (many identical sequence). Only
the first was kept in nr.tsv.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $blast, $annotations);
my $verbose = 0;

GetOptions(
   'blast=s'       => \$blast,
   'annotations=s' => \$annotations,
   'verbose' 	    => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$blast) or die "Can't open $blast\n";
while(<IN>){
   chomp;
   my($gi) = $_ =~ m/(gi\|\d+\|)/;
   $hash{$gi} = $_;
}
close(IN);

open(IN, "<".$annotations) or die "Can't open $annotations\n";
while(<IN>){
   chomp;
   #my($gi) = $_ =~ m/gi\|\d+\|)/;
   my @row = split(/\t/, $_);
   my $gi = $row[0];
   my $annotation = $row[1];
   if(exists $hash{$gi}){
      print STDOUT $hash{$gi}."\t".$annotation."\n"; 
   }
}
close(IN);







