#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
extractAnnotationsFromGeneList.pl

PURPOSE:

INPUT:
--gene_list <string>   : Sequence file.
--annotations <string> : Annotation tsv file.
				
OUTPUT:
STDOUT                 : Selected genes (rows).

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $annotations, $gene_list);
my $verbose = 0;

GetOptions(
   'gene_list=s'   => \$gene_list,
   'annotations=s' => \$annotations,
   'verbose' 	    => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;

open(IN, "<".$gene_list) or die "Can't open $gene_list\n";
while(<IN>){
   chomp;
   $hash{$_} = "";

}
close(IN);

open(IN, "<".$annotations) or die "Can't open $annotations\n";
while(<IN>){
   chomp;
   if($. == 1){
     print STDOUT $_."\n";
   }
   my @row = split(/\t/, $_);
   if(exists $hash{$row[1]}){
      print STDOUT $_."\n";
   }

}
close(IN);
