#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
parseCanopies.pl

PURPOSE:

INPUT:
--infile <string> : Cluster tsv files output of canopy
--n <int>         : first n canopies to keep.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $n);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'n=i'       => \$n,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(IN, "<".$infile) or die "Can't open $infile\n";
my %hash;
my $counter = 0;
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   
   if(exists $hash{$row[0]}){
      $hash{$row[0]} = $row[1];
   }else{
      $hash{$row[0]} = $row[1];
      $counter++;
   }

   last if($counter > $n);
   print STDOUT $_."\n";
}

