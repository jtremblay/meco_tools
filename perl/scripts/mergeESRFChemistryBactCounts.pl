#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
mergeESRFChemistryBactCounts.pl

PURPOSE:

INPUT:
--infile_chemistry <string> : spreadsheet chemistry.
--infile_counts <string>    : spreadsheet bacterial counts.
				
OUTPUT:
STDOUT                      : Same as chemistry spreadsheet, but with added bact. counts.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_chemistry, $infile_counts);
my $verbose = 0;

GetOptions(
   'infile_chemistry=s' => \$infile_chemistry,
   'infile_counts=s'    => \$infile_counts,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
die "Missing --infile_chemistry\n" unless $infile_chemistry;
die "Missing --infile_counts\n" unless $infile_counts;

open(IN, "<".$infile_chemistry) or die "Can't open $infile_chemistry\n";
open(IN_COUNTS, "<".$infile_counts) or die "Can't open $infile_counts\n";

my %hash;

while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $id = $row[7];
   $hash{$id} = $_;   
}
close(IN);

print STDERR Dumper(\%hash);

while(<IN_COUNTS>){
   chomp;
   my @row = split(/\t/, $_);
   my $id = $row[4];
   my $cell_count = $row[2];
   print STDERR "CELL COUNT: ".$cell_count."\n";
   print STDERR "ID: ".$id."\n";
   
   if(exists $hash{$id}){
      my @row2 = split(/\t/, $hash{$id});
      print STDOUT "$row[0]\t$row[1]\t$row[2]\tcell_count\t$row[4]\t$row[5]\t$cell_count\t$row[7]\t\n";
   }
}
close(IN_COUNTS);
