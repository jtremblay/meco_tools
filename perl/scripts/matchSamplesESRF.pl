#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
matchSamplesESRF.pl

PURPOSE:
From 2 files having only one column with no headers. Determine the intersection between the 2s.

INPUT:
--infile1 <string> : tsv file
--infile2 <string> : tsv file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile1, $infile2);
my $verbose = 0;

GetOptions(
   'infile1=s' => \$infile1,
   'infile2=s' => \$infile2,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash1;
my %hash2;
open(IN, "<".$infile1) or die "Can't open $infile1\n";
while(<IN>){
   chomp;
   $hash1{$_} = $_; 
}
close(IN);

open(IN, "<".$infile2) or die "Can't open $infile2\n";
while(<IN>){
   chomp;
   $hash2{$_} = $_; 
}
close(IN);

# Then compare hashes
my $j=0;

for my $key (keys %hash1){
   if(exists $hash2{$key}){
         
   }else{
      print STDOUT "Present in $infile1 but absent in $infile2: $key\n";
      $j++;
   }
}
print STDOUT "$j samples are present in $infile1 but absent in $infile2\n";

$j = 0;
for my $key (keys %hash2){
   if(exists $hash1{$key}){
         
   }else{
      print STDOUT "Present in $infile2 but absent in $infile1: $key\n";
      $j++;
   }
}
print STDOUT "$j samples are present in $infile2 but absent in $infile1\n";

