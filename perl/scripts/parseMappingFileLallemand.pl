#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

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
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp $_;

   if($. == 1){
      print STDOUT $_."\n";
   }

   my @row = split(/\t/, $_);
   my $number = $row[0];
   my $letter = $row[1];
   for(my $i=2; $i<6; $i++){
      print STDOUT $number."v".$i."\.".$letter."\t".join("\t", @row[1..(@row-1)])."\n";
   }
}
close(IN);
exit;
