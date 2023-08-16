#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateKeggModuleMatrix.pl

PURPOSE:

INPUT:
--infile <string> : parsed blastp output
				
OUTPUT:
STDOUT

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
   my $module_id = row[8];


}
