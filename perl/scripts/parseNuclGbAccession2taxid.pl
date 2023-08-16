#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseNuclGbAccession2taxid.pl

PURPOSE:
Parse nucl_gb_accession2taxid file. bring 4 columns to 2 columns.
Keep accession.version and taxid columns only.

INPUT:
--infile <string> : Input file
				
OUTPUT:
STDOUT <string>   : Parsed file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   print STDOUT $row[1]."\t".$row[2]."\n";

}
close(IN);
exit;
