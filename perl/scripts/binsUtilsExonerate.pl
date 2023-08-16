#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
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
   'indir=s'              => \$indir,
   'targeted_chunk_sie=s' => \$targeted_chunk_size, 
   'verbose'              => \$verbose,
   'help'                 => \$help
);
if ($help) { print $usage; exit; }

## MAIN


open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);

}
