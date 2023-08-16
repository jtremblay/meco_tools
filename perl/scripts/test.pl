#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
test.pl

PURPOSE:

INPUT:
--infile <string>  : infile
--prefix <string>  : string
OUTPUT:
--outfile <string> : outfile

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $prefix);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'outfile=s' => \$outfile,
   'prefix=s'  => \$prefix,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
open(OUT, ">".$outfile) or die "Can't open $outfile\n";
while(<IN>){
   chomp;
   print OUT $_."\n";
}
close(IN);
close(OUT);
