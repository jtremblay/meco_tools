#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Storable;

my $usage=<<'ENDHERE';
NAME:
storeNcbiTaxonomyInHash.pl

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
my ($help, $infile, $outfile);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'outfile=s' => \$outfile,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN


