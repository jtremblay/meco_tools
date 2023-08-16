#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use POSIX;

my $usage=<<'ENDHERE';
NAME:
splitFastaByFileSize.pl

PURPOSE:

INPUT:
--infile <string>            : Sequence file
--targeted_chunk_size <int>  : Approx file size of each chunk to target (Mb).
				
OUTPUT:
STDOUT <string>              : outfile were the number of required file chunks will
                               be written. 

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $targeted_chunk_size, $outdir);
my $verbose = 0;

GetOptions(
   'infile=s'              => \$infile,
   'targeted_chunk_size=s' => \$targeted_chunk_size, 
   'outdir=s'              => \$outdir,
   'verbose'               => \$verbose,
   'help'                  => \$help
);
if ($help) { print $usage; exit; }

## MAIN

# size in byte
my $size =  (stat $infile)[7];
my $chunk_size = $targeted_chunk_size * 1000000;

my $number_of_files = ceil($size/$chunk_size);
#print STDERR "Number of files: $number_of_files\n";

print STDOUT $number_of_files; 

