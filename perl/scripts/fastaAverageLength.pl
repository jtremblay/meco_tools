#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(sum);
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Statistics::Descriptive;

my $usage=<<'ENDHERE';
NAME:
fastaAverageLength.pl

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
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }


## MAIN
my $stat = Statistics::Descriptive::Sparse->new();
my @array;
my $counter = 0;
my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fastq file, $infile\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
    my $length = length($curr->seq);
    push(@array, $length);
    $stat->add_data(length($curr->seq()));
    $counter++;
}

my $mean = sum(@array)/@array;
my $std = $stat->standard_deviation();
print STDOUT $mean." +/- ".$std."\n";

exit;
