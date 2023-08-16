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
--infile_R1 <string> : Sequence file
--infile_R2 <string> : Sequence file
				
OUTPUT:
STDERR / STDOUT msgs

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2);
my $verbose = 0;

GetOptions(
   'infile_R1=s' => \$infile_R1,
   'infile_R2=s' => \$infile_R2,
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fastq_db1 = Iterator::FastqDb->new($infile_R1) or die("Unable to open Fastq file, $infile_R1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($infile_R2) or die("Unable to open Fastq file, $infile_R2\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
    my $curr2 =  $ref_fastq_db2->next_seq();
    my $base1 = $curr1->base;
    my $base2 = $curr2->base;
    if($base1 ne $base2){
        print STDERR $base1."\n";
        print STDERR $base2."\n";
        die("base1 not equal base2\n");
    }

}

