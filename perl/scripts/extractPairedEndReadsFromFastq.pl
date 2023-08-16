#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
extractPairedEndReadsFromFastq.pl

PURPOSE:
Intended for fastq holding empty fastq entries. Assumes all fastq 
entries are in the same orders for both R1 and R2.

INPUT:
--infile_R1 <string>  : Sequence file
--infile_R2 <string>  : Sequence file
				
OUTPUT:
--outfile_R1 <string> : Sequence file
--outfile_R2 <string> : Sequence file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2, $outfile_R1, $outfile_R2);
my $verbose = 0;

GetOptions(
   'infile_R1=s'  => \$infile_R1,
   'infile_R2=s'  => \$infile_R2,
   'outfile_R1=s' => \$outfile_R1,
   'outfile_R2=s' => \$outfile_R2,
   'verbose'      => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT1, ">".$outfile_R1) or die "Can't open $outfile_R1\n";
open(OUT2, ">".$outfile_R2) or die "Can't open $outfile_R2\n";

my $ref_fastq_db1 = Iterator::FastqDb->new($infile_R1) or die("Unable to open Fastq file, $infile_R1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($infile_R2) or die("Unable to open Fastq file, $infile_R2\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
    my $curr2 = $ref_fastq_db2->next_seq();

    if((length($curr1->seq) > 0) && (length($curr2->seq) > 0)){
        print OUT1 $curr1->output;
        print OUT2 $curr2->output;
    }else{
        print STDERR $curr1->output;
        print STDERR length($curr1->seq)."\n";
        print STDERR $curr2->output;
        print STDERR length($curr2->seq)."\n";
    }
}
close(OUT1);
close(OUT2);
