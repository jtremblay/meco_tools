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
--infiles <string> : Sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles);
my $verbose = 0;

GetOptions(
   'infiles=s' => \$infiles,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @infiles = split(/,/, $infiles);
my $R1 = $infiles[0]; #R1
my $R2 = $infiles[1]; #R2
my $R3 = $infiles[2]; #Index

my $ref_fastq_db_R1 = Iterator::FastqDb->new($R1) or die("Unable to open Fastq file, $R1\n");
my $ref_fastq_db_R2 = Iterator::FastqDb->new($R2) or die("Unable to open Fastq file, $R2\n");
my $ref_fastq_db_R3 = Iterator::FastqDb->new($R3) or die("Unable to open Fastq file, $R3\n");

while( my $curr1 = $ref_fastq_db_R1->next_seq() ) {
   my $curr2 = $ref_fastq_db_R2->next_seq();
   my $curr3 = $ref_fastq_db_R3->next_seq();

   my $barcode = $curr3->seq;

   print STDOUT $curr1->header.$barcode."\n";
   print STDOUT $curr1->seq."\n+\n".$curr1->qual."\n";

   my $header2 = $curr2->header;
   $header2 =~ s/ 3:/ 2:/;
   print STDOUT $header2.$barcode."\n";
   print STDOUT $curr2->seq."\n+\n".$curr2->qual."\n";

}
exit;
