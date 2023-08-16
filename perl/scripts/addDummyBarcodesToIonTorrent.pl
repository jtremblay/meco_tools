#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
addDummyBarcodesToIonTorrent.pl

PURPOSE:

INPUT:
--infile <string>  : Sequence file fastq
--barcode <string> : string of barcode
				
OUTPUT:
STDOUT

NOTES:
Simple script to add a dummy barcode sequence in the header of barcodeless
iontorrent fastq data. We assume that fastq contains only one sample of course...

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $barcode);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'barcode=s' => \$barcode,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fastq file, $infile\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
   my $header = $curr->header;
   my $new_entry = $header."/1#".$barcode."\n".$curr->seq."\n+\n".$curr->qual."\n";
   print STDOUT $new_entry;
}
exit;
