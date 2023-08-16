#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
correctSilva.pl

PURPOSE:
Remove g__BacteriaKI and g__ArchaeaKI

INPUT:
--infile_fasta <string>                             : Sequence file.
 
OUTPUT:
<STDOUT>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'                              => \$infile_fasta,
   'verbose'                                     => \$verbose,
   'help'                                        => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    if(!($header =~ m/g__BacteriaKI/) && 
        !($header =~ m/g__ArchaeaKI/) && 
        !($header =~ m/proteobacteria-gammaproteobacteria-gammaproteobacteria incertae sedis/) &&
        !($header =~ m/RickettsialesOR/)) {
        print STDOUT $header."\n".$curr->seq."\n";
    }
}
exit;
