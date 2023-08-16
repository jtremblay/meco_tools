#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
rnammerToBed.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT            : Bed file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
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
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    if(length($curr->seq) == 0){
        next;
    }
    #k121_149544 arc_ssu coords:61-546
    #my ($header, $start, $end) = $curr->header =~ /^>.*_(contig-\d+)_(\d+)-(\d+)_.*/;
    my ($header, $start, $end) = $curr->header =~ /^>(\S+) .*coords:(\d+)-(\d+)$/;
    print STDOUT "$header\t$start\t$end\n";
}
exit;
