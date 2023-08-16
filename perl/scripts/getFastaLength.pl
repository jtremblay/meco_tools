#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
getFastaLength.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file				
--gc              : if set will compute GC percentage.

OUTPUT:
STDOUT            : header\tlength\n and gc if applicable.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $gc);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'gc'       => \$gc,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header();
    if($header =~ m/^(\S+)/){
        $header = $1;
    }
    $header =~ s/>//;
    if($gc){
        my $length = length($curr->seq);
        my $gcCount = 0;
        $gcCount += ($curr->seq =~ tr/gGcC/gGcC/);
        my $gcPerc = sprintf( "%0.2f", ($gcCount / $length) * 100);
        
        print STDOUT $header."\t".length($curr->seq)."\t".$gcPerc."\n";

    }else{
        print STDOUT $header."\t".length($curr->seq)."\n";
    }
}

exit;
