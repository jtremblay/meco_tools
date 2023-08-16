#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
filterFastaByHeader.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file.
--terms <string>  : text file with one entry per line to match in fasta headers.
				
OUTPUT:
<string>          : Output sequence file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $terms);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'terms=s'  => \$terms,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

my @array;
open(IN, "<".$terms) or die "Can't open $terms\n";
while(<IN>){
   chomp;
   push(@array, $_);
}
close(IN);

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $found = "false";

    foreach my $el (@array){
        if($curr->header =~ m/$el/ && $found eq "false"){
            print STDOUT $curr->output;
            my $found = "true";
        }
    }
}
exit;

