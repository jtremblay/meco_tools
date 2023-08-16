#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
filterFastaByLength.pl

PURPOSE:

INPUT:
--infile <string>  : Sequence file
--length <int>     : length

OUTPUT:
STDOUT             : higher or equal than <length>
--lower            : lower than <length>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $lower, $length);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'length=i' => \$length,
   'lower=s'  => \$lower,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(LOWER, ">".$lower) or die "Can't open $lower\n" if($lower);

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   if(length($curr->seq) >= $length){
      print STDOUT $curr->header."\n".$curr->seq."\n";
   }else{
      print LOWER $curr->header."\n".$curr->seq."\n" if($lower);
   }
}

