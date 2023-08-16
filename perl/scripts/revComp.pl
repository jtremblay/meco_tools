#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
revComp.pl

PURPOSE:
Output to STDOUT a revomped sequences 

INPUT:
--infile <string> : Sequence file
--rev_only        : set if only reverse is to be calculated.
--rna             : set if rna

OUTPUT:
STDOUT            : recomped sequences in fasta format.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $rev_only, $rna);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'rna'      => \$rna,
   'rev_only' => \$rev_only,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

sub reverse_complement_IUPAC {
   my $dna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $revcomp;
}

sub reverse_complement_IUPAC_rna {
   my $rna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($rna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   $revcomp =~ s/T/U/g;
   $revcomp =~ s/t/u/g;
   return $revcomp;
}


## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   if($rev_only){ 
      my $rev = reverse($curr->seq);
      print STDOUT $curr->header."\n".$rev."\n";
   }else{
      if($rna){
        my $revcomp = reverse_complement_IUPAC_rna($curr->seq);
        print STDOUT $curr->header."\n".$revcomp."\n";
      }else{
        my $revcomp = reverse_complement_IUPAC($curr->seq);
        print STDOUT $curr->header."\n".$revcomp."\n";
      }
   }
}

exit
