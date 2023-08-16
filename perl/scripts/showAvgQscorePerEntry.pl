#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use List::Util qw(sum);
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
tagsQC.pl

PURPOSE:
Filter amplicon reads.

INPUT:
--infile <string>          : Sequence file. Assumes it is fwd reads for single-end reads data or
                             merged/assembled amplicon fastq.

OUTPUT:

OPTIONAL:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
  'infile=s'                   => \$infile,
  'verbose'                    => \$verbose,
  'help'                       => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--infile $infile might be empty or wrong file path?\n") if((!-e $infile) and (!-s $infile));

my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fastq file, $infile\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
    my @qual = unpack("C*", $curr->qual);
    my $prob_sum = 0;
    
    foreach(@qual){
       $prob_sum = $prob_sum + ($_ - 33);
    }
  
    my $average = $prob_sum/@qual;
  
    my $N_count = 0;
    
    #foreach(@seq){
    #   $N_count++ if uc($_) eq "N";
    #}
    print STDOUT $curr->header."\t".$average."\n";
}


sub reverse_complement_IUPAC {
   my $dna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $revcomp;
}

sub complement_IUPAC {
   my $dna = shift;


   # complement the reversed DNA sequence
   my $comp = $dna;
   $dna =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $comp;
}

