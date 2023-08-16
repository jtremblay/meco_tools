#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use POSIX;

my $usage=<<'ENDHERE';
NAME:
breakLongFastaSeqs.pl

PURPOSE:

INPUT:
--infile <string>  : Sequence file
--seq_length <int> : Seq length in files.

OUTPUT:
STDOUT            : Sequence file with seq lengths no larger than --seq_length

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $seq_length);
my $verbose = 0;

GetOptions(
   'infile=s' 	  => \$infile,
   'seq_length=i'=> \$seq_length,
   'verbose' 	  => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $counter = 1;
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $length = length($curr->seq);
   #my $header = $curr->header;
   my($header) = $curr->header =~ /^(>\S+) /;
   my $curr_seq = $curr->seq;

   if($length >= $seq_length){
      my $j=0;
      my $num_chunks = ceil($length / $seq_length);
      for(my $i=1;$i<=$num_chunks;$i++){
          print STDOUT $header."=part$i=\n";
          print STDOUT substr($curr_seq, $j, ($seq_length - 1))."\n";
          $j = $j + $seq_length;
      
      }
   }else{
      print STDOUT $curr->output;
   }
}
exit;
