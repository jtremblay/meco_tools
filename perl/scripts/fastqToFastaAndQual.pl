#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fastqToFastaAndQual.pl

PURPOSE:
The name says it all.

INPUT:
--fastq <fasta_infile>

OUTPUT:
--fasta <fasta_outfile> 
--qual <qual_outfile>

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@nrc-cnrc.gc.ca
ENDHERE

## OPTIONS
my ($help, $fastq, $fasta, $qual);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'fastq=s'   => \$fastq,
  'fasta=s'   => \$fasta,
  'qual=s'    => \$qual,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fastq fastq file required\n") unless $fastq;
die("--fasta fasta file required\n") unless $fasta;
die("--qual qual file required\n") unless $qual;

## MAIN
open(OUT_F, ">".$fasta) or die "Can't open file ".$!."\n";
open(OUT_Q, ">".$qual) or die "Can't open file ".$!."\n";

my $db = new Iterator::FastqDb($fastq);
while(my $curr=$db->next_seq) {
   my $header = $curr->header;
   #$header =~ s/#/-/g;
   $header =~ s/^@// if(substr($header,0,1) eq "@");
   print OUT_F ">".$header."\n".$curr->seq()."\n";
   my $qual = $curr->qual;
   my @qual = unpack("C*", $qual);
   print OUT_Q ">".$header."\n";
   my $length = @qual;
   my $i = 0;
   foreach my $curr_qual (@qual){
      if($i < $length - 1){
         print OUT_Q $curr_qual." "; 
      }else{
         print OUT_Q $curr_qual;
      }
      $i++;
   }
   print OUT_Q "\n";
}
close(OUT_F);
close(OUT_Q);



exit;
