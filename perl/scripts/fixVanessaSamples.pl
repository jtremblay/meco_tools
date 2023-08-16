#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fixVanessaSamples.pl

PURPOSE:

INPUT:
--infile_R1 <string>  : Sequence file
--infile_R2 <string>  : Sequence file
--fwd_primer <string> : Forward primer fasta file
				
OUTPUT:
--outfile_R1 <string> : fastq
--outfile_R2 <string> : fastq

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2, $outfile_R1, $outfile_R2, $fwd_primer);
my $verbose = 0;

GetOptions(
   'infile_R1=s'  => \$infile_R1,
   'infile_R2=s'  => \$infile_R2,
   'outfile_R1=s'	=> \$outfile_R1,
   'outfile_R2=s'	=> \$outfile_R2,
   'fwd_primer=s' => \$fwd_primer,
   'verbose'      => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT1, ">".$outfile_R1) or die "Can't open $outfile_R1\n";
open(OUT2, ">".$outfile_R2) or die "Can't open $outfile_R2\n";


my $primer_seq;
my $ref_fasta_db = Iterator::FastaDb->new($fwd_primer) or die("Unable to open Fasta file, $fwd_primer\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   $primer_seq = $curr->seq;
}
$primer_seq = substr($primer_seq, 3, (length($primer_seq) -6));
print STDERR "Primer seq to use: ".$primer_seq."\n";

my $ref_fastq_db1 = Iterator::FastqDb->new($infile_R1) or die("Unable to open Fastq file, $infile_R1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($infile_R2) or die("Unable to open Fastq file, $infile_R2\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
   my $curr2 = $ref_fastq_db2->next_seq();

   if(substr($curr1->seq, 0, 40) =~ m/$primer_seq/){
      print OUT1 $curr1->output;
      print OUT2 $curr2->output;
   }
}

close(OUT1);
close(OUT2);

