#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
demultiplex454.pl

PURPOSE:

INPUT:
--infile <string>    : Sequence file fastq
--barcodes <string>  : fasta
--outdir <string>    : outdirectory

OUTPUT:

NOTES:
demultiplex454.pl --infile ./seq.fastq,./seq2.fasta --barcodes ./barcodes1.fa,./barcodes2.fa --outdir ./outdir

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $barcodes, $outdir);
my $verbose = 0;

GetOptions(
   'infile=s' 	 => \$infile,
   'barcodes=s' => \$barcodes,
   'outdir=s'   => \$outdir,
   'verbose' 	 => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## SUBS
sub revcomp {
        my $dna = shift;

   # reverse the DNA sequence
        my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


## MAIN
system("mkdir -p $outdir");

my @infiles = split(/,/, $infile);
my @barcodes = split(/,/, $barcodes);


foreach my $file(@infiles){
   my $barcode = shift(@barcodes);

   # Open fastq
   my $ref_fastq_db = Iterator::FastqDb->new($file) or die("Unable to open Fastq file, $file\n");
   while( my $curr_fastq = $ref_fastq_db->next_seq() ) {
      my $fastq_seq = $curr_fastq->seq;
      my $fastq_qual = $curr_fastq->qual;
      my $fastq_header = $curr_fastq->header;
   
      # Then loop through fasta barcodes and print were appropriate.
      my $ref_fasta_db = Iterator::FastaDb->new($barcode) or die("Unable to open Fasta file, $barcode\n");
      while( my $curr = $ref_fasta_db->next_seq() ) {
         my $barcode_seq = $curr->seq;
         my $header = $curr->header;
         $header =~ s/>//;

         # Then check if match
         if(substr($fastq_seq, 0, 40) =~ m/$barcode_seq/i){
             my $last_match_pos = $-[0];
             my $matched_seq = substr($fastq_seq, $last_match_pos);
             my $matched_qual = substr($fastq_qual, $last_match_pos);

             my $outfile = "./".$outdir."/".$header.".fastq";
             open(OUT, ">>".$outfile) or die "Can't open $outfile\n";
             print OUT ">".$fastq_header."-".$header."\n".$matched_seq."\n+\n".$matched_qual."\n";
             close(OUT);

         }else{
            my $barcode_rc = revcomp($barcode_seq);
            if( substr($fastq_seq, length($fastq_seq) - 80) =~ m/$barcode_rc/i){
               
                my $first_match_pos = $+[0];
                my $matched_seq = substr($fastq_seq, 0, length($fastq_seq) - $first_match_pos);
                my $matched_qual = substr($fastq_qual, 0, length($fastq_seq) - $first_match_pos);
   
                my $outfile = "./".$outdir."/".$header.".fastq";
                open(OUT, ">>".$outfile) or die "Can't open $outfile\n";
                print OUT ">".$fastq_header."-".$header."\n".$matched_seq."\n+\n".$matched_qual."\n";
                close(OUT);
            }
         }
      }
   }
}


