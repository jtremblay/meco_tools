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
--suffix <string>    : R1 or R2

OUTPUT:

NOTES:
demultiplex454.pl --infile ./seq.fastq,./seq2.fastq --barcodes ./barcodes1.fa,./barcodes2.fa --suffix R1 --outdir ./outdir

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2, $barcodes, $outdir, $suffix);
my $verbose = 0;

GetOptions(
   'infile_R1=s'  => \$infile_R1,
   'infile_R2=s'  => \$infile_R2,
   'suffix=s'     => \$suffix,
   'barcodes=s'   => \$barcodes,
   'outdir=s'     => \$outdir,
   'verbose'      => \$verbose,
   'help'         => \$help
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

#my @infiles1 = split(/,/, $infile_R1);
#my @infiles2 = split(/,/, $infile_R2);
#my @barcodes = split(/,/, $barcodes);

#foreach my $file(@infiles){
    #my $barcode = shift(@barcodes);

   # Open fastq
my $ref_fastq_db1 = Iterator::FastqDb->new($infile_R1) or die("Unable to open Fastq file, $infile_R1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($infile_R2) or die("Unable to open Fastq file, $infile_R2\n");
while( my $curr_fastq1 = $ref_fastq_db1->next_seq() ) {
   my $curr_fastq2 = $ref_fastq_db2->next_seq();
   
   my $fastq_seq1 = $curr_fastq1->seq;
   my $fastq_qual1 = $curr_fastq1->qual;
   my $fastq_header1 = $curr_fastq1->header;
   my $fastq_seq2 = $curr_fastq2->seq;
   my $fastq_qual2 = $curr_fastq2->qual;
   my $fastq_header2 = $curr_fastq2->header;

   # Then loop through fasta barcodes and print were appropriate.
   my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my $barcode_seq = $curr->seq;
      my $header = $curr->header;
      $header =~ s/>//;

      # Then check if match R1
      my $found_R1 = "false";
      my $found_R2 = "false";
      my ($matched_seq1, $matched_seq2, $matched_qual1, $matched_qual2);
      if(substr($fastq_seq1, 0, 12) =~ m/$barcode_seq/i){
          my $last_match_pos = $+[0];
          #print STDERR $last_match_pos."\n";
          $matched_seq1 = substr($fastq_seq1, $last_match_pos);
          $matched_qual1 = substr($fastq_qual1, $last_match_pos);
          $found_R1 = "true";
      }#else{
      #   $matched_seq1 = $fastq_seq1;
      #    $matched_qual1 = $fastq_qual1;
      #    $found_R1 = "false";
      #}

      if(substr($fastq_seq2, 0, 12) =~ m/$barcode_seq/i){
          my $last_match_pos = $+[0];
          #print STDERR $last_match_pos."\n";
          $matched_seq2 = substr($fastq_seq2, $last_match_pos);
          $matched_qual2 = substr($fastq_qual2, $last_match_pos);
          $found_R2 = "true";
      }#else{
      #   $matched_seq2 = $fastq_seq1;
      #    $matched_qual2 = $fastq_qual1;
      #    $found_R2 = "false"
      #}

      if($found_R1 eq "true" && $found_R2 eq "true"){
          my $outfile1 = "./".$outdir."/".$header."_R1.fastq";
          open(OUT1, ">>".$outfile1) or die "Can't open $outfile1\n";
          print OUT1 "".$fastq_header1."#".$barcode_seq."/1\n".$matched_seq1."\n+\n".$matched_qual1."\n";
          close(OUT1);
          my $outfile2 = "./".$outdir."/".$header."_R2.fastq";
          open(OUT2, ">>".$outfile2) or die "Can't open $outfile2\n";
          print OUT2 "".$fastq_header2."#".$barcode_seq."/2\n".$matched_seq2."\n+\n".$matched_qual2."\n";
          close(OUT2);
          next;
      }
      #else{
      #die "Did not match in 5'...";
         #my $barcode_rc = revcomp($barcode_seq);
         #if( substr($fastq_seq, length($fastq_seq) - 40) =~ m/$barcode_rc/i){
         #   
         #    my $first_match_pos = $+[0];
         #    my $matched_seq = substr($fastq_seq, 0, length($fastq_seq) - $first_match_pos);
         #    my $matched_qual = substr($fastq_qual, 0, length($fastq_seq) - $first_match_pos);
         # 
         #    my $outfile = "./".$outdir."/".$header.".fastq";
         #    open(OUT, ">>".$outfile) or die "Can't open $outfile\n";
         #    print OUT "".$fastq_header."#".$barcode_seq."/".$suffix."\n".$matched_seq."\n+\n".$matched_qual."\n";
         #    close(OUT);
         #}
         #}
   }
}
#}


