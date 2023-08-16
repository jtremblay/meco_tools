#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;

my $usage=<<'ENDHERE';
NAME:
removeUnpairedReads.pl

PURPOSE:

INPUT:
--infile <string>          : Sequence file
				
OUTPUT:
STDOUT                     : Fastq in paired interleaved format.
--unpaired_reads1 <string> : Unpaired reads 1
--unpaired_reads2 <string> : Unpaired reads 2

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $reads1, $reads2);
my $verbose = 0;

GetOptions(
   'infile=s' 	        => \$infile,
   'unpaired_reads1=s' => \$reads1,
   'unpaired_reads2=s' => \$reads2,
   'verbose' 	        => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(FASTQ_OUT_1, '>'.$reads1) or die "Can't open $reads1\n";
open(FASTQ_OUT_2, '>'.$reads2) or die "Can't open $reads2\n";

my $tmpdir = File::Temp->newdir(
    "tmpdir-mergePairs-XXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

my $pipe1 = "$tmpdir/reads1.pipe";
system("mkfifo $pipe1");
system("gunzip -c ".$infile." > $pipe1 &");

my $pair_1=0;
my $pair_2=0;
my $counter=0;
my $seq_1;
my $seq_2;
my $seq_1_base;
my $seq_2_base;

my $base;
my $barcode;
my $pair;
my $header;
my $seq;
my $qual;
my $output;

my $ref_fastq_db = Iterator::FastqDb->new($pipe1) or die("Unable to open Fastq file, $pipe1\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
   $base = "@".$curr->base;   
   $barcode = $curr->barcode if(defined($curr->barcode));   
   $pair = $curr->pair;
   if(defined($curr->barcode)){
      $header = $base."#".$barcode."/".$pair;
   }else{
      $header = $base."/".$pair;
   }
   $seq = $curr->seq;
   $qual = $curr->qual;
   
   ## Perform task here
   unless($pair){
      print STDERR "================\n";
      print STDERR $header."\n".$seq."\n+\n".$qual."\n";
      print STDERR "================\n";
      die("No paring information is available in the fastq header\n");
   }
   $output = $header."\n".$seq."\n+\n".$qual."\n";

   if( !defined($pair) ) {
      die($seq->header." does not having pairing information\n");
      next;
   }elsif($pair == 1 ) {
      $pair_1++;
      $seq_1 = $output;
      $seq_1_base = $base;
      $counter++;
   }elsif($pair == 2) {
      $pair_2++;
      $seq_2 = $output;
      $seq_2_base = $base;
      $counter++;
   }else {
      die($header." has invalid pair ID, ".$pair."\n");
   }
   
   #validate   
   if($pair_1 == 1 && $pair_2 == 1 && $counter == 2 && $seq_1_base eq $seq_2_base){ #GOOD PAIRS!
      print STDOUT $seq_1.$seq_2;
      $pair_1=0;
      $pair_2=0;
      $counter=0;
      $seq_1="";
      $seq_2="";
      $seq_1_base=0;
      $seq_2_base=0;
   }elsif($pair_1 == 1 && $pair_2 == 1 && $counter == 2 && $seq_1_base ne $seq_2_base){ #GOOD PAIRS, but base id does not correspond
      print FASTQ_OUT_1 $seq_1;
      print FASTQ_OUT_2 $seq_2;
      $pair_1=0;
      $pair_2=0;
      $counter=0;
      $seq_1="";
      $seq_2="";
      $seq_1_base=0;
      $seq_2_base=0;
   }elsif($pair_1 == 1 && $pair_2 == 0 && $counter == 2){ #Not sure this can happen...
      print FASTQ_OUT_1 $seq_1;
      $pair_1=0;
      $pair_2=0;
      $counter=0;
      $seq_1="";
      $seq_1_base=0;
      $seq_2_base=0;
      next;   
   }elsif($pair_1 == 0 && $pair_2 == 1 && $counter == 2){ #Not sure this can happen...
      print FASTQ_OUT_2 $seq_2;
      $pair_1=0;
      $pair_2=0;
      $counter=0;
      $seq_2="";
      $seq_1_base=0;
      $seq_2_base=0;
      next;   
   }elsif($pair_1 == 0 && $pair_2 == 1 && $counter == 1){ #Problem, because pair 2 comes before pair 1. Means that there is no pair one.
      print FASTQ_OUT_2 $seq_2;
      $pair_1=0;
      $pair_2=0;
      $counter=0;
      $seq_2="";
      #$seq_1_base=0;
      $seq_2_base=0;
      next;
   #}elsif($pair_1 == 1 && $pair_2 == 0 && $counter == 1){
   #   print FASTQ_OUT_1 $seq_1;
   #   $pair_1=0;
   #   $pair_2=0;
   #   $counter=0;
   #   $seq_1="";
   #   $seq_1_obj=0;
   #   #$seq_2_obj=0;
   #   next;      
   }elsif($pair_1 == 2 && $counter == 2){
      print FASTQ_OUT_1 $seq_1;
      $pair_1=1;
      $pair_2=0;
      $counter=1;
      $seq_2="";
      #$seq_1_base=0;
      $seq_2_base=0;
      next;   
   }else{
      #DEBUGGING
      if($verbose){
         print "PAIR 1: ".$pair_1."\n";
         print "PAIR 2: ".$pair_2."\n";
         print "COUNTER: ".$counter."\n";
         print $seq_1."\n";
         print $seq_2."\n";
      }
   
   }
}
close FASTQ_OUT_1;
close FASTQ_OUT_2;
system("rm $pipe1");

exit;

