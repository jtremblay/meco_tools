#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
validateMapAndBarcodes.pl

PURPOSE:

INPUT:
--infile_barcodes <string>     : Sequence file
--infile_mapping_file <string> : mapping file in Qiime format
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_barcodes, $infile_mapping);
my $verbose = 0;

GetOptions(
   'infile_barcodes=s'     => \$infile_barcodes,
   'infile_mapping_file=s' => \$infile_mapping,
   'verbose'               => \$verbose,
   'help'                  => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(IN, "<".$infile_mapping) or die "Can't open $infile_mapping\n";
my $numberOfColumns = 0;
my $numberOfSamplesMap = 0;

my %hash;
my $redFlag = 0;
while(<IN>){
   chomp;
   if($. == 1){
      # First line
      my @row = split(/\t/, $_);

      my $first_el = $row[0];
      if($first_el ne "#SampleID"){
         die "first column header MUST be labeled #SampleID. Found value was: $first_el\n";
      }
      
      #Check if header ends with $ and not a \t char!
      if($_ =~ m/.*\t$/g){
         die "Header line ends with a \\t char! Please remove the \\t char\n";
      }

      $numberOfColumns = @row;
      foreach my $header(@row){
         if($header =~ m/\s+/){
            print STDERR "Field '$header' in header of mapping file contains a whitespace! Fix that before launching the pipeline!\n";
            $redFlag = 1;
         }
      }
   }elsif($_ =~ m/^#/){
       next;
   }else{
      my @row = split(/\t/, $_);
      if($numberOfColumns != @row){
         print STDERR "Number of columns for sample '".$row[0]."' does not equal the number of field in the header (".$numberOfColumns." fields...)\n";
         $redFlag = 1;
      }
      if(exists $hash{$row[0]}){
         print STDERR "Sample ".$row[0]." is present multiple times in your mapping file...\n";
         $redFlag = 1;
      }else{
          $hash{$row[0]} = $_;
      }
      #if($row[0] =~ m/-+/){
      #   print STDERR "Sample labeled '$row[0]' in mapping file contains a '-' character! Replace that with a '_' or '.' before launching the pipeline!\n";
      #   $redFlag = 1;
      #}

      $numberOfSamplesMap++;
   }
}
close(IN);

my %hash_barcodes;
my $numberOfSamplesBarcodes = 0;
my $ref_fasta_db = Iterator::FastaDb->new($infile_barcodes) or die("Unable to open Fasta file, $infile_barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/^>//;
   #if($header =~ m/-+/){
   #   print STDERR "Sample labeled '$header' in barcodes file contains a '-' character! Replace that with a '_' or '.' before launching the pipeline!\n";
   #   $redFlag = 1;
   #}
   
   if(exists $hash_barcodes{$header}){
      print STDERR "Sample ".$header." is present multiple times in your barcodes file...\n";
      $redFlag = 1;
   }else{
      $hash_barcodes{$header} = $curr->seq;
   }
   
   # Check if sample is in mapping file.
   if(exists $hash{$header}){

   }else{
      print STDERR "In barcodes files, sample ".$header." is not present in corresponding mapping_file...\n";
      $redFlag = 1;
   }

   $numberOfSamplesBarcodes++;
}

# Then check if sample from mapping file exists in barcode file

for my $key (keys %hash){
   if(exists $hash_barcodes{$key}){
      #print STDERR
   }else{
      print STDERR "Entry ".$key." exists in mapping file, but not in barcodes file...\n";
      $redFlag = 1;
   }
}

if($redFlag == 1){
   die "Something is wrong with your mapping file. Please fix and rerun pipeline wrapper.\n";
}else{
   print STDERR "All seems good, proceed and run the pipeline wrapper...\n";
}
exit;
