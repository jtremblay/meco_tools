#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
wheatMappingFile.pl

PURPOSE:

INPUT:
--infile_meta <string>      : metadata provided by saskatoon.
--infile_barcodes <string>  : barcodes after preprocessing.
				
OUTPUT:
--outfile_barcodes <string> : outfile barcodes
--outfile_map <string>      : outfile mapping_file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_meta, $infile_barcodes, $outfile_barcodes, $outfile_map);
my $verbose = 0;

GetOptions(
   'infile_meta=s'      => \$infile_meta,
   'infile_barcodes=s'  => \$infile_barcodes,
   'outfile_barcodes=s' => \$outfile_barcodes,
   'outfile_map=s'      => \$outfile_map,
   'verbose' 	         => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

open(OUT_BARCODES, ">$outfile_barcodes") or die "Can't open $outfile_barcodes\n";
open(OUT_MAP, ">$outfile_map") or die "Can't open $outfile_map\n";

my $ref_fasta_db = Iterator::FastaDb->new($infile_barcodes) or die("Unable to open Fasta file, $infile_barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   $hash{$header} = $curr->seq;
}

open(META, "<".$infile_meta) or die "Can't open $infile_meta\n";
while(<META>){
   chomp;
   my $date;
   my @row = split(/\t/, $_);
   
   $date = "jul2013" if($row[3] eq "2013-07-03");
   $date = "aug2013" if($row[3] eq "2013-08-26");
   $date = "jul2014" if($row[3] eq "2014-07-07");
   $date = "sep2014" if($row[3] eq "2014-09-03");
   
   my $mtl_sample_name = $row[4].".".$row[1].$date;
   my $sask_sample_name = $row[0];
   $sask_sample_name =~ s/[-|_|-|\s]//g;

   if(exists $hash{$mtl_sample_name}){
      print OUT_BARCODES ">".$sask_sample_name."\n".$hash{$mtl_sample_name}."\n";
      print OUT_MAP $sask_sample_name."\t".$row[1]."\t".$row[2]."\t".$row[3]."\t".$row[4]."\t".$row[5]."\n";   
   }else{
      print STDERR "Sask sample name does not exists in mtl barcodes: ".$sask_sample_name." ... ".$mtl_sample_name." \n";
   
   }
}
close(META);


