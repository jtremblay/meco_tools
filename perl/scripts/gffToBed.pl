#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
gffToBed.pl

PURPOSE:

INPUT:
--infile <string> : gff file (output from Prodigal)

OUTPUT:
STDOUT            : Bed file.

NOTES:
2017-09-03 JT: modified it to adapt to output of Prodigal instead of metagenemark

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   my $start = $row[3];
   my $end = $row[4];
   #my($contig_id) = $row[0] =~ m/(contig-\d+) /;
   my $contig_id = "";
   if($row[0] =~ m/(contig-\d+) /){
      $contig_id = $1; 
   }else{
      if($row[0] =~ m/(\S+\d+)/){$contig_id = $1;}
   }
   
   #my($gene_id) = $row[8] =~ m/gene_id=(\d+),/;
   my $gene_id = "";
   if( $row[8] =~ m/gene_id_(\d+);/){
      $gene_id = $1; 
   }else{
      if($row[8] =~ m/ID=(\S+);/){$gene_id = $1;}
   }

   print STDOUT $contig_id."\t".$start."\t".$end."\tgene_id_".$gene_id."\n";
}
close(IN);
exit;
