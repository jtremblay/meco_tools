#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
rdpToTaxonomy.pl

PURPOSE:

INPUT:
--infile <string>  : Sequence file.
--cutoff <float>   : RDP score cutoff to reconstitute lineages.
                     lineages will be reconstructed up to the 
                     --cutoff value. Default = 0.5.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $cutoff);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'cutoff=f'  => \$cutoff,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

$cutoff = 0.5 unless($cutoff);

my %hash;

## MAIN

print STDOUT "contig_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   #my ($id) = $row[0] =~ m/(rRNA_contig-\d+)/;
   my $id = $row[0];
   # 5,6,7 = k__bact, Kingdom, score
   # 8,9,10 = phy
   # 11,12,13 = class
   # 14,15,16 = order
   # 17,18,19 = family
   # 20,21,22 = genus
   # 23,24,25 = specie
   
   my $lineage = "";
   if($row[7] >= $cutoff){
      $lineage .= "$row[5]";
   }else{
      $lineage .= "";
   }

   if($row[10] >= $cutoff){
      $lineage .= "\t$row[8]";
   }else{
      $lineage .= "\t";
   }

   if($row[13] >= $cutoff){
      $lineage .= "\t$row[11]";
   }else{
      $lineage .= "\t";
   }

   if($row[16] >= $cutoff){
      $lineage .= "\t$row[14]";
   }else{
      $lineage .= "\t";
   }

   if($row[19] >= $cutoff){
      $lineage .= "\t$row[17]";
   }else{
      $lineage .= "\t";
   }

   if($row[22] >= $cutoff){
      $lineage .= "\t$row[20]";
   }else{
      $lineage .= "\t";
   }

   if($row[25] >= $cutoff){
      $lineage .= "\t$row[23]";
   }else{
      $lineage .= "\t";
   }
   next if $lineage eq "";

   $id =~ s/rRNA_//;
   # Also remove k__, p__ etc prefixes.
   $lineage =~ s/k__//g;
   $lineage =~ s/p__//g;
   $lineage =~ s/c__//g;
   $lineage =~ s/o__//g;
   $lineage =~ s/f__//g;
   $lineage =~ s/g__//g;
   $lineage =~ s/s__//g;
   print STDOUT $id."\t".$lineage."\n";

}

exit;
