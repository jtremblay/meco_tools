#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseTaxonomySummary.pl

PURPOSE:

INPUT:
--infile <string>     : Taxonomy summary file as provided by Qiime.
--deepest <string>    : For instance genus. either 'specie', 'genus', 'family', 'order', 'class', 'phylum' or 'kingdom'
--shallowest <string> : For instance phylum.

OUTPUT:

NOTES:
k__Archaea;p__Crenarchaeota;c__Thermoprotei;o__Thermoproteales;f__Thermoproteaceae;g__Thermoproteus
with args:
  --deepest genus
  --shallowest phylym
would be transformed into this:
[p__Crenarchaeota];g__Thermoproteus

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $deepest, $shallowest);
my $verbose = 0;

GetOptions(
   'infile=s' 	   => \$infile,
   'deepest=s'    => \$deepest,
   'shallowest=s' => \$shallowest,
   'verbose' 	   => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## MAIN
die "--infile missing..." unless($infile);
die "--deepest missing..." unless($deepest);
die "--shallowest missing..." unless($shallowest);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   if($. == 1){
      print STDOUT $_."\n";
      next;
   }

   my @row = split(/\t/, $_);
   #my $lineage = $row[0];
   my $lineage = shift(@row);
   my @lineage = split(/;/, $lineage);

   my $curr_deep = "";
   my $curr_shallow = "";

   # Deep
   if($deepest eq "specie"){
      $curr_deep = $lineage[6];
   }elsif($deepest eq "genus"){
      $curr_deep = $lineage[5];
   }elsif($deepest eq "family"){
      $curr_deep = $lineage[4];
   }elsif($deepest eq "order"){
      $curr_deep = $lineage[3];
   }elsif($deepest eq "class"){
      $curr_deep = $lineage[2];
   }elsif($deepest eq "phylum"){
      $curr_deep = $lineage[1];
   }elsif($deepest eq "kingdom"){
      $curr_deep = $lineage[0];
   }else{
      die "--deepest not a valid value, see --help\n";
   }

   # Shallow 
   if($shallowest eq "specie"){
      $curr_shallow = $lineage[6];
   }elsif($shallowest eq "genus"){
      $curr_shallow = $lineage[5];
   }elsif($shallowest eq "family"){
      $curr_shallow = $lineage[4];
   }elsif($shallowest eq "order"){
      $curr_shallow = $lineage[3];
   }elsif($shallowest eq "class"){
      $curr_shallow = $lineage[2];
   }elsif($shallowest eq "phylum"){
      $curr_shallow = $lineage[1];
   }elsif($shallowest eq "kingdom"){
      $curr_shallow = $lineage[0];
   }else{
      die "--shallowest not a valid value, see --help\n";
   }

   print STDOUT "[$curr_shallow];$curr_deep\t".join("\t", @row)."\n";

}
close(IN);
exit;




