#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateGCtable.pl

PURPOSE:

INPUT:
--infile <string>    : gc table.
--minLength <int>    : contig min length to include in the json file.
				
OUTPUT:
STDOUT               : same as gc table but in json format.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $minLength);
my $verbose = 0;

GetOptions(
   'infile=s' 	  => \$infile,
   'minLength=s' => \$minLength,
   'verbose' 	  => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

my $nlines = 0;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){ 
   chomp;
   if($. == 1){
      my @row = split(/\t/, $_);
   }else{
      my @row = split(/\t/, $_);
      $nlines++ if $row[1] >= $minLength;
   }
}
close(IN);

open(IN, "<".$infile) or die "Can't open $infile\n";
my @header;
print STDOUT '{'."\n".'   "nodes":['."\n";
my $j = 0;
while(<IN>){
   chomp;
   if($. == 1){
      my @row = split(/\t/, $_);
      @header = @row;
   }else{
      my @row = split(/\t/, $_);
      next if $row[1] < $minLength;

      my $i = 0;
      print STDOUT '      {';
      foreach my $el (@row){
         if($i < @header - 1 ){

            if($i == 3){ # KEGG PATHWAYS
               print STDOUT '"'.$header[$i].'":[';
               my @keggs = split(/==/, $el);
               my $w = 0;
               foreach my $kegg (@keggs){
                  if($w < @keggs - 1){
                     print STDOUT '"'.$kegg.'",';
                 }else{
                     print STDOUT '"'.$kegg.'"';
                  }
                  $w++;
               }
               print STDOUT '], ';
            
            }else{
               print STDOUT '"'.$header[$i].'":"'.$el.'", ';
            }
         }else{
            print STDOUT '"'.$header[$i].'":"'.$el.'"';
         }
         $i++;
      }
      if($j < $nlines - 1){
         print STDOUT '},';
      }else{
         print STDOUT '}'; 
      }
      print STDOUT "\n";
      $j++;
   }
}
print STDOUT "   ]\n}";
close(IN);
