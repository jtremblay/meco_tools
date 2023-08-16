#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:
For the project of generating synthetic bacterial community. This
script will pick <n>random genome for the specified <taxon>.

INPUT:
--infile <string> : Sequence file
--n <int>         : number of each taxon we want.
--col <int>       : column to filter for these values:
--taxon <string>  : Phylum, Class, Order, Family or Genus.
				
OUTPUT:
STDERR            : JGI/IMG ID number so we can then fetch the genome
                    from the IMG db (genepool/NERSC).

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $n, $col);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'n=i'       => \$n,
   'col=i'     => \$col,
#   'taxon=s'   => \$taxon,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
die("--infile is missing\n") unless $infile;
die("--n is missing\n") unless $n;
#die("--taxon is missing\n") unless $taxon;
die("--col is missing\n") unless $col;

my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   if($. == 1){
      next;
   }
   my @row = split(/\t/, $_);
   my $toid = $row[0];
   my $currentTaxon = $row[$col];
   
   if(exists $hash{$currentTaxon}){
      if($hash{$currentTaxon}{abundance} < $n){
         $hash{$currentTaxon}{abundance}++;
         print STDOUT "$row[0]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\n";
      }
   }else{
      $hash{$currentTaxon}{abundance} = 1;
      print STDOUT "$row[0]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\n";
   }
}
