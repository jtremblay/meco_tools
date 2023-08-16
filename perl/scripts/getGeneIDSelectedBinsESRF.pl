#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
getGeneIDSelectedBinsESRF.pl

PURPOSE:

INPUT:
--infile_gene_id <string> : Txt file with gene ids to fetch cpms.
--infile_cpm <string>     : cpm matrix
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_gene_id, $infile_cpm);
my $verbose = 0;

GetOptions(
   'infile_gene_id=s' => \$infile_gene_id,
   'infile_cpm=s' 	 => \$infile_cpm,
   'verbose' 	       => \$verbose,
   'help'             => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile_gene_id) or die "Can't open $infile_gene_id\n";
while(<IN>){
   chomp;
   #my @row = split(/\t/, $_);
   my $gene_id = $_;
   $gene_id =~ s/\s+//g;
   $hash{$gene_id} = $gene_id;
}
close(IN);

#print STDERR Dumper(\%hash);

open(IN, "<".$infile_cpm) or die "Can't open $infile_cpm\n";
while(<IN>){
   chomp;
   if($. == 1){
      print STDOUT $_."\n";
      #print STDERR $_."\n";

   }else{
      my @row = split(/\t/, $_);
      my $gene_id = $row[0];
      if(exists $hash{$gene_id}){
         print STDOUT $_."\n";
      }
   }
}
close(IN);
