#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
linkAnnotationFromGeneList.pl

PURPOSE:

INPUT:
--annotations <string>        : annotation master file.
--abundance <string>          : abundance master file. (should use normalized_significant.tsv) from DDA or DEG
                                Should contain only genes significantly differentially expressed or abundant.
--tax_group_to_match <string> : list of taxonomic groups to match. 
                                example: Alteromonadales,Oceanospirillales,Rhodobacterales,Pseudomonadales

OUTPUT:
STDOUT                        : annotations of genes matching tax_group_to_match
--outfile_abundance <string>  : abundance of genes matching tax_group_to_match

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $annotations, $abundance, $tax_group_to_match, $outfile_abundance);
my $verbose = 0;

GetOptions(
   'annotations=s'         => \$annotations,
   'abundance=s'           => \$abundance,
   'tax_group_to_match=s'  => \$tax_group_to_match,
   'outfile_abundance=s'   => \$outfile_abundance,
   'verbose' 	            => \$verbose,
   'help'                  => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(ANN, "<".$annotations) or die "Can't open $annotations\n";
open(ABUN, "<".$abundance) or die "Can't open $abundance\n";
open(OUT, ">".$outfile_abundance) or die "Can't open $outfile_abundance\n";

die "--annotations <string> missing...\n" unless($annotations);
die "--abundance <string> missing...\n" unless($abundance);
die "--tax_group_to_match <string> missing...\n" unless($tax_group_to_match);
die "--outfile_abundance <string> missing...\n" unless($outfile_abundance);

my %hash;
my %hash_significant;

$tax_group_to_match =~ s/,/|/g;

while(<ABUN>){
   chomp;
   if($_ =~ m/^#/){
      next;
   }
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   $hash_significant{$gene_id} = "";
   
}
close(ABUN);
open(ABUN, "<".$abundance) or die "Can't open $abundance\n";

while(<ANN>){
   chomp;
   if($. == 1){
      print STDOUT $_."\n";
   }else{
      my @row = split(/\t/, $_);
      my $gene_id = $row[1];

      if(exists $hash_significant{$gene_id}){
         if($_ =~ m/$tax_group_to_match/){
            print STDOUT $_."\n";
            $hash{$gene_id} = "";
         }
      }
   }
}
close(ANN);

while(<ABUN>){
   chomp;
   if($. == 1){
      print OUT $_."\n";
      next;
   }
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   
   if(exists $hash{$gene_id}){
      print OUT $_."\n"; 
   }
}
close(ABUN);
exit;

