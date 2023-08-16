#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
linkBinsToGenes.pl

PURPOSE:

INPUT:
--link <string>          : i.e. parsed bins table
--link2 <string>         : link between genes and contigs.
--selected_bins <string> : list of selected bins

OUTPUT:
STDOUT                   : tsv with 3 columns: bin, contig, gene

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $link, $link2, $selected_bins);
my $verbose = 0;

GetOptions(
   'link=s'           => \$link,
   'link2=s'          => \$link2,
   'selected_bins=s'  => \$selected_bins,
   'verbose'          => \$verbose,
   'help'             => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash_bin_to_contig;
my %hash_contig_to_gene;
my %hash_bin;

open(IN, "<".$link) or die "Can't open $link\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $bin_id = $row[0];
   my $contig_id = $row[1];
   $hash_bin_to_contig{$bin_id} .= $contig_id.";";
}
close(IN);

#print STDERR Dumper(\%hash_bin_to_contig);

open(IN, "<".$link2) or die "Can't open $link2\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   $contig_id =~ s/contig-//g;
   my $gene_id = $row[1];
   $hash_contig_to_gene{$contig_id} .= $gene_id.";";
}
close(IN);
#print STDERR Dumper(\%hash_contig_to_gene);

open(IN, "<".$selected_bins) or die "Can't open $selected_bins\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $bin_id = $row[0];
   $hash_bin{$bin_id} = "";
}
close(IN);
#print STDERR Dumper(\%hash_bin);

#exit;

foreach my $bin (keys %hash_bin_to_contig) {
   if(exists $hash_bin{$bin}){
      my $contigs = $hash_bin_to_contig{$bin};
      my @contigs = split(/;/, $contigs);

      foreach my $contig (@contigs){

         my $genes = $hash_contig_to_gene{$contig};
         my @genes = split(/;/, $genes);
   
         foreach my $gene (@genes){
            print STDOUT "$bin\t$contig\t$gene\n";
         }
      }
   }else{
       #print STDERR "Bin does not exists: ".$bin."\n";
   }
}
