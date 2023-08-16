#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--link <string>          : i.e. parsed bins table
--link2 <string>         : link between genes and contigs.
--selected_bins <string> : list of selected bins
--otu_table <string>     : otu table format file.
				
OUTPUT:
STDOUT                   : otu_table of selected bins only

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $link, $link2, $selected_bins, $otu_table);
my $verbose = 0;

GetOptions(
   'link=s' 	      => \$link,
   'link2=s' 	      => \$link2,
   'selected_bins=s' => \$selected_bins,
   'otu_table=s' 	   => \$otu_table,
   'verbose' 	      => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash_bin_to_contig;
my %hash_contig_to_gene;
my %hash_bin;
my %hash_otu_table;

open(IN, "<".$otu_table) or die "Can't open $otu_table\n";
while(<IN>){
   chomp;
   
   if($_ =~ m/^#/){
      print STDOUT $_."\n";
   }else{
      my @row = split(/\t/, $_);
      my $gene_id = $row[0];
      $hash_otu_table{$gene_id} = $_;
   }
}
close(IN);

open(IN, "<".$link) or die "Can't open $link\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $bin_id = $row[0];
   my $contig_id = $row[1];
   $hash_bin_to_contig{$bin_id} .= $contig_id.";";
   #$hash_contig_to_gene{$contig_id} = $bin_id;
}
close(IN);

open(IN, "<".$link2) or die "Can't open $link2\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $gene_id = $row[1];
   $hash_contig_to_gene{$contig_id} .= $gene_id.";";
   #$hash_contig_to_gene{$contig_id} = $bin_id;
}
close(IN);

open(IN, "<".$selected_bins) or die "Can't open $selected_bins\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $bin_id = $row[0];
   $hash_bin{$bin_id} = "";
}
close(IN);

#print STDERR Dumper(\%hash_otu_table);
#exit;
foreach my $bin (keys %hash_bin_to_contig) {

   if(exists $hash_bin{$bin}){
      print STDERR "#".$bin."\n";
      my $contigs = $hash_bin_to_contig{$bin};
      my @contigs = split(/;/, $contigs);

      foreach my $contig (@contigs){
         my $genes = $hash_contig_to_gene{$contig};
         my @genes = split(/;/, $genes);
   
         foreach my $gene (@genes){
             #print STDERR "GENE: ".$gene."\n";
            if(exists $hash_otu_table{$gene}){
               print STDOUT $hash_otu_table{$gene}."\n"; 
               print STDERR $hash_otu_table{$gene}."\n"; 
            }
         }
      }
   }else{
       #print STDERR "Bin does not exists: ".$bin."\n";
   }
}
