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
--infile <string>    : metagenomics annotation master file (final file of the pipeline)
--contigs <string>   : contigs in fasta format.
--abundance <string> : gene abundance in RPKM.
				
OUTPUT:
STDOUT               : Table containing contig.id\tcontigLengt\tGC%\tabundance\tLineage

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $contigs, $abundance);
my $verbose = 0;

GetOptions(
   'infile=s' 	  => \$infile,
   'contigs=s'   => \$contigs,
   'abundance=s' => \$abundance,
   'verbose' 	  => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
#my %geneToContig;

my $ref_fasta_db = Iterator::FastaDb->new($contigs) or die("Unable to open Fasta file, $contigs\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my ($contig_id) = $curr->header =~ /^>(\S+)/;
   my $length = length($curr->seq);
   $hash{$contig_id}{LENGTH} = $length;
   
   my $gcCount = 0;
   $gcCount += ($curr->seq =~ tr/gGcC/gGcC/);
   my $gcPerc = sprintf( "%0.2f", ($gcCount / $length) * 100);
   
   $hash{$contig_id}{GC} = $gcPerc;
   
}

my %hash_tax;
my %hash_keggp;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $gene_id = $row[1];
   my $tax_genus = $row[26];
   my $kegg_pathway = $row[7];

   $hash_keggp{$contig_id} = $kegg_pathway;
   $hash_tax{$contig_id}{$tax_genus}++;
   #$geneToContig{$gene_id} = $contig_id;
}
close(IN);

my %hash_abundance;
my $abundance_header = "";
open(IN, "<".$abundance) or die "Can't open $abundance\n";
while(<IN>){
   chomp;
   if($. == 1){
       my @row = split(/\t/, $_);
       shift(@row);
       $abundance_header = join("\t", @row);
   }else{
      my @row = split(/\t/, $_);
      my $contig_id = $row[0];
      shift(@row);
   
      my @formatted_row;
      foreach my $el (@row){
         my $rounded = sprintf("%.3f", $el);
         push @formatted_row, $rounded;
      }
   
      my $abundances = join("\t", @formatted_row);
      $hash_abundance{$contig_id} = $abundances;
   }
}
close(IN);

#print STDERR Dumper(\%hash);
#print STDERR Dumper(\%hash_tax);
#print STDERR Dumper(\%hash_abundance);

print STDOUT  "contigId\tcontigLength\tGCperc\tKEGGPathway\tLineageOrder\t$abundance_header\n";
#my $i=0;
for my $contig_id (keys %hash){
   # Extract taxonomy with higest occurence for each contigs
   my $tax;
   my $max = 0;

   if(exists $hash_tax{$contig_id}){
      for my $tax_key (keys %{ $hash_tax{$contig_id} }){
         #print STDERR $tax_key."\n";
         #print STDERR $hash_tax{$contig_id}{$tax_key}."\n";
         if($hash_tax{$contig_id}{$tax_key} > $max){
            $tax = $tax_key;
            $max = $hash_tax{$contig_id}{$tax_key};
         }
      }
   }else{
      $tax = "NA";
   }
   
   my $abundance_value = "NA";
   if(exists $hash_abundance{$contig_id}){
      $abundance_value = $hash_abundance{$contig_id};
   }

   my $kegg_pathway_value = "NA";
   if(exists $hash_keggp{$contig_id}){
      $kegg_pathway_value = $hash_keggp{$contig_id};
   }

   print STDOUT $contig_id."\t".$hash{$contig_id}{LENGTH}."\t".$hash{$contig_id}{GC}."\t".$kegg_pathway_value."\t".$tax."\t".$abundance_value."\n";
   
   #exit if($i >= 100);
   #$i++; 
}

