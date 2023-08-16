#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
validateGroopmBins.pl

PURPOSE:

INPUT:
--genes_to_contigs <string>    : gene_to_contigs tsv file
--fasta <string>               : fasta file of contigs (contigs belonging to one bin)
--taxonomy <string>            : taxonomy.tsv file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $genes_to_contigs, $fasta, $taxonomy);
my $verbose = 0;

GetOptions(
   'genes_to_contigs=s' => \$genes_to_contigs,
   'fasta=s'            => \$fasta,
   'taxonomy=s'         => \$taxonomy,
   'verbose' 	         => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--genes_to_contigs is missing\n" unless($genes_to_contigs);

## MAIN
my %hash_gene_to_contig;
my %hash_contig_to_gene;

# Parse genes_to_contigs
open(IN, "<".$genes_to_contigs) or die "Can't open $genes_to_contigs\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   my $gene_id = "gene_id_".$row[0];
   my $contig_id = "contig-".$row[1];

   $hash_gene_to_contig{$gene_id} = $contig_id;
   if(exists $hash_contig_to_gene{$contig_id}){
      $hash_contig_to_gene{$contig_id} .= "\t".$gene_id; 
   }else{
      $hash_contig_to_gene{$contig_id} = $gene_id;
   }
}
close(IN);

my %hash_tax;
open(IN, "<".$taxonomy) or die "Can't open $taxonomy\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   my $species = join(".", $row[8]);
   my $lineage = join("\t", $row[2], $row[3], $row[4], $row[5], $row[6], $row[7], $species);
   
   if(exists $hash_gene_to_contig{$gene_id}){
      my $contig_id = $hash_gene_to_contig{$gene_id};
      $hash_tax{$contig_id}{$lineage}++;
   }#else{
   #  $hash_tax{$contig_id} = 1;
   #}
}
close(IN);

# Then parse fasta.
my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my ($contig_id) = $curr->header =~ /^>(\S+)/;
   #$contig_id = $curr->header;
   $contig_id =~ s/^>//;
   my $tax;
   my $max = 0;

   if(exists $hash_tax{$contig_id}){   
      foreach my $lineage (keys %{ $hash_tax{$contig_id} }){
         if($hash_tax{$contig_id}{$lineage} > $max){
            $tax = $lineage;
            $max = $hash_tax{$contig_id}{$lineage};
         }
      } 
   }else{
      $tax = "NA";
   }
   print STDOUT "$contig_id\t$tax\n";
}

exit;
