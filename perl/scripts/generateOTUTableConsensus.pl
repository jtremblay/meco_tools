#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateOTUTableConsensus.pl

PURPOSE:
Generate a consensus lineage based on the hits of contig sequences vs NCBI nt and gene sequences vs NCBI nr.

INPUT:
--infile_tax_contigs <string>     : contigs taxonomy 
--infile_tax_genes <string>       : genes taxonomy
--infile_abundance <string>       : abundance matrix genes.
--infile_gff <string>             : file to link contigs to genes.
--infile_blastn <string>          : Blastn results of Contigs vs NCBI nt 
--infile_blastp <string>          : Diamond Blastp results of Genes vs NCBI nr

OUTPUT:
STDOUT <string>                   : Consensus OTU table
--outfile_taxonomy <string>       : Consensus taxonomy file.


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_tax_contigs, $infile_tax_genes, $infile_abundance, $infile_gff, $infile_blastn, $infile_blastp, $taxonomy_outfile);
my $verbose = 0;

GetOptions(
   'infile_tax_contigs=s' => \$infile_tax_contigs,
   'infile_tax_genes=s'   => \$infile_tax_genes,
   'infile_abundance=s'   => \$infile_abundance,
   'infile_gff=s'         => \$infile_gff,
   'infile_blastn=s'      => \$infile_blastn,
   'infile_blastp=s'      => \$infile_blastp,
   'outfile_taxonomy=s'   => \$taxonomy_outfile,
   'verbose'              => \$verbose,
   'help'                 => \$help
);
if ($help) { print $usage; exit; }

## MAIN

# Parse Prodigal's gff output.
my $megahit_flag = 0;
my $megahit_id = "";
my %hash_gene_to_contig;
open(IN, "<".$infile_gff) or die "Can't open $infile_gff\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];

   # Hack if assembly is done with ray or megahit
   if($contig_id =~ m/^contig-(\d+)/){
      $contig_id = $1; 
   }elsif($contig_id =~ m/^(k\d+_)(\d+) /){
      $contig_id = $2;
      $megahit_flag = 1;
      $megahit_id = $1;
   }

   #1 to 7
   my $attributes = $row[8];

   my @attributes = split(/;/, $attributes);
   my $gene_id = $attributes[0];
   $hash_gene_to_contig{$gene_id} = $contig_id;
}
close(IN);
#print STDERR Dumper(\%hash_gene_to_contig);


my %hash_contigs;
open(IN, "<".$infile_tax_contigs) or die "Can't open $infile_tax_contigs\n";
while(<IN>){
   chomp;
   next if($. == 1);
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $accession = $row[1];

   my $k = "k__".$row[2].";";
   my $p = "p__".$row[3].";";
   my $c = "c__".$row[4].";";
   my $o = "o__".$row[5].";";
   my $f = "f__".$row[6].";";
   my $g = "g__".$row[7].";";
   my $s = "s__".$row[8]."";

   my $lineage = $k.$p.$c.$o.$f.$g.$s;
   $hash_contigs{$contig_id}{lineage} = $lineage;     
   $hash_contigs{$contig_id}{accession} = $accession;     
}
close(IN);
#print STDERR Dumper(\%hash_contigs);

my %hash_genes;
open(IN, "<".$infile_tax_genes) or die "Can't open $infile_tax_genes\n";
while(<IN>){
   chomp;
   next if($. == 1);
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   my $accession = $row[1];

   my $k = "k__".$row[2].";";
   my $p = "p__".$row[3].";";
   my $c = "c__".$row[4].";";
   my $o = "o__".$row[5].";";
   my $f = "f__".$row[6].";";
   my $g = "g__".$row[7].";";
   my $s = "s__".$row[8]."";

   my $lineage = $k.$p.$c.$o.$f.$g.$s;
   $hash_genes{$gene_id}{lineage} = $lineage;     
   $hash_genes{$gene_id}{accession} = $accession;     
}
close(IN);
#print STDERR Dumper(\%hash_genes);

#Parse blast_nt
my %hash_blastn;
open(IN, "<".$infile_blastn) or die "Can't open $infile_blastn\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $subject = $row[1];
   #$subject =~ m/accession\|(\d+)\|/;
   #my $accession = $subject;
   #$hash_accession{$accession} = ""; 
   $hash_blastn{$row[0]} = $row[10];
   #$hash_lengths{$row[0]} = $row[3];
}
close(IN);
#print STDERR Dumper(\%hash_blastn);

#Parse blast_nr
my %hash_blastp;
open(IN, "<".$infile_blastp) or die "Can't open $infile_blastp\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $subject = $row[1];
   #$subject =~ m/accession\|(\d+)\|/;
   #my $accession = $subject;
   #$hash_accession{$accession} = ""; 
   $hash_blastp{$row[0]} = $row[10];
   #$hash_lengths{$row[0]} = $row[3];
}
close(IN);
#print STDERR Dumper(\%hash_blastp);

open(IN, "<".$infile_abundance) or die "Can't open $infile_tax_genes\n";
open(OUT, ">".$taxonomy_outfile) or die "Can't open $taxonomy_outfile\n";


print OUT "gene_id\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
while(<IN>){
   chomp;
   if($. == 1){
      my @row = split(/\t/, $_);
      shift @row;
      my $new_row = join("\t", @row);
      print STDOUT "#OTU ID\t".$new_row."\ttaxonomy\n";
	  next;
   }
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   my $contig_id = $hash_gene_to_contig{$gene_id}; #it has to exist!

   my $tax_contig = "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL";
   my $accession_contig = "NA";
   my $accession_gene = "NA";
   if(exists $hash_contigs{$contig_id}){
      $tax_contig = $hash_contigs{$contig_id}{lineage};
      $accession_contig = $hash_contigs{$contig_id}{accession};
   }
   my $tax_gene = "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL";
   if(exists $hash_genes{$gene_id}){
      $tax_gene = $hash_genes{$gene_id}{lineage};
      $accession_gene = $hash_genes{$gene_id}{accession};
   }
   
   my $final_tax = "";
   my $accession = "";
 
   if($tax_gene eq "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL"     && $tax_contig eq "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL"){
      $final_tax = $tax_contig; # Does not matter which actually.
      $accession = $accession_contig;
   }elsif($tax_gene eq "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL" && $tax_contig ne "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL"){
      $final_tax = $tax_contig;
      $accession = $accession_contig;
   }elsif($tax_gene ne "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL" && $tax_contig eq "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL"){
      $final_tax = $tax_gene;    
      $accession = $accession_gene;
   }elsif($tax_gene ne "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL" && $tax_contig ne "k__NULL;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL"){
      
      if($hash_blastp{$gene_id} < $hash_blastn{$contig_id}){
		 $final_tax = $tax_gene;
         $accession = $accession_gene;
	  }elsif($hash_blastp{$gene_id} > $hash_blastn{$contig_id}){
         $final_tax = $tax_contig;
         $accession = $accession_contig;
	  }else{
	     $final_tax = $tax_contig;
         $accession = $accession_contig;
      }
   }
   print STDOUT $_."\t".$final_tax."\n";
   #my @final_tax = split(/;/, $final_tax);
   my $final_tax2 = $final_tax;
   #print STDERR $final_tax2."\n";
   $final_tax2 =~ s/k__/\t/;
   $final_tax2 =~ s/;p__/\t/;
   $final_tax2 =~ s/;c__/\t/;
   $final_tax2 =~ s/;o__/\t/;
   $final_tax2 =~ s/;f__/\t/;
   $final_tax2 =~ s/;g__/\t/;
   $final_tax2 =~ s/;s__//;

   print OUT $gene_id."\t".$accession.$final_tax2."\n";
}
close(IN);
close(OUT);

