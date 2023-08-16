#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
use Cwd 'abs_path';
use Statistics::Descriptive qw(:all);

my $usage=<<'ENDHERE';
NAME:
summarizeMGData.pl

PURPOSE:
Only keep contigs that were binned and generate reduced load datasets from
annotations, abundance and checkm tables.

INPUT:
--annotations <string>       : metagenomics annotation master file (final file of the pipeline)
--abundance_contigs <string> : Contigs abundance.
--abundance_genes <string>   : Genes abundance.
--checkm <string>            : check outfile.
--parsed_bin_table <string>  : output of parseBins.pl	

OUTPUT:
--outdir <string>            : Directory where data will be generated. In this directory, there will
                               be 
                               - annotations (bin, contig, gene, taxonomy and kegg annotation)
                               - abundance of genes only present in contigs.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $annotations, $abundance_contigs, $abundance_genes, $checkm, $parsed_bin_table, $outdir, $contigs_indir);
my $verbose = 0;

GetOptions(
   'annotations=s'       => \$annotations,
   'abundance_contigs=s' => \$abundance_contigs,
   'abundance_genes=s'   => \$abundance_genes,
   'checkm=s'            => \$checkm,
   'parsed_bin_table=s'  => \$parsed_bin_table,
   'contigs_indir=s'     => \$contigs_indir,
   'outdir=s'            => \$outdir,
   'verbose' 	          => \$verbose,
   'help'                => \$help
);
if ($help) { print $usage; exit; }

die("--annotations missing\n") unless $annotations;
die("--abundance_contigs missing\n") unless $abundance_contigs;
die("--abundance_genes missing\n") unless $abundance_genes;
die("--parsed_bin_table missing\n") unless $parsed_bin_table;
die("--checkm missing\n") unless $checkm;

#print STDERR "outdir: $outdir\n";
#$outdir = abs_path($outdir);
#print STDERR "outdir: $outdir\n";
system("mkdir -p $outdir");

## MAIN
my %hash_bins_contigs;
my %hash_contigs_bins;

open(IN, "<".$parsed_bin_table) or die "Can't open $parsed_bin_table\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $bin_id = $row[0];
   my $contig_id = $row[1];

   $hash_bins_contigs{$bin_id}{CONTIG_ID} = $contig_id;
   $hash_contigs_bins{$contig_id} = $bin_id;
}
close(IN);

# Calculate GC from contigs

foreach my $key (keys %hash_bins_contigs){
   my $curr_bin = $key;
   print STDERR "Curr bin: $curr_bin\n";
   my $fasta = $contigs_indir."/".$curr_bin.".fa";

   my $stat = Statistics::Descriptive::Full->new();
   my @gcPerc;
   my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my $header = $curr->header;
      $header =~ s/>//;
      my $length = length($curr->seq);

      # Compute GC
      my $gcCount = 0;
      $gcCount += ($curr->seq =~ tr/gGcC/gGcC/);
      my $gcPerc = sprintf( "%0.2f", ($gcCount / $length) * 100);
      #push(@gcPerc, $gcPerc);
      $stat->add_data($gcPerc); 
   
   }
   my $mean = $stat->mean();
   #my $mean = mean( @gcPerc);
   $hash_bins_contigs{$key}{GC} = $mean;
}

print STDERR Dumper(\%hash_bins_contigs);
exit;

# Annotations
my %hash_genes_bins_contigs;
open(OUT, ">$outdir/parsed_annotations.tsv") or die "Can't open $outdir/bins_annotations.tsv";
print OUT "bin_id\tcontig_id\tgene_id\ttax_k\ttax_p\ttax_c\ttax_o\ttax_f\ttax_g\ttax_s\tkegg_pathway\tblastp\n";
open(IN, "<".$annotations) or die "Can't open $annotations\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $gene_id = $row[1];
   my $tax_k = $row[26];
   my $tax_p = $row[27];
   my $tax_c = $row[28];
   my $tax_o = $row[29];
   my $tax_f = $row[30];
   my $tax_g = $row[31];
   my $tax_s = $row[32];
   my $kegg_pathway = $row[4];
   my $blastp = $row[25];

   if(exists($hash_contigs_bins{$contig_id})){
      print OUT $hash_contigs_bins{$contig_id}."\t".$contig_id."\t".$gene_id."\t".$tax_k."\t".$tax_p."\t".$tax_c."\t".$tax_o."\t".$tax_f."\t".$tax_g."\t".$tax_s."\t".$kegg_pathway."\t".$blastp."\n";
      $hash_genes_bins_contigs{$gene_id}{bin} = $hash_contigs_bins{$contig_id};
      $hash_genes_bins_contigs{$gene_id}{contig} = $contig_id;
   }

}
close(IN);
close(OUT);

# Parse abundance files and only keep contigs having bins.
open(OUT, ">$outdir/parsed_abundance_contigs.tsv") or die "Can't open $outdir/bins_abundance_contigs.tsv";

# Annotations
open(IN, "<".$abundance_contigs) or die "Can't open $abundance_contigs\n";
while(<IN>){
   chomp;
   if($. == 1){
      print OUT $_."\n";
   }
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   if(exists($hash_contigs_bins{$contig_id})){
       #print OUT $hash_contigs_bins{$contig_id}."\n";
       print OUT $_."\n";
   }
}
close(IN);
close(OUT);

# Parse abundance files and only keep contigs having bins.
open(OUT, ">$outdir/parsed_abundance_genes.tsv") or die "Can't open $outdir/bins_abundance_genes.tsv";

open(IN, "<".$abundance_genes) or die "Can't open $abundance_genes\n";
while(<IN>){
   chomp;
   if($. == 1){
      print OUT $_."\n";
   }
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   if(exists($hash_genes_bins_contigs{$gene_id})){
       print OUT $_."\n";
   }
}
close(IN);
close(OUT);

exit;
