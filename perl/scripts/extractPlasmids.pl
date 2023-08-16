#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
extractPlasmids.pl

PURPOSE:

INPUT:
--infile_blast <string>      : Blast file of best plasmids hits.
--infile_fna <string>        : Fasta file of contigs.
--infile_faa <string>        : Fasta (aa) file of predicted genes.
--infile_pfam  <string>      : Hmmscan of genes against pfam db.
--infile_gff <string>        : Output of metagenemark (or prodigal).
				
OUTPUT:
--outfile_tsv <string>       : Tsv file of validated plasmids. To be integrated in final GFF file later.
--outfile_faa <string>       : outfile of each potential plasmid encoded genes in faa format.
--outfile_gene_list <string> : Simple tsv file having 2 columns. contigs-ID\tgene-ID
STDOUT                       : fasta file of potential plasmids contigs.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_blast, $infile_fna, $infile_faa, $infile_pfam, $infile_gff, 
    $outfile_tsv, $outfile_faa, $outfile_gene_list);
my $verbose = 0;

GetOptions(
   'infile_blast=s' 	      => \$infile_blast,
   'infile_fna=s' 	      => \$infile_fna,
   'infile_faa=s' 	      => \$infile_faa,
   'infile_pfam=s' 	      => \$infile_pfam,
   'infile_gff=s' 	      => \$infile_gff,
   'outfile_tsv=s' 	      => \$outfile_tsv,
   'outfile_faa=s' 	      => \$outfile_faa,
   'outfile_gene_list=s'   => \$outfile_gene_list,
   'verbose' 	            => \$verbose,
   'help'                  => \$help
);
if ($help) { print $usage; exit; }

## MAIN

die("--infile_blast <string> required\n") unless($infile_blast);
die("--infile_fna <string> required\n") unless($infile_fna);
die("--infile_faa <string> required\n") unless($infile_faa);
die("--infile_pfam <string> required\n") unless($infile_pfam);
die("--infile_gff <string> required\n") unless($infile_gff);
die("--outfile_tsv <string> required\n") unless($outfile_tsv);
die("--outfile_faa <string> required\n") unless($outfile_faa);
die("--outfile_gene_list <string> required\n") unless($outfile_gene_list);

# Open output in tsv
open(OUT, ">".$outfile_tsv) or die "Can't open $outfile_tsv\n";
open(OUT_FAA, ">".$outfile_faa) or die "Can't open $outfile_faa\n";
open(OUT_LIST, ">".$outfile_gene_list) or die "Can't open $outfile_gene_list\n";

my %hashPfam;
my %hash_gene_to_contig;

# 1- loop through gene prediciton to associate each gene_id with contig_id.
open(IN, "<".$infile_gff) or die "Can't open $infile_gff\n";
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    next if($_ =~ m/^$/);
    my($contigId) = $_ =~ m/^(contig-\d+)/;
    my($geneNum) = $_ =~ m/gene_id\=(\d+)/;
    my $geneId = "gene_id_".$geneNum;
    $hash_gene_to_contig{$geneId} = $contigId;
}
close(IN);
print STDERR "Contigs to genes linking done...\n";

# 2- loop through Pfam Hmmscan results to associate each gene_id a pfam id.
open(IN, "<".$infile_pfam) or die "Can't open $infile_pfam\n";
my %seen;
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\s+/, $_);
   my $gene_id = $row[2];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   $gene_id = "gene_id_".$gene_id;
   if(!exists $seen{$gene_id}){
      # get contig id first...
      my $contig_id = $hash_gene_to_contig{$gene_id};

      $hashPfam{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME} = $row[0];
      $hashPfam{$contig_id}{$gene_id}{PFAM_ACCESS} = $row[1];

      # Then get desc.
      my @desc = @row[18 .. (@row - 1)];
      my $desc = join(" ", @desc);
      $hashPfam{$contig_id}{$gene_id}{PFAM_DESC} = $desc;

      $seen{$gene_id} = 1;
   }
}
print STDERR "Pfam parsing done...\n";

# 3- Loop through blast best hits and only keep hits having:
# ID% > 0.95; hit lenth > 89 nt; query hit length at least 90% of the read length.
my %hash;
open(IN, "<".$infile_blast) or die "Can't open $infile_blast\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $queryId = $row[0];
   my $targetId = $row[1];
   my $percId = $row[2];
   my $hitLength = $row[3];
   my $evalue = $row[10];

   # If curr hit meeting filtering criteria
   if($percId >= 90.00 && $hitLength >= 89){
      $hash{$queryId}{targetId} = $targetId;
      $hash{$queryId}{percId} = $percId;
      $hash{$queryId}{hitLength} = $hitLength;
      $hash{$queryId}{evalue} = $evalue;
   }
}
close(IN);
print STDERR "First blast hits filtering done...\n";

# 4- Check if last criteria "query hit length at least 90% of the read length" is OK.
# If not, delete from hash.
my %hash_gene_to_contig_found;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fna) or die("Unable to open Fasta file, $infile_fna\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my($header) = $curr->header =~ m/^>(\S+-\d+) .*/;

   if(exists $hash{$header}){
      my $seq = $curr->seq;
      my $seqLength = length($seq);
      my $fraction = $hash{$header}{hitLength} / $seqLength;

      # If fraction > 0.90, store seq and gene_id in hash.
      if($fraction < 0.90){
         delete $hash{$header};
      }else{
   
         # If meets all criteria, then print
         print STDOUT $curr->output;

         #print STDERR "fraction >= .90\n";
         $hash{$header}{seq} = $seq;

         foreach my $geneId (keys %{ $hashPfam{$header} }) {
            print OUT "$header\t$hash{$header}{targetId}\t$hash{$header}{percId}\t$hash{$header}{hitLength}\t$hash{$header}{evalue}\t$hashPfam{$header}{$geneId}{PFAM_PRODUCT_NAME}\t$hashPfam{$header}{$geneId}{PFAM_ACCESS}\t$hashPfam{$header}{$geneId}{PFAM_DESC}\n";
            $hash_gene_to_contig_found{$geneId} = $header
         }

      }
   }   
}
close(OUT);

print OUT_LIST "#contig_ID\tgene_ID\n";
$ref_fasta_db = Iterator::FastaDb->new($infile_faa) or die("Unable to open Fasta file, $infile_faa\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my($geneId) = $curr->header =~ m/^>(gene_id_\d+)/;
   if(exists $hash_gene_to_contig_found{$geneId}){
      print OUT_FAA ">".$geneId."===".$hash_gene_to_contig{$geneId}."\n".$curr->seq."\n";
      print OUT_LIST $hash_gene_to_contig{$geneId}."\t".$geneId."\n";
   }
}
close(OUT_FAA);
close(OUT_LIST);

