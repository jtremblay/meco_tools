#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use Cwd 'abs_path';
use Statistics::Descriptive qw(:all);

my $usage=<<'ENDHERE';
NAME:
addReadsCountPerBinsCheckM.pl

PURPOSE:

INPUT:
--checkm <string>             : CheckM (tab table) or CheckM2 (quality_report.tsv) output.
--parsed_bins <string>        : Parsed bins output. To get contig id.
--contigs_abundance <string>  : Contigs abundance per sample matrix obtained right before metabat step.
--contigs_indir <string>      : in directory where are the binned fastas.
				
OUTPUT:
STDOUT                        : CheckM/CheckM2 output with number of reads per bins + estimated cov.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $checkm, $contigs_abundance, $parsed_bins, $contigs_indir);
my $verbose = 0;

GetOptions(
   'checkm=s' 	          => \$checkm,
   'parsed_bins=s' 	      => \$parsed_bins,
   'contigs_indir=s' 	  => \$contigs_indir,
   'contigs_abundance=s'  => \$contigs_abundance,
   'verbose'              => \$verbose,
   'help'                 => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

# CheckM or CheckM2 quality table. The header is taken from the file itself, so any column
# layout works: CheckM's tab table (first column "Bin Id") or CheckM2's quality_report.tsv
# (first column "Name"). Column names are made shell/R friendly: spaces become underscores
# and a leading "#" becomes "no" ("# genomes" -> "no_genomes").
my $checkm_header = "";
open(IN, "<".$checkm) or die "Can't open $checkm\n";
while(<IN>){
   chomp;
   s/\r$//;
   next if($_ =~ m/^\s*$/);

   my @row = split(/\t/, $_);
   if($checkm_header eq ""){
       shift(@row);
       foreach my $col (@row){
           $col =~ s/^#\s*/no_/;
           $col =~ s/\s+/_/g;
       }
       $checkm_header = join("\t", "Bin_Id", @row);
       next;
   }

   # CheckM2 reports bins without their file extension; strip one anyway in case it is kept.
   my $bin_id = shift(@row);
   $bin_id =~ s/\.(fa|fna|fasta)$//;
   $hash{$bin_id} = join("\t", @row);
}
close(IN);
#print STDERR Dumper(\%hash);


# Parsed bins
my %contigs_to_bin;
my %bin_to_contigs;
my %tax;
open(IN, "<".$parsed_bins) or die "Can't open $parsed_bins\n";
while(<IN>){
   chomp;

   my @row = split(/\t/, $_);
   my $bin_id = $row[0];
   my $contig_id = $row[1];
   shift(@row);
   shift(@row);
   my $lineage = join("\t", @row);

   # Basically, here do a survey of the number different lineages for each contig of a given bin.
   # These lineages can old some NULL values that will be used later.
   $contigs_to_bin{$contig_id} = $bin_id;
   if(exists $bin_to_contigs{$bin_id}){
      $bin_to_contigs{$bin_id} .= ":".$contig_id;
      $tax{$bin_id}{$lineage}++;
   }else{
      $bin_to_contigs{$bin_id} = $contig_id;
      $tax{$bin_id}{$lineage} = 1;
   }
}
close(IN);
#print STDERR Dumper(\%tax);
#print STDERR Dumper(\%bin_to_contigs);

# abundance
my %contigs_abundance;
open(IN, "<".$contigs_abundance) or die "Can't open $contigs_abundance\n";
while(<IN>){
   chomp;

   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $length = $row[1];
   my $coverage = $row[2];
   $contigs_abundance{$contig_id}{LENGTH} = $length;
   $contigs_abundance{$contig_id}{COV} = $coverage;
}
close(IN);
#print STDERR Dumper(\%contigs_abundance);

my %final_hash;
# for each bins, sum up the contig length + contig cov.
for my $bin_id (keys %hash){
   if(exists $bin_to_contigs{$bin_id}){
      my @contigs = split(/:/, $bin_to_contigs{$bin_id});
      #print STDERR join(":", @contigs)."\n";
      foreach my $contig_id (@contigs){
         if(exists $final_hash{$bin_id}){
            $final_hash{$bin_id}{LENGTH} += $contigs_abundance{$contig_id}{LENGTH};
            $final_hash{$bin_id}{COV} += $contigs_abundance{$contig_id}{COV};
            $final_hash{$bin_id}{NO} += 1
         }else{
            $final_hash{$bin_id}{LENGTH} = $contigs_abundance{$contig_id}{LENGTH};
            $final_hash{$bin_id}{COV} = $contigs_abundance{$contig_id}{COV};
            $final_hash{$bin_id}{NO} = 1
         }
      }
   }else{
       print STDERR "Bin ".$bin_id." does not exist\n";
   }
}

# Then calculate GC%
foreach my $bin_id (keys %final_hash){
   my $fasta = $contigs_indir."/".$bin_id.".fa";

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
      $stat->add_data($gcPerc); 
   
   }
   my $mean = $stat->mean();
   $final_hash{$bin_id}{GC} = $mean;
}

# Sort tax
my %tax2;
foreach my $bin_id (keys %tax) {
   my $max = 0;
   my $min = 999999;
   # Sorted so that ties between lineages are resolved the same way on every run.
   foreach my $lineage (sort keys %{ $tax{$bin_id} }) {
   
      my @lineage = split(/\t/, $lineage);

      ## Here skip taxonomy if the majority is set to NULL.
      my $x=0;
      foreach my $taxon (@lineage){
         $x++ if($taxon =~ m/NULL/);
      }

      # If we encounter 3 or more NULLs, reset the popularity (occurence) count to one.
      if($x >= 3){
         $tax{$bin_id}{$lineage} = 1;
      }

      ## If not NULL; (well not enough NULLs)...
      if($tax{$bin_id}{$lineage} > $max){
         $tax2{$bin_id} = $lineage;
         $max = $tax{$bin_id}{$lineage};
      }
   }
}
#print STDERR "Tax2 hash\n";
#print STDERR Dumper(\%tax2);
#print STDERR Dumper(\%final_hash);


for my $bin_id (keys %final_hash){
   $final_hash{$bin_id}{COV_AVG} = sprintf("%.3f" , $final_hash{$bin_id}{COV} / $final_hash{$bin_id}{NO})
}

print STDOUT $checkm_header."\tnumber_of_contigs\ttotal_contig_length\tavg_coverage\tGC%\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
my @bins = sort { $final_hash{$b}{LENGTH} <=> $final_hash{$a}{LENGTH} } keys %final_hash;
for my $bin_id (@bins){
   print STDOUT $bin_id."\t".$hash{$bin_id}."\t".$final_hash{$bin_id}{NO}."\t".$final_hash{$bin_id}{LENGTH}."\t".$final_hash{$bin_id}{COV_AVG}."\t".sprintf("%.3f", $final_hash{$bin_id}{GC})."\t".$tax2{$bin_id}."\n";
}
exit;
