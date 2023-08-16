#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum);

my $usage=<<'ENDHERE';
NAME:
filterCanopy.pl

PURPOSE:

INPUT:
--infile_clusters <string>    : canopy clusters
--infile_profiles <string>    : canopy profiles
--min_number_of_genes <int>   : min number of genes. Canopy having lower genes number than this value will be discarded
--high_abundance_perc <float> : percentage of high abundance profiles
--high_abundance_freq <int>   : frequencies of high abundance profiles
				
OUTPUT:
--outfile_clusters <string>    : canopy clusters filtered
--outfile_profiles <string>    : canopy profiles filtered

NOTES:
--high_abundance_perc 0.90 --high_abundance_freq 3 means that canopies having 3 samples representing 90% of abundance signal
will be discarded. Note that canopies having 1 or 2 samples representing 90% of abundance signal will be discarded as well.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_clusters, $infile_profiles, $outfile_clusters, $outfile_profiles, $min_number_of_genes, $high_abundance_perc, $high_abundance_freq, $outfile_abundance);
my $verbose = 0;

GetOptions(
   'infile_clusters=s'           => \$infile_clusters,
   'infile_profiles=s'           => \$infile_profiles,
   'outfile_clusters=s'          => \$outfile_clusters,
   'outfile_profiles=s'          => \$outfile_profiles,
   'outfile_canopy_abundance=s'  => \$outfile_abundance,
   'min_number_of_genes=i'       => \$min_number_of_genes,
   'high_abundance_perc=f'       => \$high_abundance_perc,
   'high_abundance_freq=i'       => \$high_abundance_freq,
   'verbose' 	                  => \$verbose,
   'help'                        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT_CLUSTERS, ">".$outfile_clusters) or die "Can't open $outfile_clusters\n";
open(OUT_PROFILES, ">".$outfile_profiles) or die "Can't open $outfile_profiles\n";
open(OUT_ABUNDANCE, ">".$outfile_abundance) or die "Can't open $outfile_abundance\n";

my %hash_gene_count;
my $pre_qc_canopies = 0;
my $pre_qc_genes = 0;

open(IN, "<".$infile_clusters) or die "Can't open $infile_clusters\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $canopy_id = $row[0];
   my $read = $row[1];

   if(exists $hash_gene_count{$canopy_id}){
      $hash_gene_count{$canopy_id}++;
   }else{
      $hash_gene_count{$canopy_id} = 0;
   }
   $pre_qc_genes++;
}
close(IN);

open(IN, "<".$infile_profiles) or die "Can't open $infile_profiles\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $canopy_id = shift(@row);
   $pre_qc_canopies++;
   next if($hash_gene_count{$canopy_id} < $min_number_of_genes);

   my @sorted_row = sort {$b <=> $a} @row;
   my $sum = sum(@sorted_row);
   my $curr_sum = 0;
   
   print STDERR join(",", @sorted_row)."\n" if($verbose);

   my $higher = 0;
   print STDERR "sum :".$sum."\n" if($verbose);
   for(my $i=0; $i<$high_abundance_freq; $i++){
      $curr_sum = $curr_sum + $sorted_row[$i];
      my $curr_perc = $curr_sum / $sum;
      print STDERR "Curr perc: ".$curr_perc."\n" if($verbose);
      if($curr_perc >= $high_abundance_perc){
         $higher = 1;
      }
   }  

   if($higher == 0){
      print OUT_PROFILES $_."\n";
   }
}
close(OUT_PROFILES);
close(IN);


## Finally, print qc passed clusters.
my %hash_qc_passed;
my $canopyCounts = 0;
open(IN, "<".$outfile_profiles) or die "Can't open $outfile_profiles\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $canopy_id = $row[0];

   $hash_qc_passed{$canopy_id} = 0;
   $canopyCounts++;
}
close(IN);

my $geneCounts = 0;
open(IN, "<".$infile_clusters) or die "Can't open $infile_clusters\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $canopy_id = $row[0];
   my $read = $row[1];

   if(exists $hash_qc_passed{$canopy_id}){
       print OUT_CLUSTERS "$canopy_id\t$read\n";
       $geneCounts++;
       $hash_qc_passed{$canopy_id}++;
   }
}
close(IN);
close(OUT_CLUSTERS);

print STDOUT "Number of initial (pre-qc) canopies: ".commify($pre_qc_canopies)."\n";
print STDOUT "Representing a total of # genes: ".commify($pre_qc_genes)."\n";
print STDOUT "Number of qc passed canopies: ".commify($canopyCounts)."\n";
print STDOUT "Representing a total of # genes: ".commify($geneCounts)."\n";
print STDOUT "Number of discarded canopies: ".commify($pre_qc_canopies - $canopyCounts)."\n";
print STDOUT "Number of discarded genes: ".commify($pre_qc_genes - $geneCounts)."\n";

my @keys = sort { $a cmp $b } keys %hash_qc_passed;
for my $key (@keys){
   print OUT_ABUNDANCE $key."\t".commify($hash_qc_passed{$key})."\n"; 
}
close(OUT_ABUNDANCE);

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

exit;
