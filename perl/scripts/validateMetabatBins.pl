#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
validateMetabatBins.pl

PURPOSE:
Takes indirectory of metabin and generates a tsv file with each bins + percentages of each 
different order in each bins. If --split_bins_by_taxonomy <string> is set, only to up to
3 most abundant bins of that particular <string> will be kept.

INPUT:
--indir <string>            : indir containing fasta files.
--taxonomy <string>         : taxonomy.tsv file (contigs tax file.)
--split <string>            : Either 'k__', 'p__', 'c__', 'o__', 'f__', 'g__' or 's__' (default: 'o__')
--outdir <string>           : Outdir where to write new bins. 
--keep_most_abundant <int>  : Keep most n abundant bins according to taxonomy.

OUTPUT:
STDOUT
--outfile_link_raw <string> : File were will be written bin_id\tcontig_id\ttaxonomy\n

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $taxonomy, $split, $outdir, $n, $outfile_link_raw);
my $verbose = 0;

GetOptions(
   'indir=s'              => \$indir,
   'outdir=s'             => \$outdir,
   'outfile_link_raw=s'   => \$outfile_link_raw,
   'taxonomy=s'           => \$taxonomy,
   'split=s'              => \$split,
   'verbose'              => \$verbose,
   'keep_most_abundant=i' => \$n,
   'help'                 => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--indir is missing\n" unless($indir);
die "--outdir is missing\n" unless($outdir);
die "--taxonomy is missing\n" unless($taxonomy);

$n = 3 unless($n);

## MAIN
$indir = abs_path($indir);
$outdir = abs_path($outdir);
system("mkdir -p $outdir");

opendir(DIR, $indir);
my @files = grep(/\.fa$/,readdir(DIR));
closedir(DIR);

my @prefixes;
my @infiles;
foreach (@files){
   my $prefix = $_;
   $prefix =~ s{.*/}{};      # removes path  
   $prefix =~ s{\.[^.]+$}{}; # removes extension
   push(@prefixes ,$prefix);
   push(@infiles, $indir."/".$_);
}

my %hash_tax;
open(IN, "<".$taxonomy) or die "Can't open $taxonomy\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $species = join(".", $row[7]);
   my $lineage = join("\t", $row[1], $row[2], $row[3], $row[4], $row[5], $row[6], $species);
   
   $hash_tax{$contig_id} = $lineage;

}
close(IN);


# Parse fasta.
my %hash_raw;
for my $fasta (@infiles){
   my $prefix = shift(@prefixes);

   next if($prefix eq "all");

   my %hash;
   my %hash_fh;

   my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my ($contig_id) = $curr->header =~ /^>(\S+)/;
      $contig_id =~ s/^>//;
      
      if(exists $hash_tax{$contig_id}){

         my $lineage = $hash_tax{$contig_id};
         
         # Populate hash_raw
         $hash_raw{$prefix}{$contig_id} = $lineage;
         
         my @row = split(/\t/, $lineage);
         my $k = $row[0]; 
         my $p = $row[1]; 
         my $c = $row[2]; 
         my $o = $row[3]; 
         $o =~ s/ /_/g; # remove whitespaces
         my $f = $row[4]; 
         my $g = $row[5]; 
         my $s = $row[6];
         if(exists $hash_fh{$o}){
            $hash_fh{$o}{abundance}++;
            
         }else{
            $hash_fh{$o}{file} = ">>".$outdir."/".$prefix."-".$o.".fa";
            $hash_fh{$o}{abundance} = 1;

         }
      }
   }

   # Loop to keep only $n most abundant;
   my $k=1;
   my @keys = sort { $hash_fh{$b}{abundance} <=> $hash_fh{$a}{abundance} } keys(%hash_fh);
   foreach my $key ( @keys ) {
       #print STDERR $key."\t".$hash_fh{$key}{abundance}."\n";
      if($k <= $n){
         # means most it is in the most n abundants.
      }else{
         delete($hash_fh{$key});
      }

      $k++;
   }

   #print STDERR Dumper(\%hash_fh);
   #exit;
   
   $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my ($contig_id) = $curr->header =~ /^>(\S+)/;
      $contig_id =~ s/^>//;

      if(exists $hash_tax{$contig_id}){
         my $lineage = $hash_tax{$contig_id};
         my @row = split(/\t/, $lineage);
         my $k = $row[0]; 
         my $p = $row[1]; 
         my $c = $row[2]; 
         my $o = $row[3]; 
         $o =~ s/ /_/g; # remove whitespaces
         my $f = $row[4]; 
         my $g = $row[5]; 
         my $s = $row[6];

         if(exists $hash_fh{$o}){
            open(OUT, $hash_fh{$o}{file}) or die "Can't open $hash_fh{$o}\n";
            print OUT $curr->output;
            close(OUT);
            my $curr_j = $hash_fh{$o}{j};
            $hash{$prefix}{$o}{$contig_id} = $lineage;
         }
      }
   }
   #print STDERR Dumper(\%hash);
   for my $k1 (sort keys %hash) {
      for my $k2 (sort keys %{ $hash{$k1} }){
         for my $k3 (sort keys %{ $hash{$k1}{$k2} }){
            print STDOUT $k1."-".$k2."\t".$k3."\t".$hash{$k1}{$k2}{$k3}."\n";
         }
      }
   }
}

if(defined($outfile_link_raw)){
    open(OUT, ">".$outfile_link_raw) or die "Can't open $outfile_link_raw\n";
    for my $k1 (sort keys %hash_raw) {
       for my $k2 (sort keys %{ $hash_raw{$k1} }){
        print OUT $k1."\t".$k2."\t".$hash_raw{$k1}{$k2}."\n";
       
       }
    }
    close(OUT);
}
exit;
