#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
convertGffToAnvio.pl

PURPOSE:

INPUT:
--infile <string> : Called genes in gff
				
OUTPUT:
STDOUT <string>   : Called genes in tsv for anvio

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
print STDOUT "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n";
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);

   my ($gene_id) = $_ =~ m/gene_id=(\d+)/;
   $gene_id =~ s/=/_/;

   my @row = split(/\t/, $_);
   my ($contig_id) = $row[0] =~ m/^(\S+)/;
   my $start = $row[3];
   my $end = $row[4];
   my $start1 = $start - 1; # To make compatible with Anvio which calculate sequence length differently than Metagenemark
   
   my $direction = "";
   if($end > $start){
      $direction = "f";
   }elsif($start > $end){
      $direction = "r";
   }else{
      die "Start ($start) and End ($end) are of equal values! - This should not happen, something is wrong in the gene calling process.\n";
   }

   my $partial = "1"; # 0 for complete call and 1 for incomplete call.

   my $source = "metagenemark";
   my $version = "v1.0";

   print STDOUT "$gene_id\t$contig_id\t$start1\t$end\t$direction\t$partial\t$source\t$version\n";
}
close(IN);
