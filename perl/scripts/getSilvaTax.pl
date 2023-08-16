#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
getSilvaTax.pl

PURPOSE:

INPUT:
--infile_blastn <string>    : Sequence file
--infile_accession <string> : Accession file
--infile_abundance <string> : Abundance file
--evalue <float>            : evalue cutoff
--length <length>           : length cutoff

OUTPUT:
<STDOUT>                    : OTU table format abundance+taxa file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_blastn, $infile_accession, $infile_abundance, $evalue, $length);
my $verbose = 0;

GetOptions(
   'infile_blastn=s'    => \$infile_blastn,
   'infile_accession=s' => \$infile_accession,
   'infile_abundance=s' => \$infile_abundance,
   'evalue=f'           => \$evalue,
   'length=i'           => \$length,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

$evalue = 1e-10 unless($evalue);
$length = 90 unless($length);

## MAIN
my %hash_tax;
open(IN, "<".$infile_accession) or die "Can't open $infile_accession\n";
while(<IN>){
   chomp;
   #my @row = split(/\t/, $_);
   my $silva_id = "";
   my $taxonomy = "";

   if($_ =~ m/^>(\S+) (.*)$/){
      $silva_id = $1;
      $taxonomy = $2;
      $silva_id =~ s/^>//;
      $hash_tax{$silva_id} = $taxonomy;
   }else{
      die "Can't parse $_\n";
   }
    
}
close(IN);

my %hash;
open(IN, "<".$infile_blastn) or die "Can't open $infile_blastn\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $silva_id = $row[1];
   my $evalue = $row[10];
   my $length = $row[3];
    
   #get taxonomy
   my $taxonomy = $hash_tax{$silva_id};

   if(!exists $hash{$contig_id}){
       $hash{$contig_id}{silva_id} = $silva_id;
       $hash{$contig_id}{taxonomy} = $taxonomy;
       $hash{$contig_id}{evalue} = $evalue;
       $hash{$contig_id}{length} = $length;
   }
}
close(IN);

open(IN, "<".$infile_abundance) or die "Can't open $infile_abundance\n";
while(<IN>){
   chomp;
   if($. == 1){
     print STDOUT $_."\ttaxonomy\n";
   }else{

       my @row = split(/\t/, $_);
       my $contig_id = shift(@row);
       
       my $taxonomy;
       if(exists $hash{$contig_id}){
          if($hash{$contig_id}{evalue} <= $evalue && $hash{$contig_id}{length} >= $length){
             $taxonomy = $hash{$contig_id}{taxonomy};
          }else{
             $taxonomy = "undef";
          }
       }else{
          $taxonomy = "undef";
       }

       print STDOUT $contig_id."\t".join("\t", @row)."\t".$taxonomy."\n";
   }
}
close(IN);


