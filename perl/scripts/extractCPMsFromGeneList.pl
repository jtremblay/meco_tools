#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
extractCPMsFromGeneList.pl

PURPOSE:

INPUT:
--infile <string>    : cpm matrix (or other abundance matrix)
--gene_list <string> : txt file having one gene_id_12345 entry per line
--invert             : set flag if to exclude entries instead of keeping them.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $gene_list, $invert, $number_of_header_lines);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'gene_list=s' => \$gene_list,
   'invert'      => \$invert,
   'number_of_header_lines=s' => \$number_of_header_lines,
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN

$number_of_header_lines = 1 unless($number_of_header_lines);

my %hash;
open(IN, "<".$gene_list) or die "Can't open $gene_list\n";
while(<IN>){
   chomp;
   my $entry = $_;
   $entry =~ s/\t//;
   $hash{$entry} = $entry;
}
close(IN);


if($invert){
    open(IN, "<".$infile) or die "Can't open $infile\n";
    while(<IN>){
       chomp;
       my @row = split(/\t/, $_);
       my $gene_id = $row[0];

       #if($. == 1){
       if($. <= $number_of_header_lines){
          print STDOUT $_."\n";
       }else{
           if(!exists $hash{$gene_id}){
              print STDOUT $_."\n";
           }
       }
    }
    close(IN);

}else{
    open(IN, "<".$infile) or die "Can't open $infile\n";
    while(<IN>){
       chomp;
       my @row = split(/\t/, $_);
       my $gene_id = $row[0];

       if($. <= $number_of_header_lines){
          print STDOUT $_."\n";
       }else{
           if(exists $hash{$gene_id}){
              print STDOUT $_."\n";
           }
       }
    }
    close(IN);
}
