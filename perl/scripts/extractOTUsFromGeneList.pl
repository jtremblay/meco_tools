#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
extractOTUsFromGeneList.pl

PURPOSE:

INPUT:
--infile <string>    : OTU table
--gene_list <string> : gene list
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $gene_list);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'gene_list=s' => \$gene_list,
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$gene_list) or die "Can't open $gene_list\n";
while(<IN>){
   chomp;
   $hash{$_} = "";
}
close(IN);

print STDERR Dumper(\%hash);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   if($_ =~ m/^#/){
       print STDOUT $_."\n";
   }else{
      my @row = split(/\t/, $_);
      if(exists $hash{$row[0]}){
         print STDOUT $_."\n";
      }
   }
}
close(IN);
