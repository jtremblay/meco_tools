#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
generateRpoBMatrix.pl

PURPOSE:

INPUT:
--infile <string>    : cpm matrix (or other abundance matrix)
--cog_term <string>  : A COG ID. For rpoB, COG = COG0085
--cog_file <string>  : COG annotation file

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $cog_file, $cog_term, $number_of_header_lines);
my $verbose = 0;

GetOptions(
   'infile=s'                 => \$infile,
   'cog_file=s'               => \$cog_file,
   'cog_term=s'               => \$cog_term,
   'number_of_header_lines=s' => \$number_of_header_lines,
   'verbose'                  => \$verbose,
   'help'                     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

$number_of_header_lines = 1 unless($number_of_header_lines);

my %hash;
open(IN, "<".$cog_file) or die "Can't open $cog_file\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];
    if($row[12] =~ m/$cog_term/){
        $hash{$gene_id} = $row[12];
    }
}
close(IN);


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
