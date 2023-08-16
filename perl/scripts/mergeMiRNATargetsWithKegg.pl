#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
mergeMiRNATargetsWithKegg.pl

PURPOSE:

INPUT:
--infile <string>          : mirna table
--kegg_functions <string>  : parse kegg table.
				
OUTPUT:
STDOUT <string>            : 

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $kegg_functions);
my $verbose = 0;

GetOptions(
   'infile=s'         => \$infile,
   'kegg_functions=s' => \$kegg_functions,
   'verbose'          => \$verbose,
   'help'             => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
open(IN, "<".$kegg_functions) or die "Can't open $kegg_functions\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $key = shift(@row);
   my $new_line = join("\t", @row);
   $hash{$key} = $new_line;
}
close(IN);

print STDERR Dumper(\%hash);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $key = $row[1]."_".$row[15]."_".$row[16]."_1";
   #my @new_row = @row;
   #print STDERR $key."\n";
   #shift(@new_row);
   if(exists $hash{$key}){
      print STDOUT join("\t", @row)."\t".$hash{$key}."\n";
   }
}
close(IN);
