#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
searchAnnotations.pl

PURPOSE:
To search annotations.tsv file generated the metagenomics pipeline.

INPUT:
--infile <string>    : Annotation file.
--column <int>       : column to search for term. (1st column = 0)
--term_file <string> : File containing terms to search for. (1 per line)
--case               : Boolean - if set will be case insensitive.
				
OUTPUT:
STDOUT               : Rows of annotation file matching terms in selected column.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $column, $termFile, $case);
my $verbose = 0;

GetOptions(
   'infile=s' 	   => \$infile,
   'column=i'     => \$column,
   'term_file=s'  => \$termFile,
   'case'         => \$case,
   'verbose' 	   => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## MAIN
die "--infile missing...\n" unless($infile);
die "--column missing...\n" unless($column);
die "--termFile missing...\n" unless($termFile);
#die "--case missing...\n" unless($case);

$case = 0 unless($case);

# Construct term string.
open(TERM, "<".$termFile) or die "Can't open $termFile\n";
my @terms;
while(<TERM>){
   chomp;
   push(@terms, $_);
}
my $terms = join("|", @terms);
$terms = "(".$terms.")";
print STDERR "Will search for $terms\n";
close(TERM);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   if($_ =~ m/^#/){
      print STDOUT $_."\n";
      print STDERR "... in column labeled $row[$column]\n";
      next;
   }

   if($case){
      if($row[$column] =~ m/$terms/i){
         print STDOUT "$_\n";
      }
   }else{
      if($row[$column] =~ m/$terms/){
         print STDOUT "$_\n";
      }
   }
}
exit;
