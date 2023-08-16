#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $binid = $row[0];
    my $stats = $row[1];
    my @row2 = split(/,/, $stats);
    if($. == 1){
        print STDOUT "bin_id\ttranslation_table\tGC_std\tno_ambiguous_bases\tgenome_size\tlongest_contig\tN50_scaffolds\tmean_scaffold_length\tcontigs\tscaffolds\tpredicted_genes\tlongest_scaffold\tGC\tN50_contigs\tcoding_density\tmean_contig_length\n";
        print STDOUT "$binid";
        foreach my $el(@row2){
            my @row3 = split(/:/, $el);
            print STDOUT "\t".$row3[1];
        }
        print STDOUT "\n";
    }else{
        print STDOUT "$binid";
        foreach my $el(@row2){
            my @row3 = split(/:/, $el);
            print STDOUT "\t".$row3[1];
        }
        print STDOUT "\n";
    }
}
