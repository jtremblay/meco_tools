#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseUsearchPangenome.pl

PURPOSE:

INPUT:
--infile <string> : usearch .uc output file.
				
OUTPUT:

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
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

# do a first pass to identify cluster "C"
my %hash;
my %hash2;
my %hash_bins;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    my $seed = $row[0];
    my $gene_bin = $row[8];

    if($seed eq "S"){
        my @gene_bin = split(/;/, $gene_bin);
        my $gene = $gene_bin[0];
        my $bin = $gene_bin[1];
        $hash{$gene} = $gene;
        #$hash{$gene}{gene} = $gene;
        $hash_bins{$bin} = $bin;
        $hash2{$gene}{$bin} = 1;
    }
}
close(IN);

# Second pass where we want to match what genes belong to each cluster representative.
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    my $seed = $row[0];
    my $gene_bin_q = $row[8];
    my $gene_bin = $row[9];

    if($gene_bin ne "*"){
        my @gene_bin = split(/;/, $gene_bin);
        my $gene = $gene_bin[0];
        my $bin = $gene_bin[1];

        if(exists $hash{$gene}){
            my @gene_bin_q = split(/;/, $gene_bin_q);
            my $gene_q = $gene_bin_q[0];
            my $bin_q = $gene_bin_q[1];

            $hash2{$gene}{$bin_q} = 1;
        
        }
    }
}
close(IN);

print STDERR Dumper(\%hash2);
print STDERR Dumper(\%hash_bins);

# Then a 3rd pass to add zeros to absent bins.
for my $gene (sort keys %hash2){
    foreach my $bin (sort keys %{ $hash2{$gene} }) {
        
        for my $ref_bin (sort keys %hash_bins){
            if(!exists $hash2{$gene}{$ref_bin}){
                $hash2{$gene}{$ref_bin} = 0;
            }
        }
    }
}


print STDERR Dumper(\%hash2);

# Finally print in the form of a matrix.
print STDOUT "gene_id";
for my $ref_bin (sort keys %hash_bins){
    print STDOUT "\t".$ref_bin;
}
print STDOUT "\n";


for my $gene (sort keys %hash2){
    print STDOUT $gene;
    foreach my $bin (sort keys %{ $hash2{$gene} }) {
        print STDOUT "\t".$hash2{$gene}{$bin};
    }
    print STDOUT "\n";
}
