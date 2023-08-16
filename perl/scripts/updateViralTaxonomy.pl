#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
updateViralTaxonomy.pl

PURPOSE:

INPUT:
--infile <string>        : matrix for which the first col = feature id to match with the ones into viral contigs file.
--viral_contigs <string> : Contigs files. One entry of viral contig header per line.
				
OUTPUT:
STDOUT                   : updated table.
STDERR                   : Some statistics.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $viral_contigs);
my $verbose = 0;

GetOptions(
   'infile=s'       => \$infile,
   'viral_contigs=s' => \$viral_contigs,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
open(IN, "<".$viral_contigs) or die "Can't open $viral_contigs\n";
while(<IN>){
   chomp;
   #my @row = split(/\t/, $_);
   $hash{$_} = $_;
}
close(IN);

#print STDERR Dumper(\%hash);
#exit;

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#/){
        print STDOUT $_."\n";
        next;
    }
    my @row = split(/\t/, $_);
    my $contig_id = $row[0];

    if(exists $hash{$contig_id}){
        my $taxonomy = pop(@row);
        my @taxa = split(/;/, $taxonomy);
        $taxa[0] = "k__Viruses";
        my $new_taxonomy = join(";", @taxa);

        print STDOUT join("\t", @row)."\t".$new_taxonomy."\n";
        print STDERR $new_taxonomy."\t".$taxonomy."\n";
   }
}
close(IN);
