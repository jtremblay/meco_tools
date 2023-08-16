#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateTaxonomyFileForAnvio.pl

PURPOSE:

INPUT:
--infile_taxonomy    <string>   : taxonomic annotations file (shotgunMG pipeline).
--infile_gene_caller <string>   : infile gene callers (the one used to build anvio db).
                                  Mainly to associate gene ids to contig ids.

OUTPUT:
<STDOUT>                        : taxonomy file compatible with anvio file import.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxonomy, $infile_gene_calls);
my $verbose = 0;

GetOptions(
   'infile_taxonomy=s'     => \$infile_taxonomy,
   'infile_gene_calls=s'   => \$infile_gene_calls,
   'verbose'               => \$verbose,
   'help'                  => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

print STDOUT join("\t", ("gene_caller_id", "t_kingdom", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_specie")); print STDOUT "\n";
open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){
    next if($_ =~ m/^#/);
    chomp;
    my @row = split(/\t/, $_);
    my $contig_id = $row[0];
    my $lineage = $row[1]."\t".$row[2]."\t".$row[3]."\t".$row[4]."\t".$row[5]."\t".$row[6]."\t".$row[7];
    $hash{$contig_id} = $lineage;


}
close(IN);

open(IN, "<".$infile_gene_calls) or die "Can't open $infile_gene_calls\n";
while(<IN>){
    chomp;
    if($. == 1){next;}
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];
    $gene_id =~ s/gene_id_//;
    my $contig_id = $row[1];
    print STDOUT $gene_id."\t".$hash{$contig_id}."\n";
}
close(IN);
