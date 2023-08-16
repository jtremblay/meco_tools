#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateKeggFileForAnvio.pl

PURPOSE:

INPUT:
--infile_kegg_annotations <string>   : kegg annotations file (shotgunMG pipeline).
--infile_kegg_diamond_blastp <string : kegg diamond blastp results.

OUTPUT:
<STDOUT>                             : function file compatible with anvio file import.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_kegg_annotations, $infile_kegg_diamond_blastp);
my $verbose = 0;

GetOptions(
   'infile_kegg_annotations=s'    => \$infile_kegg_annotations,
   'infile_kegg_diamond_blastp=s' => \$infile_kegg_diamond_blastp,
   'verbose'                      => \$verbose,
   'help'                         => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile_kegg_diamond_blastp) or die "Can't open $infile_kegg_diamond_blastp\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];
    $gene_id =~ s/gene_id_//;
    my $evalue = $row[10];
    $hash{$gene_id} = $evalue;
}
close(IN);

print STDOUT join("\t", ("gene_callers_id", "source", "accession", "function", "e_value")); print STDOUT "\n";
open(IN, "<".$infile_kegg_annotations) or die "Can't open $infile_kegg_annotations\n";
while(<IN>){
    next if($_ =~ m/^#/);
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];
    $gene_id =~ s/gene_id_//;
    my $kegg_gene = $row[1];
    my $KO = $row[2];
    my $name = $row[3];
    if($name eq $row[2]){
        $name = "-";
    }
    my $definition = $row[5];

    print STDOUT $gene_id."\tKEGG_KO\t".$KO."\t".$definition."\t1e-05\n";

    my $pathways = $row[6];
    my $pathways_desc = $row[7];
    my $modules = $row[8];
    my $modules_desc = $row[9];
    
    my @pathways = split(/==/, $pathways);
    my @pathways_desc = split(/==/, $pathways_desc);
    for(my $i=0; $i<@pathways;$i++){
        print STDOUT $gene_id."\tKEGG_pathway\t".$pathways[$i]."\t".$pathways_desc[$i]."\t".$hash{$gene_id}."\n";
    }
    
    my @modules = split(/==/, $modules);
    my @modules_desc = split(/==/, $modules_desc);
    
    for(my $i=0; $i<@modules;$i++){
        print STDOUT $gene_id."\tKEGG_module\t".$modules[$i]."\t".$modules_desc[$i]."\t".$hash{$gene_id}."\n";
    }
}
close(IN);
