#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile_taxonomy <string> : Sequence file
--infile_link <string>     : link tsv file

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxonomy, $infile_link);
my $verbose = 0;

GetOptions(
   'infile_taxonomy=s' => \$infile_taxonomy,
   'infile_link=s'     => \$infile_link,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
open(IN, "<".$infile_link) or die "Can't open $infile_link\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];
    my $contig_id = $row[1];

    $hash{$gene_id} = $contig_id;
}
close(IN);

open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];

    if($. == 1){
        print STDOUT "contig_id\t$_\n";
    }

    if(exists($hash{$gene_id})){
        print STDOUT $hash{$gene_id}."\t".$_."\n";    
    }
}
close(IN);
