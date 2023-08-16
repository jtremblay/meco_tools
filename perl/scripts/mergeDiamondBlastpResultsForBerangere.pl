#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
mergeDiamondBlastpResultsForBerangere.pl

PURPOSE:

INPUT:
--infile_diamond <string>     : Sequence file.
--infile_desc <string>        : Description of each fastq entry.
--infile_annotations <string> : annotations file provided by the pipeline.

OUTPUT:
<STDOUT>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_diamond, $infile_desc, $infile_annotations);
my $verbose = 0;

GetOptions(
   'infile_diamond=s' => \$infile_diamond,
   'infile_desc=s' => \$infile_desc,
   'infile_annotations=s' => \$infile_annotations,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;

open(IN, "<".$infile_desc) or die "Can't open $infile_desc\n";
while(<IN>){
   chomp;
   my @row = split(/ /, $_);
   my $id = shift(@row);
   my $annotation = join(" ", @row);
   $hash{$id} = $annotation;
}
close(IN);

my %hash_genes_contigs;
open(IN, "<".$infile_annotations) or die "Can't open $infile_annotations\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = shift(@row);
   my $gene_id = shift(@row);
   $hash_genes_contigs{$gene_id} = $contig_id."\t".join("\t", @row)."";
}
close(IN);

#print STDERR Dumper(\%hash);

open(IN, "<".$infile_diamond) or die "Can't open $infile_diamond\n";
while(<IN>){
    chomp;
    my @row = split(/\s/, $_);
    my $id = $row[1];
    my $gene_id = $row[0];
    if(exists $hash{$id}){
        print STDOUT $_."\t".$hash{$id}."\t".$hash_genes_contigs{$gene_id}."\n";
    }else{
        print STDOUT $_."\t"."Undefined hit"."\n";
    }
}
close(IN);



