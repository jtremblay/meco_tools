#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile_ko <string>          : KEGG list K12345\nK12333\netc. 
--infile_ko_to_genes <string> : KEGG link file KO\tgene_id
--infile_genes <string>       : KEGG faa genes fasta file
				
OUTPUT:
--outdir <string>             : outdir where will be written one multi fasta file per KO.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_ko, $infile_ko_to_genes, $infile_genes, $outdir);
my $verbose = 0;

GetOptions(
   'infile_ko=s'          => \$infile_ko,
   'infile_ko_to_genes=s' => \$infile_ko_to_genes,
   'infile_genes=s'       => \$infile_genes,
   'outdir=s'             => \$outdir,
   'verbose'              => \$verbose,
   'help'                 => \$help
);
if ($help) { print $usage; exit; }

my %hash;
open(IN, "<".$infile_ko) or die "Can't open $infile_ko\n";
while(<IN>){
    chomp;
    if($_ =~ m/(K\d\d\d\d\d)/){
        $hash{$1}{KO} = $1;
    }
}
close(IN);
print STDERR "Done parsing KO list\n";

my %hash_ref;
open(IN, "<".$infile_ko_to_genes) or die "Can't open $infile_ko_to_genes\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $ko = $row[0];
    my $kegg_gene = $row[1];
    $ko =~ s/ko://;
    if(exists($hash{$ko})){
        $hash_ref{$kegg_gene} = $ko;
    }
}
close(IN);
print STDERR "Done parsing KO Genes link\n";

## MAIN
my %final_hash;
my $counter = 0;
my $ref_fasta_db = Iterator::FastaDb->new($infile_genes) or die("Unable to open Fasta file, $infile_genes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    my $kegg_gene = "";
    if($header =~ m/^>(\S+).*/){
        $kegg_gene = $1;
        if(exists $hash_ref{$kegg_gene}){
            my $ko = $hash_ref{$kegg_gene};
            $final_hash{$ko}{$kegg_gene}{GENE} = $kegg_gene;
            $final_hash{$ko}{$kegg_gene}{HEADER} = $curr->header;
            $final_hash{$ko}{$kegg_gene}{SEQ} = $curr->seq;
        }
    }
    $counter++;
}


print STDERR Dumper(\%final_hash);

foreach my $ko (keys %final_hash) {
    open(OUT, ">>".$outdir."/".$ko.".faa");
    foreach my $kegg_gene (keys %{ $final_hash{$ko} }) {
        print OUT $final_hash{$ko}{$kegg_gene}{HEADER}."\n";
        print OUT $final_hash{$ko}{$kegg_gene}{SEQ}."\n";

    }
    close(OUT);
}
