#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Cwd 'abs_path';
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
translateBinFnaToFaa.pl

PURPOSE:

INPUT:
--bin_ids <string>        : Text file containing one bin id per line.
--link <string>           : File linking bin_id with contig_id and gene_id
--infile_faa_ref <string> : Fasta file containing all genes from metagenome (gene_id).
--outdir <string>         : Directory where fasta faa files are written.
				
OUTPUT:
<STDOUT> <string>         : Fasta file (faa) containing all genes of a given bin
                            supplied into --bin_ids <string> arg.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $bin_ids, $link, $infile_faa_ref, $outdir);
my $verbose = 0;

GetOptions(
   'bin_ids=s'        => \$bin_ids,
   'link=s'           => \$link,
   'infile_faa_ref=s' => \$infile_faa_ref,
   'outdir=s'         => \$outdir,
   'verbose'          => \$verbose,
   'help'             => \$help
);
if ($help) { print $usage; exit; }

## MAIN
$outdir = abs_path($outdir);
system("mkdir -p ".$outdir);

my %hash_bin_to_gene;
my %hash_gene_to_bin;

open(IN, "<".$bin_ids) or die "Can't open $bin_ids\n";
while(<IN>){
    chomp;
    my $bin_id = $_;
    $hash_bin_to_gene{$bin_id}{bin_id} = $bin_id;     
}
close(IN);

open(IN, "<".$link) or die "Can't open $link\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $bin_id = $row[0];
    my $contig_id = $row[1];
    my $gene_id = $row[2];

    if(exists $hash_bin_to_gene{$bin_id}){
        #$hash_bin_to_gene{$bin_id}{gene_ids} .= $gene_id.";";
        $hash_gene_to_bin{$gene_id} = $bin_id;
    }
}
close(IN);

#:w
#print STDERR Dumper(\%hash_gene_to_bin);
#exit;

#my $i = 0;
my $ref_fasta_db = Iterator::FastaDb->new($infile_faa_ref) or die("Unable to open Fasta file, $infile_faa_ref\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    #print STDERR $header."\n";
    my ($gene_id) = $header =~ m/^>(.*)$/g;
    #print STDERR $gene_id."\n";
    if(exists $hash_gene_to_bin{$gene_id}){
        open(OUT, ">>".$outdir."/".$hash_gene_to_bin{$gene_id}.".faa");
        print OUT "".$header.";".$hash_gene_to_bin{$gene_id}."\n".$curr->seq."\n";
        close(OUT);
    }
    #$i++;
    #$if($i == 10){exit;}
}


