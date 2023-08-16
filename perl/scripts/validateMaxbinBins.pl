#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
validateGroopmBins.pl

PURPOSE:

INPUT:
--fasta <string>               : fasta file of contigs (contigs belonging to one bin)
--taxonomy <string>            : taxonomy.tsv file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $fasta, $taxonomy);
my $verbose = 0;

GetOptions(
   'fasta=s'    => \$fasta,
   'taxonomy=s' => \$taxonomy,
   'verbose' 	 => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--fasta is missing\n" unless($fasta);
die "--taxonomy is missing\n" unless($taxonomy);

## MAIN

# Parse fasta.
my %hash_bin;

my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my ($contig_id) = $curr->header =~ /^>(\S+)/;
   #$contig_id = $curr->header;
   $contig_id =~ s/^>//;

   $hash_bin{$contig_id} = "";

}

my %hash_tax;
open(IN, "<".$taxonomy) or die "Can't open $taxonomy\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $species = join(".", $row[8]);
   my $lineage = join("\t", $row[2], $row[3], $row[4], $row[5], $row[6], $row[7], $species);
   
   if(exists $hash_bin{$contig_id}){
       print STDOUT $contig_id."\t".$lineage."\t".$species."\n";
   }
}
close(IN);


exit;
