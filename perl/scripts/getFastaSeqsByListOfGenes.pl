#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
getFastaSeqsByListOfGenes.pl

PURPOSE:

INPUT:
--infile    <string>   : Sequence file
--gene_list <string>   : Term to match.

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $gene_list);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'gene_list=s' => \$gene_list,
   'verbose' 	 => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$gene_list) or die "Can't open $gene_list\n";
while(<IN>){
   chomp;
   $hash{$_} = $_;
}
print STDERR Dumper(\%hash);
#exit;
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   my $parsed_header;
   if($header =~ m/(^\S+) /){
      $parsed_header = $1;
      #print STDERR $parsed_header."\n";
   }elsif($header =~ m/(^\S+)/){
      $parsed_header = $1;
      #print STDERR $parsed_header."\n";
   }
   
   if(exists $hash{$parsed_header}){
      print STDOUT $curr->output;
   }
}
exit;
