#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
filterForAlignment.pl

PURPOSE:
Filter the raw FEATURE or ASV fasta file to include only FEATURE/ASV IDs that are present in the provided --infile_feature_table <string> file.

INPUT:
--infile_feature_table <string> : FEATURE table in Qiime format (qiime.org)
--infile_fasta <string>         : Fasta file
        
OUTPUT:
STDOUT

NOTES:


BUGS/LIMITATIONS:

 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile_feature_table, $infile_fasta);
my $verbose = 0;

GetOptions(
  'infile_feature_table=s' => \$infile_feature_table,
  'infile_fasta=s'         => \$infile_fasta,
  'verbose'                => \$verbose,
  'help'                   => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile_feature_table <string> required\n" unless($infile_feature_table);
die "--infile_fasta <string> required\n" unless($infile_fasta);

## MAIN
my %hash_feature_table;
open(FEATURE_IN, "<".$infile_feature_table) or die "Can't open feature table ".$infile_feature_table."\n";
while(<FEATURE_IN>){
  chomp;
  next if($_ =~ m/#/);
  my @row = split(/\t/, $_);
  $hash_feature_table{$row[0]} = $_;  
}
close(FEATURE_IN);
#print STDERR Dumper(\%hash_feature_table);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  my $header = $curr->header();
  $header =~ s/>//;
  if(exists $hash_feature_table{$header}){
    print STDOUT ">".$header."\n".$curr->seq."\n";
  }
}
exit;
