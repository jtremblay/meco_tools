#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseUsearch.pl

PURPOSE:

INPUT:
--infile_uc <string>    : usearch file format
--infile_fasta <string> : Sequence file
				
OUTPUT:
STDOUT                  : Sequence fasta file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_uc, $infile_fasta);
my $verbose = 0;

GetOptions(
   'infile_uc=s' 	   => \$infile_uc,
   'infile_fasta=s' 	=> \$infile_fasta,
   'verbose' 	      => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;

# Parse usearch
open(IN, "<".$infile_uc) or die "Can't open $infile_uc\n";
while(<IN>){
   chomp;
   if($_ =~ m/^S/){
      my @row = split(/\t/, $_);
      my $header= $row[8];
      my @headers = split(/\s+/, $header);
      my $id = $headers[0];
      $hash{$id} = $header;
   }
}
close(IN);

print STDERR Dumper(\%hash);

# Parse fasta file to extract only seed sequences.
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   my @headers = split(/\s+/, $header);
   my $id = $headers[0];
   $id =~ s/>//;

   if(exists $hash{$id}){
      print STDOUT $curr->header."\n".$curr->seq."\n";
   }
}
exit;

