#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseProbesForChad.pl

PURPOSE:

INPUT:
--infile_blast <string> : Sequence file
--infile_fasta <string> : Sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_blast, $infile_fasta);
my $verbose = 0;

GetOptions(
   'infile_blast=s'	=> \$infile_blast,
   'infile_fasta=s'	=> \$infile_fasta,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//g;
   $hash{$header}{SEQ} = $curr->seq;
   $hash{$header}{SUBJECT_LENGTH} = length($curr->seq);
   $hash{$header}{MAX} = 0;
}

open(IN, "<".$infile_blast) or die "Can't open $infile_blast\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $ali_length = $row[3];
   my $query = $row[0];
   my $subject = $row[1];

   if(exists $hash{$subject}){
       #print STDERR "exists!!\n";
      if($ali_length > $hash{$subject}{MAX} ){
         $hash{$subject}{QUERY} = $query;
         $hash{$subject}{MAX} = $ali_length;
      }
   }
}
close(IN);

print STDOUT "subject\tquery\tsubject_length\tquery_length\n";
for my $key (keys %hash){
   print STDOUT "$key\t$hash{$key}{QUERY}\t$hash{$key}{SUBJECT_LENGTH}\t$hash{$key}{MAX}\n";
}

print STDERR Dumper(\%hash);

