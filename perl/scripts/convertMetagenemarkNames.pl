#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
convertMetagenemarkNames.pl

PURPOSE:

INPUT:
--gff <string>         : gff infile
--fna <string>         : fasta nucl. file

OUTPUT:
STDOUT                 : renamed fna seqs.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $gff, $fna);
my $verbose = 0;

GetOptions(
   'gff=s' 	   => \$gff,
   'fna=s' 	   => \$fna,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
my $counter = 0;
open(GFF, "<".$gff) or die "Can't open $gff\n";
while(<GFF>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   # Replace all " characters with ' chars. HtSeq has a problem with " chars...
   $_ =~ s/\"/\'/g;
   my @row = split(/\t/, $_);
   my $originalId = $row[0];
   my @field = split(/,/, $row[8]);
   my $mgmId = $field[0];
   $mgmId =~ s/=/_/;
   $originalId =~ s/\s/\./g;
   $hash{$mgmId} = $originalId."=".$mgmId;
   $counter++;
}
close(GFF);

my $ref_fasta_db = Iterator::FastaDb->new($fna) or die("Unable to open Fasta file, $fna\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   if(exists $hash{$header}){
      print STDOUT ">".$hash{$header}."\n".$curr->seq."\n";
   }
}

