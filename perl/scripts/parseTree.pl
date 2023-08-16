#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseTree.pl

PURPOSE:

INPUT:
--infile <string>        : fasta file (multialignemnt file)
--link <string>          : tsv file bin-contig-gene				
--selected_bins <string>  : selected_bins

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $link, $selected);
my $verbose = 0;

GetOptions(
   'infile=s' 	      => \$infile,
   'link=s'          => \$link,
   'selected_bins=s' => \$selected,
   'verbose' 	      => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash_selected;
open(IN, "<".$selected) or die "Can't open $selected\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   $hash_selected{$row[0]} = $row[0];
}
close(IN);

my %hash;
open(IN, "<".$link) or die "Can't open $link\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   if(exists $hash_selected{$row[0]}){
      $hash{$row[2]} = $row[0];
   }
}
close(IN);

my %hash2;
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   if(exists $hash{$header}){
       #print STDOUT ">".$header."-".$hash{$header}."\n".$curr->seq."\n";
       $hash2{$hash{$header}."-".$header} = ">".$header."-".$hash{$header}."\n".$curr->seq."\n"
   }
}

for my $key (sort keys %hash2){
   print STDOUT $hash2{$key};
}




