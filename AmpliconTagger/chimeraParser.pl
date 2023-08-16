#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
chimeraParser.pl

PURPOSE:

INPUT:
--infile <string>       : uchime output file (tab delimited).
--infile_fasta <string> : infile fasta file.
--cutoff <float>        : alignment perc cutoff.
				
OUTPUT:
STDOUT                  : outfile fasta. Fasta entries having %align higher than 
                          cutoff.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $infile_fasta, $cutoff);
my $verbose = 0;

GetOptions(
   'infile=s'       => \$infile,
   'infile_fasta=s' => \$infile_fasta,
   'cutoff=f'       => \$cutoff,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $perc_id = $row[8];
   #my ($otu_id) = $row[1] =~ m/(\d+);/;
   my $otu_id = $row[1];
   if($perc_id eq "*"){
     $hash{$otu_id} = $perc_id;
   }elsif($perc_id >= $cutoff){
     $hash{$otu_id} = $perc_id;
   }
}
close(IN);
print STDERR Dumper(\%hash);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my($header) = $curr->header =~ m/^>(.*)$/;
    $header =~ $1;
    #print STDERR $header."\n";;
    if(exists $hash{$header}){
        print STDOUT ">".$header."\n".$curr->seq."\n";
    }else{
        print STDERR ">".$header."\n".$curr->seq."\n";
    }
}

