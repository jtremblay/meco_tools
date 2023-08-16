#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseSilvaDb.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#open(IN, "<".$infile) or die "Can't open $infile\n";
#while(<IN>){
#   chomp;
#   $hash{$_} = $_;
#}
#close(IN);

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    if($header =~ m/g__Bacillus/){
        print STDOUT $curr->header."\n".$curr->seq."\n";
    }
    if($header =~ m/g__Lactobacillus/){
        print STDOUT $curr->header."\n".$curr->seq."\n";
    }
    if($header =~ m/Paenibacillus/){
        print STDOUT $curr->header."\n".$curr->seq."\n";
    }
}

exit;
