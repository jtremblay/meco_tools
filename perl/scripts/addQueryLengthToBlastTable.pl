#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile_fasta <string>    : Sequence file
--infile_blastout <string> : Sequence file
				
OUTPUT:
STDOUT

NOTES:
--infile_blastout can be any file as long as the first column represents the contig id or
same fasta header entry as the ones in the fasta file.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_blastout);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'    => \$infile_fasta,
   'infile_blastout=s' => \$infile_blastout,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/>//g;
    $hash{$header} = length($curr->seq);
}
 
open(IN, "<".$infile_blastout) or die "Can't open $infile_blastout\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    my $contig_id = $row[0];
    if(exists $hash{$contig_id}){
        print STDOUT $_."\t".$hash{$contig_id}."\n";
    }else{
        print STDERR "Something went wrong, did not find $contig_id corresponding fasta entry...\n";
    }
}
close(IN);
