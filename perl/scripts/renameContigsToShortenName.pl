#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Basename;
#use Cwd;
#use File::Spec;
use Path::Class qw(dir);

my $usage=<<'ENDHERE';
NAME:
renameContigsToShortenName.pl

PURPOSE:

INPUT:
--infile <string> : infile
				
OUTPUT:
STDOUT

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
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#my %hash;
#open(IN, "<".$gene_list) or die "Can't open $gene_list\n";
#while(<IN>){
#    chomp;
#    $hash{$_} = $_;
#}
#close(IN);

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my ($header) = $curr->header =~ m/^(\S+)/;
    my $seq = $curr->seq;
    print STDOUT $header."\n".$curr->seq."\n";
    
}  

