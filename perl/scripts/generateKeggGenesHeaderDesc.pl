#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
generateKeggGenesHeaderDesc.pl

PURPOSE:
From the genes.pep file (KEGG), generate a tab file linking
the gene_id to product name specified in the sequence headers.

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT <string>   : Tab delimited outfile.

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
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;

    #my $gene_id = "";
    #my $desc = "";
    #print STDERR "header:\n$header\n";
    if($header =~ m/^(>\S+)\s+(.*)$/){
         my $gene_id = $1;
         my $desc = $2;
         #print STDERR "1: $1\n";
         #print STDERR "2: $2\n";
         $gene_id =~ s/^>//;
         print STDOUT $gene_id."\t".$desc."\n";
    }else{
        print STDERR $header."\n";
        die "Could not parse headers...\n";
    }
    

}
exit;


