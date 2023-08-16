#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseFastaDefLines.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT <string>   : Sequence file with corrected def lines.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   my ($corrected_header) = $header =~ m/^(>\S+)/;

   print STDOUT $corrected_header."\n";
   my $seq = $curr->seq;
   my $old = $seq;
   my $new="";
   while (length($old)> 80) {
       $new .= substr($old,0, 80)."\n";
       $old = substr($old, 80);
   }
   $new .= $old."\n";
   print STDOUT $new;
}

