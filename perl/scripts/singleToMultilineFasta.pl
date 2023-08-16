#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
singleToMultilineFasta.pl

PURPOSE:

INPUT:
--infile <string>         : Sequence file
--num_char_per_line <int> : default = 80 

OUTPUT:
STDOUT <string>   : Sequence file with corrected def lines.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $num_char);
my $verbose = 0;

GetOptions(
   'infile=s'            => \$infile,
   'num_char_per_line=i' => \$num_char,
   'verbose'             => \$verbose,
   'help'                => \$help
);
if ($help) { print $usage; exit; }

$num_char = 80 unless($num_char);

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   my ($corrected_header) = $header =~ m/^(>\S+)/;

   print STDOUT $corrected_header."\n";
   my $seq = $curr->seq;
   my $old = $seq;
   my $new="";
   while (length($old)> $num_char) {
       $new .= substr($old,0, $num_char)."\n";
       $old = substr($old, $num_char);
   }
   $new .= $old."\n";
   print STDOUT $new;
}

