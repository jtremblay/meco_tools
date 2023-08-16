#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
concatenateContigsWithNewIDs.pl

PURPOSE:

INPUT:
--infiles <string> : Sequence files separated by a <,>
--names <string>   : Names separated by a <,>				

OUTPUT:

NOTES:
Was written for merging Ray contigs.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $names);
my $verbose = 0;

GetOptions(
   'infiles=s' => \$infiles,
   'names=s' 	=> \$names,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @infiles = split(/,/, $infiles);
my @names = split(/,/, $names);


foreach my $infile (@infiles){

   my $name = shift(@names);

   my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my($header) = $curr->header() =~ m/^(>\S+)/;
      $header = $header."-".$name;
      print STDOUT $header."\n";

      my $splitseq = join ("\n",  ( $curr->seq =~ /.{1,80}/gs ));
      print STDOUT "$splitseq\n";
   }
}
exit;
