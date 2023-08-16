#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
fastaToFastq.pl

PURPOSE:
The name says it all.

INPUT:
--fasta <fasta_infile>

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $fasta);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'fasta=s'   => \$fasta,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fasta fasta required\n") unless $fasta;

## MAIN

my $db = new Iterator::FastaDb($fasta);
while(my $seq=$db->next_seq) {
  my $header = $seq->header;
  #$header =~ s/#/-/g;
  $header =~ s/^>/@/;# if(substr($header,0,1) eq ">");
  my $qual = "G" x length($seq->seq);
  print STDOUT $header."\n".$seq->seq()."\n+\n".$qual."\n";
}

exit;
