#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
keepUblastBestHits.pl

PURPOSE:
Takes keeps best n hits from ublast table.

INPUT:
--infile <string>  : Sequence file
--e <float>        : max evalue
--length <int>     : min alignment length.
--perc <float>     : min perc id.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

Genomics and Microbiomes

ENDHERE

## OPTIONS
my ($help, $infile, $e, $length, $perc);
my $verbose = 0;

GetOptions(
  'infile=s'  => \$infile,
  'e=f'       => \$e,
  'length=i'  => \$length,
  'perc=f'    => \$perc,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

die("--infile missing...\n") unless($infile);
$e = "1e-05" unless($e);
## MAIN

my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
print STDERR "Processing $infile\n";
while(<IN>){
  chomp;
  if($_ =~ m/^#/){
    next;
  }

  my @row = split(/\t/, $_);
  my $curr_id = $row[0];
  my $curr_target = $row[1];
  my $curr_evalue = $row[10];
  my $curr_length = $row[3];
  my $curr_percid = $row[2];

  if($curr_evalue < $e && $curr_length >= $length && $curr_percid >= $perc){
    print STDOUT $_."\n";
  }
}
exit;
