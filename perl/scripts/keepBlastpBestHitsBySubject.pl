#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
keepBlastpBestHitsBySubject.pl

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
   my $curr_length = $row[4];
   my $curr_percid = $row[2];
 
   #print STDERR $curr_evalue."\t".$e."\n";
   #print STDERR $curr_length."\t".$length."\n";
    #print STDERR $curr_percid."\t".$perc."\n";
        #  ksgA        gene_id_12345   1e-10
    $hash{$curr_target}{$curr_id} = $curr_evalue;
 
    #if($curr_evalue < $e && $curr_length >= $length && $curr_percid >= $perc){
    # print STDOUT $_."\n";
    #}
}

foreach my $curr_target (keys %hash) {
   print STDOUT "#$curr_target\n";
   foreach my $curr_id (keys %{ $hash{$curr_target} }) {
      if($hash{$curr_target}{$curr_id} <= $e){
         print STDOUT "$curr_id\t$curr_target\t$hash{$curr_target}{$curr_id}\n";
      }
   }
}

exit;
