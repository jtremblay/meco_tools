#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string> : domtblout or tblout concatenated files
--e_value <flaot> : hits having e-value higher than this will be rejected.
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $evalue);
my $verbose = 0;

GetOptions(
   'infile=s'   => \$infile,
   'e_value=f'  => \$evalue,
   'verbose'    => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $header = "target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value\tscore\tbias\tnumber\tof\tc-Evalue\ti-Evalue\tscore\tbias\tfrom\tto\tfrom\tto\tfrom\tto\tacc\tdescription of target";
print STDOUT $header."\n";
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   if($_ =~ m/^#/ && $. > 3){
      next;
   }elsif($_ =~ m/^#/){
      next;
   }else{
       my @row = split(/ {1,}/, $_);
       my @row2 = @row[0..21];
       my $length = scalar @row;
       my @desc = @row[22..($length-1)];
       my $desc = join(" ", @desc);
       my $row = join("\t", @row2);
       print STDOUT $row."\t";
       print STDOUT $desc."\n";
   }

   #exit if($. == 20);

}
close(IN);
