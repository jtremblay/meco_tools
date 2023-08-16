#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
clwToSto.pl

PURPOSE:

INPUT:
--infile <string> : Alignment file in clustal format.
				
OUTPUT:
STDOUT            : Alignment in Stockholm format.

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

open(IN, "<".$infile) or die "Can't open $infile\n";
print STDOUT "# STOCKHOLM 1.0\n";
while(<IN>){
   chomp;
   if($. == 1){ next; }

   if($_ =~ m/^\S+/){
      my $line = $_;
      $line =~ s/-/\./g;
      print STDOUT $line."\n";
   }

   if($_ =~ m/^\s/){
      print STDOUT "\n";
   }

}

print STDOUT "//";
