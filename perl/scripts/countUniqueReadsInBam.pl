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
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
    'infile=s' 	=> \$infile,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
  chomp;
  my @row = split(/\t/, $_);
  my $id = $row[0];
  #if(exists $hash{$id}){
    
  #}
  $hash{$id} = $id;
  
}
close(IN);

my $i=0;
for my $key(keys %hash){
  $i++;
}
print STDOUT "There are ".$i." reads mapped.\n";
exit;
