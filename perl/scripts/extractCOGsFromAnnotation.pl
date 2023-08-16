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
--infile <string>   : Sequence file
--COG_list <string> : List of COGs (i.e. one COGxxxx per line).

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $KO_list);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'COG_list=s' => \$KO_list,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$KO_list) or die "Can't open $KO_list\n";
while(<IN>){
    chomp;
    $hash{$_} = "";
}
close(IN);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    if($. == 1){
       print STDOUT $_."\n"; 
    }else{
        my @row = split(/\t/, $_);
        if(exists $hash{$row[16]}){
            print STDOUT $_."\n";
        }
    }
}
close(IN);
