#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(uniq);

my $usage=<<'ENDHERE';
NAME:
getKeggModulesKOLink.pl

PURPOSE:

INPUT:
--infiles <string> : Modules html files seperated by a ,
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles);
my $verbose = 0;

GetOptions(
   'infiles=s' => \$infiles,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @infiles = split(/,/, $infiles);

foreach my $infile (@infiles){
    my ($basename) = $infile =~ m/(M\d{5})/;
    open(IN, "<".$infile) or die "Can't open $infile\n";
    
    my @matches;
    while(<IN>){
       chomp;
       my @tmp_matches = $_ =~ m/(K\d{5})/g;
       push(@matches, @tmp_matches);
    }
    @matches = uniq @matches;
    foreach my $match (@matches){
        print STDOUT $basename."\t".$match."\n";
    }
    close(IN);
}

exit;
