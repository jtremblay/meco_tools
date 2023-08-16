#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
bookChapter1.pl

PURPOSE:

INPUT:
--infile <string>    : list of strings 
--reference <string> : KEGG ref pathways file
				
OUTPUT:
STDOUT               : list of KO and ko separated by a comma

NOTES:
Written in the context of my oil degradation book chapter.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $ref);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'reference=s' => \$ref,
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $string = $row[0];

    print $string."====\n";
    open(IN2, "<".$ref) or die "Can't open $ref\n";
    my %hash;
    while(<IN2>){
        chomp;
        my @row = split(/\t/, $_);
        my $ko =  $row[4];
        my $KO =  $row[6];
        if($_ =~ m/$string/i){
            $hash{$ko}{ko} = $ko;
            $hash{$KO}{KO} = $KO;
        }
    }
    for my $key (%hash){
        if(exists $hash{$key}{KO}){
            print $hash{$key}{KO},",";
        }
    }
    print "\n";
    for my $key (%hash){
        if(exists $hash{$key}{ko}){
            print $hash{$key}{ko},",";
        }
    }
    print "\n";
    print "========\n";
    close(IN2);
}
close(IN);

