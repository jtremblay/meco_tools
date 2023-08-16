#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseNcbiTaxonomy.pl

PURPOSE:

INPUT:
--nodes <string> : nodes.dmp
--names <string> : names.dmp
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $nodes, $names);
my $verbose = 0;

GetOptions(
   'nodes=s'   => \$nodes,
   'names=s'   => \$names,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;
open(IN, "<".$nodes) or die "Can't open $nodes\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    
    my $tax_id = $row[0];
    my $parent_tax_id = $row[1];
    my $rank = $row[2];

}



#open(IN, "<".$names) or die "Can't open $names\n";
#while(<IN>){
#   chomp;
#   my @row = split(/\t/, $_);
#
#}
