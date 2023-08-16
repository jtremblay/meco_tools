#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
extractUniprot.pl

PURPOSE:

INPUT:
--infile_uniprot <string> : tab file having all 6 char uniprot ID 
                            in the 4th column.
--infile_embl <string>    : Huge file from EMBL combining all 
                            annotations.

OUTPUT:
STDOUT                    : embl file having entries only for 
                          : for uniprot id in --infile_uniprot.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_uniprot, $infile_embl);
my $verbose = 0;

GetOptions(
    'infile_uniprot=s' 	=> \$infile_uniprot,
	'infile_embl=s' 	=> \$infile_embl,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## MAIN
die "--infile_uniprot arg missing...\n" unless($infile_uniprot);
die "--infile_embl arg missing...\n" 	unless($infile_embl);

open(UNIPROT, '<'.$infile_uniprot) or die "Can't open file $infile_uniprot\n";
open(EMBL, '<'.$infile_embl) or die "Can't open file $infile_embl\n";

my %hash;
while(<UNIPROT>){
	chomp;
	my @row = split(/\t/, $_);
	$hash{$row[3]}{'geneid'} = $row[0];
	$hash{$row[3]}{'transcript_stable_id'} = $row[1];
}
close(UNIPROT);

while(<EMBL>){
	chomp;
	my ($id) = split(/\t/, $_);
	if(exists $hash{$id}){
		print STDOUT $_."\n";
	}
}
close(EMBL);


