#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateGoLinkPoplar.pl

PURPOSE:

INPUT:
--infile_GO <string>     : GO.
--infile_GOlink <string> : GOlink

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_golink, $infile_go);
my $verbose = 0;

GetOptions(
    'infile_golink=s' 	 => \$infile_golink,
	'infile_go=s' => \$infile_go,
    'verbose' 	         => \$verbose,
    'help' 		         => \$help
);
if ($help) { print $usage; exit; }

die "--infile_golink required\n" unless($infile_golink);
die "--infile_go required\n" unless($infile_go);

## MAIN
open(GOLINK, '<'.$infile_golink) or die "Can't open file ".$infile_golink."\n";
open(GO, '<'.$infile_go) or die "Can't open file ".$infile_go."\n";

my %hash;
while(<GO>){
	chomp;
	my @row = split(/\t/, $_);
	
	my $geneId = $row[0];
	my $symbol = $row[1];
	
	$hash{$geneId} = $symbol;
}
close(GO);

while(<GOLINK>){
	chomp;
	my @row = split(/\t/, $_);
	my $geneId = $row[0];
	my $go = $row[1];

	if(exists $hash{$geneId}){
		print STDOUT $hash{$geneId}."\t".$go."\n";
	}
}

close(GOLINK);
