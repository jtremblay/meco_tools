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
STDOUT

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
open(IN, '<'.$infile) or die "Can't open file $infile\n";

my %hash;
while(<IN>){
	chomp;
	my @row = split(/\t/, $_);
	my $geneId = $row[1];

	my %seen;	
	foreach my $el (@row){
		if($el =~ m/(GO:\d+)/){
			my $GO = $1;
			if(exists $seen{$GO}){
	
			}else{
				#print STDOUT $geneId."\t".$GO."\n";
				$hash{$geneId} = $GO;
			}
			$seen{$GO} = 1;
		}
	}	
}
close(IN);

for my $key (sort {$a cmp $b} keys %hash){
	print STDOUT $key."\t".$hash{$key}."\n";
}

exit;

