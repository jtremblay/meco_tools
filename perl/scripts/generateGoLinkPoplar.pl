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
--infile_uniprot <string>   : Uniprot ref file. Having xref ids.
--infile_idmapping <string> : infile_id_mapping file.

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_uniprot, $infile_idmapping);
my $verbose = 0;

GetOptions(
    'infile_uniprot=s' 	 => \$infile_uniprot,
	'infile_idmapping=s' => \$infile_idmapping,
    'verbose' 	         => \$verbose,
    'help' 		         => \$help
);
if ($help) { print $usage; exit; }

die "--infile_uniprot required\n" unless($infile_uniprot);
die "--infile_idmapping required\n" unless($infile_idmapping);

## MAIN
open(UNIPROT, '<'.$infile_uniprot) or die "Can't open file ".$infile_uniprot."\n";
open(IDMAPPING, '<'.$infile_idmapping) or die "Can't open file ".$infile_idmapping."\n";

my %hash;
while(<UNIPROT>){
	chomp;
	my @row = split(/\t/, $_);
	
	my $geneId = $row[0];
	my $xref = $row[3];
	
	$hash{$xref} = $geneId;
}
close(UNIPROT);

while(<IDMAPPING>){
	chomp;
	my @row = split(/\t/, $_);
	my $xref = $row[0];

	if(exists $hash{$xref}){
		foreach my $el (@row){
			if($el =~ m/(GO\:\d+)/){
				print STDOUT $hash{$xref}."\t".$1."\n";
			}
		}
	}
}

close(IDMAPPING);
