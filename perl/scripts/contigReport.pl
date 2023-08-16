#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
contigReport.pl

PURPOSE:
Intended for assemblies. Will display number of 
sequence, number of bases in each sequences. blabla. 

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
my $totalBases = 0;
my $counter = 1;

my %hash;

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $length = length($curr->seq);
	my $header = $curr->header;
	$header =~ s/>//;

	#print STDOUT "Sequence: ".$header."\t".$length." bp\n";
	$totalBases += $length;
	$hash{$header} = $length;
	$counter++;
}
$counter = ($counter - 1);
print STDOUT "Total of ".$counter." sequences\n";


foreach my $key (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
	print STDOUT $key . "\t" .$hash{$key}."\n";	
}

my ($N25, $N50, $N75);
my $cummSum = 0;
$counter = 1;
foreach my $key (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
	$cummSum += $hash{$key};
	my $ratio = int( ($cummSum / $totalBases) * 100);

	if($ratio >= 25 && !$N25){
		print STDOUT "N25 - 25% of total sequence length is contained in the ".$counter." sequence(s) having a length >= ".$hash{$key}." bp\n";
		$N25=1;
	}
	if($ratio >= 50 && !$N50){
		print STDOUT "N50 - 50% of total sequence length is contained in the ".$counter." sequence(s) having a length >= ".$hash{$key}." bp\n";
		$N50=1;
	}
	if($ratio >= 75 && !$N75){
		print STDOUT "N75 - 75% of total sequence length is contained in the ".$counter." sequence(s) having a length >= ".$hash{$key}." bp\n";
		$N75=1;
	}
	$counter++;
}

exit;

