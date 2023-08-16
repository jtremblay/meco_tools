#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fixSilvaDb.pl

PURPOSE:

INPUT:
--infile_fasta <string> : Fasta sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta);
my $verbose = 0;

GetOptions(
    'infile_fasta=s' => \$infile_fasta,
    'verbose' 	     => \$verbose,
    'help' 		     => \$help
);
if ($help) { print $usage; exit; }
die "File is empty of does not exists $infile_fasta\n" if((!-e $infile_fasta) and (!-s $infile_fasta));
die "--infile_fasta missing\n" unless($infile_fasta);

## MAIN
my %hash;

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header();
	$header =~ s/>//;
		
	next if($header =~ /bacteria|archaea/gi);

	my @header = split(/\s/, $header);
	my $id = shift(@header);
	my @lineage = split(/;/, $header[0]);
	
	my $newLineage = "";
	
	my $i=0;
	foreach(@lineage){
		$newLineage .= "k__".$_.";" if($i == 0);
		$newLineage .= "p__".$_.";" if($i == 1);
		$newLineage .= "c__".$_.";" if($i == 2);
		$newLineage .= "o__".$_.";" if($i == 3);
		$newLineage .= "f__".$_.";" if($i == 4);
		$newLineage .= "g__".$_.";" if($i == 5);
		$i++;
	}

	my $seq = $curr->seq();
	$seq =~ s/\s//g;
	$seq =~ s/U/T/gi;
		
	print STDOUT ">".$id." ".$newLineage."\n".$seq."\n";
}

