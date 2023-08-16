#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
binmode STDOUT, ":utf8";

my $usage=<<'ENDHERE';
NAME:
454ToFastq.pl

PURPOSE:

INPUT:
--fasta <string>    : Sequence files separated by a ,
--qual <string>     : Qual files separated by a ,

OUTPUT:
--outIndex <string> : Indexes (barcodes)				
STDOUT (Fastq reads).

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $fasta, $qual, $outIndex);
my $verbose = 0;

GetOptions(
    'fasta=s' 	 => \$fasta,
	'qual=s'     => \$qual,
	'outIndex=s' => \$outIndex,
    'verbose' 	 => \$verbose,
    'help' 		 => \$help
);
if ($help) { print $usage; exit; }

## MAIN
die("--outIndex missing\n") unless($outIndex);
open(OUTINDEX, ">".$outIndex) or die "Can't open $outIndex\n";
#die("--fasta and --qual must be of equal length\n") if(@fasta != @qual);

# Get barcodes
my $i=0;
my $j=0;
my @barcodes;
for(glob '{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}'){
	if($i % 30 == 0){ 
		$barcodes[$j] = $_; 
		$j++;
	}else{

	}   
	$i++;
} 

my @fasta = split(/,/, $fasta);
my @qual = split(/,/, $qual);

my $counter=0;
foreach my $fasta (@fasta){
	my %hash;
	my $qual = shift(@qual);

	print STDERR "[DEBUG] processing $fasta\n";
	print STDERR "[DEBUG] processing $qual\n";

	#open(FASTA, "<".$fasta) or die "Can't open $fasta\n";
	my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $header = $curr->header;
		$header =~ s/>//;
		$hash{$header}{fasta} = $curr->seq;
	}

	open(QUAL, "<".$qual) or die "Can't open $qual\n";
	my $tmp_title;
	while(<QUAL>){
		chomp;
    	if ($_ =~ /^>/) {
			$tmp_title = $_;
			$tmp_title =~ s/>//;
			#print $tmp_title."\n";
		} else {
			$hash{$tmp_title}{qual} .= $_." ";
		}          
	}
	close(QUAL);

	# THEN print as fastq
	for my $key (sort{$a cmp $b} keys %hash){
       #print STDERR "KEY: ".$key."\n";
		my $asciiQual = "";
		for( split /\s+/, $hash{$key}{qual} ){
          #print STDERR $_." ";
			$asciiQual .= chr($_ + 33);
		}
      #print STDERR "\n";
		my $header = $key;
		$header =~ s/ /_/g;
		$header .= "#".$barcodes[$counter]."/1";

		print STDOUT "@".$header."\n".$hash{$key}{fasta}."\n+\n".$asciiQual."\n";
	}

	print OUTINDEX $fasta."\t".$barcodes[$counter]."\n";

	$counter++;
}

close(OUTINDEX);
exit;     
 
