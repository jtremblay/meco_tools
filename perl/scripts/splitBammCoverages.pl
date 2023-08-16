#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;
use Text::CSV;

my $usage=<<'ENDHERE';
NAME:
splitBammCoverages.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
--outdir <string> : Output directory where to write coverage of individual files.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outdir);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'outdir=s'  => \$outdir,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#my %hash;
#open(IN, "<".$infile) or die "Can't open $infile");
#while(<IN>){
#   chomp;
#   if($_ =~ m/^#/){
#      my @row = split(/\t/, $_)
#      next;
#   }
#}

my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: ".Text::CSV->error_diag ();
open my $fh, "<:encoding(utf8)", $infile or die "Can't open $infile\n";
while ( my $row = $csv->getline( $fh ) ) {
   $row->[2] =~ m/pattern/ or next; # 3rd field should match
   push @rows, $row;
}
$csv->eof or $csv->error_diag();
close $fh;

$csv->eol ("\r\n");


open $fh, ">:encoding(utf8)", "./new.tsv" or die "./new.tsv: $!";
$csv->print ($fh, $_) for @rows;
close $fh or die "new.csv: $!";



