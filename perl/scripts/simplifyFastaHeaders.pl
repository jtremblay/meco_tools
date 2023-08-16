#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
simplifyFastaHeaders.pl

PURPOSE:
Simplify fasta headers. For anvio.

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
<STDOUT>          : Sequence file fasta.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash = (
'R' => 'A', # or G
'Y'	=> 'C', # or T
'S'	=> 'G', # or C
'W'	=> 'A', # or T
'K' => 'G', # or T
'M'	=> 'A', # or C
'B' => 'C', # or G or T
'D'	=> 'A', # or G or T
'H'	=> 'A', # or C or T
'V' => 'A' # or C or G
);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   if($_ =~ m/^>/){
       chomp;
       my ($new_header) = $_ =~ m/^(\S+)/;
       #$new_header =~ s/\./_/g;
       print STDOUT $new_header."\n";
   
   }else{
       my $seq = $_;
       $seq =~ tr/ACGTRYSWKMBDHV\.-/ACGTACGAGACAAA\.-/;
       print STDOUT uc($seq);
   }
}
exit;
