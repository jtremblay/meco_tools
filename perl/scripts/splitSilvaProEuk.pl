#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
splitSilvaProEuk.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
--outfile_pro <string> : sequence file
--outfile_euk <string> : sequence file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outfile_pro, $outfile_euk);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'outfile_pro=s' => \$outfile_pro,
   'outfile_euk=s' => \$outfile_euk,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT_PRO, ">".$outfile_pro) or die "Can't open $outfile_pro\n";
open(OUT_EUK, ">".$outfile_euk) or die "Can't open $outfile_euk\n";


my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    if($curr->header =~ m/bacteria/i || $curr->header =~ m/archaea/i){
        print OUT_PRO $curr->output;
    }else{
        print OUT_EUK $curr->output;
        
    }
}
close(OUT_PRO);
close(OUT_EUK);
exit;

