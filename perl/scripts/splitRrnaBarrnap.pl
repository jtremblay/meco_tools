#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
splitRnaBarrnap.pl

PURPOSE:

INPUT:
--infile_fasta <string> : Sequence file of contigs. rrna seqs will be extracted from this file using coords in --infile_tab
				
OUTPUT:
--fiveS <string> : ssu 5S
--SSU <string>   : ssu 16S or 18S
--LSU <string>   : lsu 23S or 28S

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $fiveS, $SSU, $LSU, $twelveS, $fiveeightS);
my $verbose = 0;

GetOptions(
   'infile_fasta=s' => \$infile_fasta,
   'fiveS=s'        => \$fiveS,
   'fiveeightS=s'   => \$fiveeightS,
   'twelveS=s'      => \$twelveS,
   'SSU=s'          => \$SSU,
   'LSU=s'          => \$LSU,
   'verbose' 	    => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#open(OUT_twelveS, ">".$twelveS)  or die "Can't open $twelveS\n";


#print STDERR Dumper(\%hash);

my $i5S = 0;
my $i58S = 0;
my $i12S = 0;
my $iSSU = 0;
my $iLSU = 0;

open(OUT_fiveS, ">".$fiveS)  or die "Can't open $fiveS\n";
close(OUT_fiveS);
open(OUT_fiveeightS, ">".$fiveeightS)  or die "Can't open $fiveeightS\n";
close(OUT_fiveeightS);
open(OUT_twelveS, ">".$twelveS)  or die "Can't open $twelveS\n";
close(OUT_twelveS);
open(OUT_SSU, ">".$SSU) or die "Can't open $SSU\n";
close(OUT_SSU);
open(OUT_LSU, ">".$LSU) or die "Can't open $LSU\n";
close(OUT_LSU);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my($header) = $curr->header =~ m/^>(\S+)/;

    if($header =~ m/5S/){
        if($i5S == 0){
            open(OUT_fiveS, ">".$fiveS)  or die "Can't open $fiveS\n";
        }else{
            open(OUT_fiveS, ">>".$fiveS)  or die "Can't open $fiveS\n";
        }
        print OUT_fiveS $curr->output;
        close(OUT_fiveS);
        $i5S++;

    }elsif($header =~ m/5.8S/){
        if($i58S == 0){
            open(OUT_fiveeightS, ">".$fiveeightS)  or die "Can't open $fiveeightS\n";
        }else{
            open(OUT_fiveeightS, ">>".$fiveeightS)  or die "Can't open $fiveeightS\n";
        }
        print OUT_fiveeightS $curr->output;
        close(OUT_fiveeightS);
        $i58S++;
    
    }elsif($header =~ m/12S/){
        if($i12S == 0){
            open(OUT_twelveS, ">".$twelveS)  or die "Can't open $twelveS\n";
        }else{
            open(OUT_twelveS, ">>".$twelveS)  or die "Can't open $twelveS\n";
        }
        print OUT_twelveS $curr->output;
        close(OUT_twelveS);
        $i12S++;

    }elsif($header =~ m/16S|18S/){
        if($iSSU == 0){
            open(OUT_SSU, ">".$SSU) or die "Can't open $SSU\n";
        }else{
            open(OUT_SSU, ">>".$SSU) or die "Can't open $SSU\n";
        }
        print OUT_SSU $curr->output;
        close(OUT_SSU);
        $iSSU++;

    }elsif($header =~ m/23S|28S/){
        if($iLSU == 0){
            open(OUT_LSU, ">".$LSU) or die "Can't open $LSU\n";
        }else{
            open(OUT_LSU, ">>".$LSU) or die "Can't open $LSU\n";
        }
        print OUT_LSU $curr->output;
        close(OUT_LSU);
        $iLSU++;
    }
}
close(OUT_fiveS); 
close(OUT_SSU);
close(OUT_LSU);

exit;
