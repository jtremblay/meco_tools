#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fixPacBioQscores.pl

PURPOSE:
Fix CCS qscores that are between 33 and 81 with PacBio CCS.
Replace Qscores higher than 40 with a value of 40.

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $size =  -s $infile or die "$! : $infile".". Cannot calculate size of file.\n";
#my $quarter = int( $size / $num_threads );
#my @starts = map $quarter * $_, 0 .. ($num_threads - 1);

open my $FASTQ, '<', $infile or die $!;
while( tell( $FASTQ ) < $size ) {
    my @lines = map scalar( <$FASTQ> ), 1 .. 4;
    chomp @lines;
    my $header = $lines[0];
    my $seq = $lines[1];
    my $qual = $lines[3];

    my @qual = split(//, $qual);
    my @newqual;
    foreach my $q (@qual){
     
        # convert qual higher than 40 to 40.
        my $ascii = ord($q);
        if($ascii > 40){
            $ascii = 33+40;
        }
        push @newqual, chr($ascii);
    }
    print STDOUT $header."\n".$seq."\n+\n".join("", @newqual)."\n";
}
