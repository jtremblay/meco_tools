#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
removePrimerSeqsForBenFlo.pl

PURPOSE:

INPUT:
--infile_primer <string> : txt file with one primer per line
--infile_fastq <string>  : fastq file

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_primer, $infile_fastq, $length, $q, $select_only, $min_length);
my $verbose = 0;

GetOptions(
   'infile_primer=s' => \$infile_primer,
   'infile_fastq=s'  => \$infile_fastq,
   'length=i'        => \$length,
   'q=i'             => \$q,
   'select_only'     => \$select_only,
   'min_length=i'    => \$min_length,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $i = 0;
my $matches = 0;
if($infile_primer){
    my $primer_seq = "";
    open(IN, "<".$infile_primer) or die "Can't open $infile_primer\n";
    while(<IN>){
        chomp;
        $primer_seq = $_;
    }
    close(IN);

    my $ref_fastq_db = Iterator::FastqDb->new($infile_fastq) or die("Unable to open Fastq file, $infile_fastq\n");
    while( my $curr = $ref_fastq_db->next_seq() ) {
        my $seq = $curr->seq;
        if($seq =~ m/($primer_seq)/){
            my $start = $-[0];
            my $end = $+[0];
            $matches++;

            print STDOUT $curr->header."\n".substr($curr->seq, $start)."\n+\n".substr($curr->qual, $start)."\n";
        }else{
            print STDOUT $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
        }
        $i++;
    }
}
print STDERR "Number of sequences: $i\n";
print STDERR "Number of sequences trimmed: $matches\n";

exit;
