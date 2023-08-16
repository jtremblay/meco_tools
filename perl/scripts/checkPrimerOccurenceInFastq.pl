#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
checkPrimerOccurenceInFastq.pl

PURPOSE:

INPUT:
--infile_primer <string> : txt file with one primer per line. Optional.
--infile_fastq <string>  : fastq file
--length <int>           : length from beginning of read to check for primer presence. Default 40
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_primer, $infile_fastq, $length);
my $verbose = 0;

GetOptions(
   'infile_primer=s' => \$infile_primer,
   'infile_fastq=s'  => \$infile_fastq,
   'length=i'        => \$length,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
if($infile_primer){
    open(IN, "<".$infile_primer) or die "Can't open $infile_primer\n";
    while(<IN>){
       chomp;
       $hash{$_} = 0;
       #my @row = split(/\t/, $_);
       $length = length($_);
    }
    close(IN);

    my $ref_fastq_db = Iterator::FastqDb->new($infile_fastq) or die("Unable to open Fastq file, $infile_fastq\n");
    while( my $curr = $ref_fastq_db->next_seq() ) {
        my $seq = $curr->seq;
        $seq = substr($seq, 0, $length);
        if(exists($hash{$seq})){
            $hash{$seq}++;
        }
    }
}else{
    my $ref_fastq_db = Iterator::FastqDb->new($infile_fastq) or die("Unable to open Fastq file, $infile_fastq\n");
    while( my $curr = $ref_fastq_db->next_seq() ) {
        my $seq = $curr->seq;
        $seq = substr($seq, 0, $length);
        $hash{$seq}++;
    }
}
foreach my $name (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
    print STDOUT $name."\t".$hash{$name}."\n";
}

print STDERR Dumper(\%hash);

exit;
