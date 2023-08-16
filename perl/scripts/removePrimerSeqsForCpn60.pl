#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
removePrimerSeqsForCpn60.p

PURPOSE:

INPUT:
--infile_primer <string> : txt file with one primer per line
--infile_fastq <string>  : fastq file
--length <int>           : length from beginning of read to check for primer presence. Default 40
--q <int>                : Quality cutoff. After cutting primer seq, read is rejected if average qual of read below that <int> value.				

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
$q = 20 unless($q);
$min_length = 1 unless($min_length);
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

    my $i = 0;
    my $ref_fastq_db = Iterator::FastqDb->new($infile_fastq) or die("Unable to open Fastq file, $infile_fastq\n");
    while( my $curr = $ref_fastq_db->next_seq() ) {
        my $seq = $curr->seq;
        if(length($seq) > $min_length){
        
            $seq = substr($seq, 0, $length);
            if(exists($hash{$seq})){
                #print STDERR $seq."\n";
                $hash{$seq}++;
                my @qual = unpack("C*", substr($curr->qual, $length));
                my $probsum = 0;
                foreach my $qual (@qual){
                    $qual = $qual - 33;
                    $probsum = $probsum + $qual;
                }
                my $average = $probsum / @qual;
                if($average >= $q){
                    # print STDERR $i."\t".$average."\n";
                    if($select_only){
                        print STDOUT $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
                    }else{
                        print STDOUT $curr->header."\n".substr($curr->seq, $length)."\n+\n".substr($curr->qual, $length)."\n";
                    }
                }
            }
        }else{
            print STDERR $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n" if($verbose);
        
        }
        $i++;
    }
}

exit;
