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
--infile_primer1 <string> : txt file with one primer per line
--infile_primer2 <string> : txt file with one primer per line
--infile_fastq1 <string>  : fastq file
--infile_fastq2 <string>  : fastq file
--length <int>            : length from beginning of read to check for primer presence. Default 40
--q <int>                 : Quality cutoff. After cutting primer seq, read is rejected if average qual of read below that <int> value.				

OUTPUT:
--outfile1 <string>
--outfile2 <string>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_primer1, $infile_primer2, $infile_fastq1, $infile_fastq2, 
    $length, $q, $select_only, $min_length, $outfile1, $outfile2);
my $verbose = 0;

GetOptions(
   'infile_primer1=s' => \$infile_primer1,
   'infile_primer2=s' => \$infile_primer2,
   'infile_fastq1=s'  => \$infile_fastq1,
   'infile_fastq2=s'  => \$infile_fastq2,
   'length=i'         => \$length,
   'q=i'              => \$q,
   'select_only'      => \$select_only,
   'min_length=i'     => \$min_length,
   'outfile1=s'       => \$outfile1,
   'outfile2=s'       => \$outfile2,
   'verbose'          => \$verbose,
   'help'             => \$help
);
if ($help) { print $usage; exit; }

## MAIN
$q = 20 unless($q);
$min_length = 1 unless($min_length);

open(OUT1, ">".$outfile1) or die "Can't open $outfile1\n";
open(OUT2, ">".$outfile2) or die "Can't open $outfile2\n";

my $primer1;
my $primer2;
open(IN, "<".$infile_primer1) or die "Can't open $infile_primer1\n";
while(<IN>){
   chomp;
   $primer1 = $_;
}
close(IN);
open(IN, "<".$infile_primer2) or die "Can't open $infile_primer2\n";
while(<IN>){
   chomp;
   $primer2 = $_;
}
close(IN);

my $i = 0;
my $ref_fastq_db1 = Iterator::FastqDb->new($infile_fastq1) or die("Unable to open Fastq file, $infile_fastq1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($infile_fastq2) or die("Unable to open Fastq file, $infile_fastq2\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
    my $curr2 = $ref_fastq_db2->next_seq();
    
    my $seq1 = $curr1->seq;
    my $seq2 = $curr2->seq;
    my $qual1 = $curr1->qual;
    my $qual2 = $curr2->qual;

    my $end1;
    my $end2;
    if(length($seq1) > $min_length && length($seq2) > $min_length){
    
        if($seq1 =~ m/($primer1)/){
            my $end1 = $+[0];

            if($seq2 =~ m/($primer2)/){
                my $end2 = $+[0];

                my $trimmed_seq1 = substr($seq1, $end1);
                my $trimmed_seq2 = substr($seq2, $end2);
                my $trimmed_qual1 = substr($qual1, $end1);
                my $trimmed_qual2 = substr($qual2, $end2);

                my @qual1 = unpack("C*", $trimmed_seq1);
                my @qual2 = unpack("C*", $trimmed_seq2);
                my $probsum1 = 0;
                my $probsum2 = 0;
                foreach my $qual (@qual1){
                    $qual = $qual - 33;
                    $probsum1 = $probsum1 + $qual;
                }
                foreach my $qual (@qual2){
                    $qual = $qual - 33;
                    $probsum2 = $probsum2 + $qual;
                }
                my $average1 = $probsum1 / @qual1;
                my $average2 = $probsum2 / @qual2;
                if($average1 >= $q && $average2 >= $q){
                    # print STDERR $i."\t".$average."\n";
                    if($select_only){
                        print OUT1 $curr1->header."\n".$curr1->seq."\n+\n".$curr1->qual."\n";
                        print OUT2 $curr2->header."\n".$curr2->seq."\n+\n".$curr2->qual."\n";
                    }else{
                        print OUT1 $curr1->header."\n".$trimmed_seq1."\n+\n".$trimmed_qual1."\n";
                        print OUT2 $curr2->header."\n".$trimmed_seq2."\n+\n".$trimmed_qual2."\n";
                    }
                }
            }
        }
    }else{
        print STDERR $curr1->header."\n".$curr1->seq."\n+\n".$curr1->qual."\n" if($verbose);
        print STDERR $curr2->header."\n".$curr2->seq."\n+\n".$curr2->qual."\n" if($verbose);
    
    }
    $i++;
}

exit;
