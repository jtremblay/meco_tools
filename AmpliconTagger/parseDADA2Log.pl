#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseDADA2Log.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

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

#open(IN, "<".$infile) or die "Can't open $infile\n";
print STDOUT "sample_id\ttotal_reads\tmerged_reads\tnumber_of_unique_pairings\n";
my %hash;
my $sample_name;
while(<>){
    chomp;
    
    if($_ =~ m/Processing: (.*)$/){
        $sample_name = $1;
        #print STDERR $sample_name."\n"; 
    }
    #55720 paired-reads (in 478 unique pairings) successfully merged out of 55880 (in 563 pairings) input.
    my $merged; my $number_of_pairings; my $total_reads;
    if($_ =~ m/^(\d+) paired-reads .*(\d+) unique pairings.* merged out of (\d+) \(.*$/){
        $merged = $1; 
        $number_of_pairings = $2;
        $total_reads = $3;
        $hash{$sample_name}{merged} = $merged;
        $hash{$sample_name}{number_of_pairings} = $number_of_pairings;
        $hash{$sample_name}{total_reads} = $total_reads;

    }
}
#close(IN);

for my $key (keys %hash){
    print STDOUT "$key\t$hash{$key}{total_reads}\t$hash{$key}{merged}\t$hash{$key}{number_of_pairings}\n";
}
exit;
