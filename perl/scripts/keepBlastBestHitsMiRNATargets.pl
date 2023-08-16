#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
parseMiRNATargets.pl

PURPOSE:
From the 9C output of ssearch36, distribute miRNA into their site types. According to 
scheme found in PMC3499661

INPUT:
--infile <string>                          : blast out file
--infile_fasta                             : fasta file (for computing query lengths).
--max_mismatch <int>                       : max number of mismatch(es).
--max_number_of_consecutive_mismatch <int> : max number of consecutive mismatch tolerated.
--max_number_of_gaps <int>                 : max number of gaps.
--min_seed_length <int>                    : min seed length 

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $infile_fasta, $max_mismatch, $max_number_of_consecutive_mismatch, $max_number_of_gaps);
my $verbose = 0;

GetOptions(
  'infile=s'                             => \$infile,
  'infile_fasta=s'                       => \$infile_fasta,
  'max_number_of_gaps=i'                 => \$max_number_of_gaps,
  'max_mismatch=i'                       => \$max_mismatch,
  'max_number_of_consecutive_mismatch=i' => \$max_number_of_consecutive_mismatch,
  'verbose'                              => \$verbose,
  'help'                                 => \$help
);
if ($help) { print $usage; exit; }

die("--infile missing...\n") unless($infile);
die("--infile_fasta missing...\n") unless($infile_fasta);

## MAIN
my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/^>//;
   $hash{$header} = length($curr->seq);
}

my %hash_res;
open(IN, "<".$infile) or die "Can't open $infile\n";
print STDERR "Processing $infile\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#/){
        next;
    }

    # using blastn out 8 fmt
    
    # for ssearch. cigar string = 13 col.
    my @row = split(/\t/, $_);
    my $curr_id       = $row[0];
    my $curr_target   = $row[1];
    my $curr_percid   = $row[2];
    my $curr_length   = $row[3];
    my $curr_mismatch = $row[4];
    my $curr_gap      = $row[5];
    my $curr_evalue   = $row[10];
    my $curr_cigar    = $row[12];

    my $passed = "false";
    if(exists($hash{$curr_id})){
        my $qlength = $hash{$curr_id};
        if($curr_length == $qlength){
            my @matches = ( $curr_cigar =~ /(\d+)M/g );
            my @insertions = ( $curr_cigar =~ /(\d+)I/g );
        
            #my $sum_matches = 0;
            #foreach my $match (@matches){
            #   $sum_matches = $sum_matches + $match;
            #}
        
            #my $sum_insertions = 0;
            #foreach my $insertion (@insertions){
            #    if($insertion <= $max_number_of_consecutive_mismatch){
            #        $passed = "true";
            #    }
            #    $sum_insertions = $sum_insertions + $insertion;
            #}
            if($curr_mismatch <= $max_mismatch){# && @insertions <= $max_number_of_gaps && $passed eq "true"){
                #print STDOUT $_."\n";
                $hash_res{$curr_id}{$curr_target} = $_;
            }
        }
    }
}
close(IN);

# Finally sort and print hash.
foreach my $key (keys %hash_res) {
    foreach my $key2 (keys %{ $hash_res{$key} }) {
        print STDOUT $hash_res{$key}{$key2}."\n";
    }
}
exit;
