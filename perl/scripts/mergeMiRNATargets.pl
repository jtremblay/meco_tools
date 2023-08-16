#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
mergeMiRNATargets.pl

PURPOSE:
merge to output of parseMiRNATargets.pl both fwd and rev into a 
single file.

INPUT:
--infile_fwd <string>  : infile fwd hits
--infile_rev <string>  : infile rev hits

OUTPUT:
STDOUT <string>        : standard output

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fwd, $infile_rev);
my $verbose = 0;

GetOptions(
  'infile_fwd=s' => \$infile_fwd,
  'infile_rev=s' => \$infile_rev,
  'verbose'      => \$verbose,
  'help'         => \$help
);
if ($help) { print $usage; exit; }

die("--infile_fwd missing...\n") unless($infile_fwd);
die("--infile_rev missing...\n") unless($infile_rev);

## MAIN
my %hash;
my $header_counter = 0;
open(IN, "<".$infile_fwd) or die "Can't open $infile_fwd\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#/ && my $header_counter == 0){
        print STDOUT $_."\n";
        $header_counter++;
        next;
    }elsif($_ =~ m/^#/){
        next;
    }
    my @row_orig = split(/\t/, $_);
    my @row = split(/\t/, $_);
    my $query = shift(@row);
    my $subject = shift(@row);
    #$hash{$query}{$subject}{fwd} = join("\t", @row_orig)."\t+";
    my $strand = $row[8];
    if($strand ne "+"){
        die("Strand value (column 10) should be '+'. Found : $strand\n");
    }
    print STDOUT join("\t", @row_orig)."\n";
}
close(IN);

#0                          1                       2                       3                       4                       5       6       7       8       9               10
#query_id	                subject_id	            match_aln	            query_aln	            subject_aln	            q_start	q_end	s_start	s_end	expect_value	strand
#ath-miR5645b_MIMAT0022409	649633106_contig_id_1	:::::. ::::::::::.. 	AUUUGAGUCAUGUCGUUAAG	AUUUGGCUCAUGUCGUUGGU	1	    19	    4239020	4239039	 3	            - 
my $rep = 1;
open(IN, "<".$infile_rev) or die "Can't open $infile_rev\n";
while(<IN>){
    next if($_ =~ m/^#/);
    chomp;
    #my @row_orig = split(/\t/, $_);
    my @row = split(/\t/, $_);
    #my $query = shift(@row);
    #my $subject = shift(@row);
    #my $aln_str = reverse(shift(@row));
    #my $query_str = reverse(shift(@row));
    #my $subject_str = reverse(shift(@row));
    #my $end = shift(@row);
    #my $start = shift(@row);
    #my $evalue = shift(@row);
    
    my $query       = $row[0];
    my $subject     = $row[1];
    my $aln_str     = reverse($row[2]);
    my $query_str   = reverse($row[3]);
    my $subject_str = reverse($row[4]);
    my $qstart      = $row[5];
    my $qend        = $row[6];
    my $sstart      = $row[7];
    my $send        = $row[8];
    my $evalue      = $row[9];
    my $strand      = $row[10];
    if($strand ne "-"){
        die("Strand value (column 10) should be '-'. Found : $strand\n");
    }
    #$hash{$query}{$subject}{rev}{1}{entry} = "$query\t$subject\t$aln_str\t$query_str\t$subject_str\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$strand";
    print STDOUT "$query\t$subject\t$aln_str\t$query_str\t$subject_str\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$strand\n";

}
close(IN);

foreach my $query (keys %hash) {
    foreach my $subject (keys %{ $hash{$query} }) {
        foreach my $strand (keys %{ $hash{$query}{$subject} }) {
            my $line = $hash{$query}{$subject}->{$strand};
            print STDOUT $line."\n";
        }
    }
}
