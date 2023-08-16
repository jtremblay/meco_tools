#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
computeDepthBreadth.pl

PURPOSE:

INPUT:
--infile <string> : Output from samtools depth
				
OUTPUT:
<STDOUT>          : output by contigs

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

push @ARGV, undef unless @ARGV;    # READ FROM STDIN IF NOT INFILES PROVIDED

if(!defined($infile)){
    $infile = shift(@ARGV);
}

## MAIN
my %hash;

my $old_id = "";
my $curr_id = "";
my $cumm_values = 0;
my $num_bases = 0;
my $num_pos_gt_1x = 0;
print STDOUT "contig_id\tcontig_length\ttotal_bases_mapped\taverage_coverage\tnumber_of_nucl_pos_gt_1X\n";
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    $curr_id = $row[0];
    if($curr_id eq $old_id){
        $cumm_values += $row[2];
        $num_bases++;
        if($row[2] > 1){
            $num_pos_gt_1x++;
        }
    }else{
        if($. != 1){
            # compute if more than 50
            print STDOUT $old_id."\t".$num_bases."\t".$cumm_values."\t".sprintf("%.2f", ($cumm_values/$num_bases))."\t".$num_pos_gt_1x."\n";
        }
        # start a new contig or gene.
        $cumm_values = $row[2];
        $num_bases = 1;
        if($row[2] > 1){
            $num_pos_gt_1x = 1;
        }
        $old_id = $curr_id;
    }
}
