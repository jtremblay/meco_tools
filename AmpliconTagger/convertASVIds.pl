#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
convertASVIds.pl

PURPOSE:

INPUT:
--infile_tsv <string>    : Sequence file
--min_reads <int>        : Optional. Will reject OTUs for which sum of reads across
                           all samples is lower than <int>.
				
OUTPUT:
--outfile_tsv <string>   : Sequence file
--outfile_fasta <string> : Sequence file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_tsv, $outfile_tsv, $outfile_fasta, $min_reads);
my $verbose = 0;

GetOptions(
   'infile_tsv=s'    => \$infile_tsv,
   'min_reads=i'     => \$min_reads,
   'outfile_tsv=s'   => \$outfile_tsv,
   'outfile_fasta=s' => \$outfile_fasta,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }
$min_reads = 0 unless($min_reads);
## MAIN
open(OUT_FNA, ">".$outfile_fasta) or die "Can't open $outfile_fasta file\n";
open(OUT_TSV, ">".$outfile_tsv) or die "Can't open $outfile_tsv file\n";

my $counter = 1;
open(IN, "<".$infile_tsv) or die "Can't open $infile_tsv\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#FEATURE_ID/){
        my $curr_line = $_;
        $curr_line =~ s/^#FEATURE_ID/#CLUSTER/;
        print OUT_TSV $curr_line."\n";
    }elsif($_ =~ m/^#/){
        next; 
    }else{
        my @row = split(/\t/, $_);
        my $id = shift(@row);
        if($min_reads > 0){
            my $total = 0;
            foreach my $el (@row){
                $total += $el;
            }
            if($total >= $min_reads){
                print OUT_TSV $counter."\t".join("\t", @row)."\n";
                print OUT_FNA ">".$counter.";size=$total\n".$id."\n";
            }

        }else{
            my $total = 0;
            foreach my $el (@row){
                $total += $el;
            }
            print OUT_TSV $counter."\t".join("\t", @row)."\n";
            print OUT_FNA ">".$counter.";size=$total\n".$id."\n";
        
        }
        $counter++;
    }
}
close(IN);
close(OUT_TSV);
close(OUT_FNA);
