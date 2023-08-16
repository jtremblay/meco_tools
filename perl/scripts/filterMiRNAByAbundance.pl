#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
~/scripts/filterMiRNAByAbundance.pl

PURPOSE:

INPUT:
--infile_tsv <string>    : Sequence file
--infile_fasta <string>  : Infile fasta
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
my ($help, $infile_tsv, $outfile_tsv, $infile_fasta, $outfile_fasta, $min_reads);
my $verbose = 0;

GetOptions(
   'infile_tsv=s'      => \$infile_tsv,
   'infile_fasta=s'    => \$infile_fasta,
   'min_reads=i'       => \$min_reads,
   'outfile_tsv=s'     => \$outfile_tsv,
   'outfile_fasta=s'   => \$outfile_fasta,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }
$min_reads = 0 unless($min_reads);

## MAIN
open(OUT_FNA, ">".$outfile_fasta) or die "Can't open $outfile_fasta file\n";
open(OUT_TSV, ">".$outfile_tsv) or die "Can't open $outfile_tsv file\n";

#my $counter = 1;
my %hash;
open(IN, "<".$infile_tsv) or die "Can't open $infile_tsv\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#CLUSTER/){
        my $curr_line = $_;
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
                print OUT_TSV $id."\t".join("\t", @row)."\n";
                $hash{$id} = $hash{$id};
            }

        }else{ # no filtering
            #my $total = 0;
            #foreach my $el (@row){
            #    $total += $el;
            #}
            #print OUT_TSV $counter."\t".join("\t", @row)."\n";
            #print OUT_FNA ">".$counter.";size=$total\n".$id."\n";
            die "No filtering to be done...\n";
        
        }
        #$counter++;
    }
}
close(IN);
close(OUT_TSV);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/^>//;
    if(exists $hash{$header}){
        my $seq = $curr->seq;
        $seq =~ s/T/U/g;
        print OUT_FNA $curr->header."\n".$seq."\n";
    }
}
close(OUT_FNA);


exit;
