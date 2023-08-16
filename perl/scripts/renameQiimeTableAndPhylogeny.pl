#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile_fasta <string>     : Sequence file fasta fmt.
--infile_otu_table <string> : OTU table file (qiime2 output).
				
OUTPUT:
--outfile_fasta <string>    : Sequence file fasta.
STDOUT <string>             : otu_table file qiime legacy format (tsv/txt).

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_otu_table, $infile_fasta, $outfile_fasta);
my $verbose = 0;

GetOptions(
   'infile_otu_table=s' => \$infile_otu_table,
   'infile_fasta=s'     => \$infile_fasta,
   'outfile_fasta=s'    => \$outfile_fasta,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my $counter = 1;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fastq file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/^>//;
   $hash{$header}{seq} = $curr->seq;
   $hash{$header}{id} = $counter;
   $counter++;
}

open(IN, "<".$infile_otu_table) or die "Can't open $infile_otu_table\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#/){
        if($_ =~ m/^#OTU ID/){
            print STDOUT $_."\n";
        }

        next;
    }

    my @row = split(/\t/, $_);
    my $id = shift(@row);
    my $lineage = pop(@row);
   
    if(exists $hash{$id}){
        print STDOUT $hash{$id}{id}."\t"; 
    }else{
        print STDERR "Warning: $id does not exists in fasta file\n";
    }

    print STDOUT join("\t", @row)."\t";

    $lineage =~ s/\s+//g;
    $lineage =~ s/D_0__/k__/;
    $lineage =~ s/D_1__/p__/;
    $lineage =~ s/D_2__/c__/;
    $lineage =~ s/D_3__/o__/;
    $lineage =~ s/D_4__/f__/;
    $lineage =~ s/D_5__/g__/;
    
    $lineage =~ s/D_6__;//;
    $lineage =~ s/D_7__;//;
    $lineage =~ s/D_8__;//;
    $lineage =~ s/D_9__;//;
    $lineage =~ s/D_10__;//;
    $lineage =~ s/D_11__;//;
    $lineage =~ s/D_12__;//;
    $lineage =~ s/D_13__;//;
    $lineage =~ s/D_14__;//;
    
    print STDOUT $lineage."\n";
    
}
close(IN);

open(OUT, ">".$outfile_fasta) or die "Can't open $outfile_fasta\n";
for my $key (keys %hash){
    print OUT ">".$hash{$key}{id}."\n".$hash{$key}{seq}."\n";
}
close(OUT);
exit;
