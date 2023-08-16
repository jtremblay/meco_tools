#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
barcodesSorted.pl

PURPOSE:

INPUT:
--infile_fastq <string>    : Fastq file (file to demultiplex)
--infile_barcodes <string> : Fasta file

OUTPUT:
--outdir <string>          : Directory where will be written fastqs

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fastq, $infile_barcodes, $outdir, $log);
my $verbose = 0;

GetOptions(
   'infile_fastq=s'    => \$infile_fastq,
   'infile_barcodes=s' => \$infile_barcodes,
   'outdir=s'          => \$outdir,
   'log=s'             => \$log,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

system("mkdir -p $outdir");

## MAIN
my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($infile_barcodes) or die("Unable to open Fasta file, $infile_barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   my $seq = $curr->seq;
   $hash{$seq} = $header;
}

my %hash_count;

my $prev_barcode = "";
my $fastq_out_fn = $outdir."/dummy.fastq";
open my $FASTQ_OUT, '>', $fastq_out_fn;
my $ref_fastq_db = Iterator::FastqDb->new($infile_fastq) or die("Unable to open Fastq file, $infile_fastq\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
    my $barcode = $curr->barcode;
    if($barcode ne $prev_barcode){
        close($FASTQ_OUT);
        $fastq_out_fn = $outdir."/".$hash{$barcode}.".fastq";
        $hash_count{$barcode}++;
        open $FASTQ_OUT, '>', $fastq_out_fn;
        print $FASTQ_OUT $curr->output;
        $prev_barcode = $barcode;
    }elsif($barcode eq $prev_barcode){
        print $FASTQ_OUT $curr->output;
    }
}

if($log){
    open(LOG, ">".$log) or die "Can't open $log\n";
    # Dump barcodes counts.
    for my $key (sort keys %hash_count){
        print LOG $key."\t".$hash_count{$key}."\n";
    }
    close(LOG);
}
exit;
