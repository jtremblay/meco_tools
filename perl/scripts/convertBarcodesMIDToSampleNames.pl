#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
#use slice;

my $usage=<<'ENDHERE';
NAME:
convertBarcodesMIDToSampleNames.pl

PURPOSE:

INPUT:
--barcodes <string>     : Sequence file
--mappingFile <string> : mapping file in tsv
--col_mid <int>         : column where MID is
--col_sampleName <int>  : column where sample name is.

OUTPUT:
STDOUT                  : barcodes
--mappingFileOut      : mapping file with sample names in first row.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $barcodes, $mappingFile, $col_mid, $col_sampleName, $mappingFileOut);
my $verbose = 0;

GetOptions(
   'barcodes=s'         => \$barcodes,
   'mappingFile=s'      => \$mappingFile,
   'col_mid=i'          => \$col_mid,
   'col_sampleName=i'   => \$col_sampleName,
   'mappingFileOut=s'   => \$mappingFileOut,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

open(MAPOUT, ">".$mappingFileOut) or die "Can't open $mappingFileOut\n";

open(MAP, "<".$mappingFile) or die "Can't open $mappingFile\n";
while(<MAP>){
    chomp;
    my @row = split(/\t/, $_);
    $hash{"MID".$row[$col_mid]} = $row[$col_sampleName];
  
    #my @out = splice @row, 2

    if($. == 1){
        print MAPOUT "#".$row[$col_sampleName]."\t".$row[$col_mid]."\t".join("\t", @row[2..(@row-1)])."\n";
     }else{
        print MAPOUT $row[$col_sampleName]."\t".$row[$col_mid]."\t".join("\t", @row[2..(@row-1)])."\n";
    }
}
close(MAP);
close(MAPOUT);

my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $mid = $curr->header;
    $mid =~ s/>//;
    if(exists $hash{$mid}){
        print STDOUT ">".$hash{$mid}."\n".$curr->seq."\n";
    }else{
        print STDERR "$mid is missing in barcodes file\n";  
    }
}

