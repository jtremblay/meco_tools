#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
getAdaptersStats.pl

PURPOSE:

INPUT:
--reads_1 <string>  : Fastq sequence file.
--reads_2 <string>  : Fastq sequence file.
--adapters <striing : Fasta sequence file. Containing only two entries.
        
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $adapters, $reads_1, $reads_2);
my $verbose = 0;

GetOptions(
  'reads_1=s'   => \$reads_1,
  'reads_2=s'   => \$reads_2,
  'adapters=s'  => \$adapters,
  'verbose'     => \$verbose,
  'help'        => \$help
);
if ($help) { print $usage; exit; }

die "--reads_1 missing\n" unless($reads_1);
die "--reads_2 missing\n" unless($reads_2);
die "--adapters missing\n" unless($adapters);

## MAIN
my %hash;
my $counter = 1;
my $ref_fasta_db = Iterator::FastaDb->new($adapters) or die("Unable to open Fasta file, $adapters\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  $hash{$counter} = $curr->seq;
  $counter++;
}
die "Adapters fasta file should contain only 2 entries. Contains $counter\n" if($counter > 3);

my %hashPos;
$hashPos{reads1}{count}=0;
$hashPos{reads1}{position}="";
$hashPos{reads2}{count}=0;
$hashPos{reads2}{position}="";
my $ref_fastq_db = Iterator::FastqDb->new($reads_1) or die("Unable to open Fastq file, $reads_1\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
  my $seq = $curr->seq;
  #my $pattern = uc($hash{1});
  my $pattern = uc(substr($hash{1}, 19));
  if($seq =~ m/($pattern)/g){
    my $pos = pos($seq) - 1;
    $hashPos{reads1}{count}++;
    $hashPos{reads1}{position} .= $pos."\t";
  }
}

print STDERR "Sequence to search for in reads 2: ".uc(substr($hash{2}, 19))."\n";
#$ref_fastq_db = undef;
my $ref_fastq_db2 = Iterator::FastqDb->new($reads_2) or die("Unable to open Fastq file, $reads_2\n");
while( my $curr = $ref_fastq_db2->next_seq() ) {
  my $seq = $curr->seq;
  my $pattern = uc(substr($hash{2}, 20));
  #print STDERR "$seq\n$pattern\n";
  #$pattern =~ tr/ACGT/TGCA/; #Just convert to complement since the reads will be reversed
  #$pattern = reverse($pattern);
  if($seq =~ m/($pattern)/g){
    my $pos = pos($seq) - 1;
    $hashPos{reads2}{count}++;
    $hashPos{reads2}{position} .= $pos."\t";
    #print STDERR "Found!!\n";
  }
}

for my $k1 (sort keys %hashPos) {
  print STDOUT $k1.":\n";
  for my $k2 (sort keys %{ $hashPos{$k1} }){
    print STDOUT $k2.":\n";
    print $hashPos{$k1}{$k2}."\n";

  }
}

exit;

