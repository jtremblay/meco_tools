#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile_barcodes <string>     : Sequence file
--infile_map <string>          : mapping file
				
OUTPUT:
--outfile_barcodes <string>    : barcode file STDOUT...

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_barcodes, $infile_map);
my $verbose = 0;

GetOptions(
   'infile_map=s' 	  => \$infile_map,
   'infile_barcodes=s' => \$infile_barcodes,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($infile_barcodes) or die("Unable to open Fasta file, $infile_barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   $hash{$header} = $curr->seq;
}

open(IN, "<".$infile_map) or die "Can't open $infile_map\n";

while(<IN>){
   next if($_ =~ m/^#/);
   chomp;
   my @row = split(/\t/, $_);
   my $old_sample_name = $row[0];
   my $new_sample_name = $row[1];
   
   if(exists $hash{$old_sample_name}){
      print STDOUT ">".$new_sample_name."\n".$hash{$old_sample_name}."\n";
      delete($hash{$old_sample_name});
   }
}
close(IN);

print STDERR Dumper(\%hash);


