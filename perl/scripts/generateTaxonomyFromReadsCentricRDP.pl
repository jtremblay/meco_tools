#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Env qw/TMPDIR/;
use Iterator::FastaDb;
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
generateReadsCentricOTUTable.pl

PURPOSE:
Split interleaved fastq in two files: one for R1 and the other for R2.

INPUT:
--infile_rdp <string>      :  Sequence file
--infile_barcodes <string> :  Sequence file

OUTPUT:
STDOUT <string>            : OTU table format file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile_rdp, $infile_barcodes, $cutoff);

my $verbose = 0;

GetOptions(
  'infile_rdp=s'      => \$infile_rdp,
  'infile_barcodes=s' => \$infile_barcodes,
  'cutoff=f'          => \$cutoff,
  'verbose'           => \$verbose,
  'help'              => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile_rdp);
die("--infile arg required\n") unless($infile_barcodes);
$cutoff = 0.5 unless($cutoff);

my %hash;
my %hash_barcodes;
my $ref_fasta_db = Iterator::FastaDb->new($infile_barcodes) or die("Unable to open Fasta file, $infile_barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  my $header = $curr->header;
  $header =~ s/>//;
  $header =~ s/\//-/; # / will interfere in the moment of writing into file.
  $hash_barcodes{$curr->seq} = $header;
}

open(IN, "<".$infile_rdp) or die "Can't open $infile_rdp\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my ($id) = $row[0] =~ m/#([ACGT]*)/;
   # 5,6,7 = k__bact, Kingdom, score
   # 8,9,10 = phy
   # 11,12,13 = class
   # 14,15,16 = order
   # 17,18,19 = family
   # 20,21,22 = genus
   # 23,24,25 = specie
  
   #print STDERR $id."\n";

   my $lineage = "";
   if($row[7] >= $cutoff){
      $lineage .= "$row[5]";
   }else{
      $lineage .= "";
   }

   if($row[10] >= $cutoff){
      $lineage .= "\t$row[8]";
   }else{
      $lineage .= "\t";
   }

   if($row[13] >= $cutoff){
      $lineage .= "\t$row[11]";
   }else{
      $lineage .= "\t";
   }

   if($row[16] >= $cutoff){
      $lineage .= "\t$row[14]";
   }else{
      $lineage .= "\t";
   }

   if($row[19] >= $cutoff){
      $lineage .= "\t$row[17]";
   }else{
      $lineage .= "\t";
   }
   
   if($row[22] >= $cutoff){
      $lineage .= "\t$row[20]";
   }else{
      $lineage .= "\t";
   }

   #if($row[25] >= $cutoff){
   #   $lineage .= "\t$row[23]";
   #}else{
   #   $lineage .= "\t";
   #}
   next if $lineage eq "";

   $hash{$lineage}{$id}++;

}
print STDERR Dumper(\%hash);
print STDERR Dumper(\%hash_barcodes);

print STDOUT "taxonomy";
for my $id (sort {$hash_barcodes{$a} cmp $hash_barcodes{$a}} keys %hash_barcodes){
    print STDOUT "\t".$hash_barcodes{$id};
}
print STDOUT "\n";

for my $key (keys %hash){
    print STDOUT $key."";

    for my $id (sort {$hash_barcodes{$a} cmp $hash_barcodes{$a}} keys %hash_barcodes){
        if(exists $hash{$key}{$id}){
            print STDOUT "\t".$hash{$key}{$id};
        }else{
            print STDOUT "\t0";
        }
    }
    print STDOUT "\n";

}

exit;
