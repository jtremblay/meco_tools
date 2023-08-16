#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
convertMetagenemarkNames.pl

PURPOSE:

INPUT:
--gff <string>         : gff infile
--fna <string>         : fasta nucl. file

OUTPUT:
STDOUT                 : renamed fna seqs.
--renamed_gff <string> : renamed gff.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $gff, $fna, $faa, $renamed_gff, $renamed_faa);
my $verbose = 0;

GetOptions(
   'gff=s'         => \$gff,
   'fna=s'         => \$fna,
   'faa=s'         => \$faa,
   'renamed_gff=s' => \$renamed_gff,
   'renamed_faa=s' => \$renamed_faa,
   'verbose'       => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(OUT_GFF, ">".$renamed_gff) or die "Can't open $renamed_gff\n";
open(OUT_FAA, ">".$renamed_faa) or die "Can't open $renamed_faa\n";

my %hash;
my %hash_geneid_old_to_new;
my $counter = 1;
open(GFF, "<".$gff) or die "Can't open $gff\n";
while(<GFF>){
   chomp;
   if($_ =~ m/^#/){
      print OUT_GFF $_."\n";
      next;
   }
   next if($_ =~ m/^$/);
   # Replace all " characters with ' chars. HtSeq has a problem with " chars...
   $_ =~ s/\"/\'/g;
   my @row = split(/\t/, $_);
   my $originalId = $row[0];
   my @field = split(/;/, $row[8]);
   my $mgmId = $field[0];
   $mgmId =~ s/ID=/gene_id=/;
   $mgmId =~ s/=/_/;
   $mgmId =~ s/gene_id_\d+_\d+/gene_id_$counter/;
   $originalId =~ s/\s/\./g;
   $hash{$mgmId} = $originalId."=".$mgmId;
   $hash_geneid_old_to_new{$field[0]} = $mgmId;

   # Then modify only gene_id in the gff output
   my $new_line = $_;
   $new_line =~ s/ID=/gene_id_/;
   $new_line =~ s/gene_id_\d+_\d+/gene_id_$counter/;
   print OUT_GFF $new_line."\n";

   $counter++;
}
close(GFF);
close(OUT_GFF);
#print STDERR Dumper(\%hash);
#exit;

my $ref_fasta_db = Iterator::FastaDb->new($fna) or die("Unable to open Fasta file, $fna\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   my ($contig_id) = $header =~ m/^(k\d+_\d\+)_\d\+ /;
   my ($old_gene_id) = $header =~ m/(ID=\d+_\d+);/;
   my $gene_id = $hash_geneid_old_to_new{$old_gene_id};
   if(exists $hash{$gene_id}){
      #print STDOUT ">".$hash{$header}."\n".$curr->seq."\n";
      #print STDOUT ">".$contig_id."=".$gene_id."\n".$curr->seq."\n";
      print STDOUT ">".$gene_id."\n".$curr->seq."\n";
   }
}

$ref_fasta_db = Iterator::FastaDb->new($faa) or die("Unable to open Fasta file, $faa\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   #my ($contig_id) = $header =~ m/^(\S+) /;
   my ($contig_id) = $header =~ m/^(k\d+_\d\+)_\d\+ /;
   my ($old_gene_id) = $header =~ m/(ID=\d+_\d+);/;
   my $gene_id = $hash_geneid_old_to_new{$old_gene_id};
   if(exists $hash{$gene_id}){
      #print OUT_FAA ">".$hash{$header}."\n".$curr->seq."\n";
      #print OUT_FAA ">".$contig_id."=".$gene_id."\n".$curr->seq."\n";
      print OUT_FAA ">".$gene_id."\n".$curr->seq."\n";
   }
}
close(OUT_FAA);
