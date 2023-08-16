#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
extractReadsFromClusters.pl

PURPOSE:
Get reads belonging to each cluster generated from AmpliconTagger's
two round clustering procedure.

INPUT:
--infile_uc1 <string>   : uc file (i.e. 97pc file)
--infile_uc2 <string>   : uc file (i.e 99pc file)
--infile_fasta <string> : fasta file (i.e. derep1.fasta file)	

OUTPUT:
--outdir <string>       : outdir containing fasta outfiles

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_uc1, $infile_uc2, $infile_fasta, $outdir);
my $verbose = 0;

GetOptions(
   'infile_uc1=s'   => \$infile_uc1,
   'infile_uc2=s'   => \$infile_uc2,
   'infile_fasta=s' => \$infile_fasta,
   'outdir=s'       => \$outdir,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash_seq;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/^>//; 
    my @row = split(/;/, $header);
    my $seq_id = $row[0];
    $hash_seq{$seq_id}{seq} = $curr->seq;
    $hash_seq{$seq_id}{orig_header} = $curr->header;
}


my %hash;
open(IN, "<".$infile_uc1) or die "Can't open $infile_uc1\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $source_id = $row[1];
   my $target_id = $row[8];
   my @target_id = split(/;/, $target_id);
   $target_id = $target_id[0];
   $hash{$target_id} = $source_id;
}
close(IN);
print STDERR Dumper(\%hash);

my %final_hash;
open(IN, "<".$infile_uc2) or die "Can't open $infile_uc2\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $source_id = $row[1];
   my $target_id = $row[8];
   my @target_id = split(/;/, $target_id);
   $target_id = $target_id[0];
   print STDERR "target_id: ".$target_id."\n";

   if(exists $hash{$source_id}){
      if(exists $hash{$target_id}){
         open(OUT, ">>".$outdir."/".$hash{$source_id}.".fasta") or die "Can't open ".$outdir."/".$hash{$source_id}.".fasta\n";
         print OUT ">".$hash_seq{$target_id}{orig_header}."\n".$hash_seq{$target_id}{seq}."\n";
         close(OUT);
      }
   }
}
close(IN)


