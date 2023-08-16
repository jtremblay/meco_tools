#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateGeneCallingFileForAnvio.pl

PURPOSE:

INPUT:
--infile_gff <string>      : Gff output file from Prodigal (the parsed one with nrc_tools).
--infile_faa <string>      : Faa output file from Prodigal (the parsed one with nrc_tools).
--infile_rrna_bac <string> : barrnap gff3 output

OUTPUT:
<STDOUT>                   : tsv file compatible with anvio.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $gff, $faa, $rrna_bac);
my $verbose = 0;

GetOptions(
   'infile_gff=s'      => \$gff,
   'infile_faa=s'      => \$faa,
   'infile_rrna_bac=s' => \$rrna_bac,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my $ref_fasta_db = Iterator::FastaDb->new($faa) or die("Unable to open Fasta file, $faa\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;
   $header =~ s/>//;
   my $gene_id = $header;
   $gene_id =~ s/gene_id_//;
   my $seq = $curr->seq;
   $seq =~ s/\*$//;
   $hash{$gene_id} = $seq;
}

my $last_gene_id;
print STDOUT "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n";
#print STDOUT "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n";
open(GFF, "<".$gff) or die "Can't open $gff\n";
while(<GFF>){
   chomp;
   if($_ =~ m/^#/){
      next;
   }
   next if($_ =~ m/^$/);
   # Replace all " characters with ' chars. HtSeq has a problem with " chars...
   $_ =~ s/\"/\'/g;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my @field = split(/;/, $row[8]);
   my $gene_id = $field[0];
   $gene_id =~ s/gene_id_//;
   my $start = $row[3] - 1;
   my $end = $row[4] + 0;
   my $direction = $row[6];
   if($direction eq "+"){ 
       $direction = "f"; 
   }elsif($direction eq "-"){ 
       $direction = "r"; 
   }else{
       die "cant figure direction ".$direction;
   }
   my $partial = $field[1];
   $partial =~ s/partial\=//;
   my $call_type = 1;
   my $soft_ver = $row[1];
   my @soft_ver = split(/_/, $soft_ver);
   my $source = $soft_ver[0];
   my $version = $soft_ver[1];
    
   print STDOUT $gene_id."\t".$contig_id."\t".$start."\t".$end."\t".$direction."\t".$partial."\t"."1"."\t".lc($source)."\t".$version."\t".$hash{$gene_id}."\n";
   #print STDOUT $gene_id."\t".$contig_id."\t".$start."\t".$end."\t".$direction."\t".$partial."\t".$source."\t".$version."\n";
   $last_gene_id = $gene_id;

}
close(GFF);

open(GFF, "<".$rrna_bac ) or die "Can't open $rrna_bac\n";
while(<GFF>){
   chomp;
   if($_ =~ m/^#/){
      next;
   }
   next if($_ =~ m/^$/);
   # Replace all " characters with ' chars. HtSeq has a problem with " chars...
   $_ =~ s/\"/\'/g;
   
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   #my $gene_id = "";
   #if($_ =~ m/Name=(\S+);.*/){
   #   $gene_id = $contig_id."_".$1;
   #}

   my @field = split(/;/, $row[8]);
   my $start = $row[3] - 1;
   my $end = $row[4] + 0;
   my $direction = $row[6];
   if($direction eq "+"){ 
       $direction = "f"; 
   }elsif($direction eq "-"){ 
       $direction = "r"; 
   }else{
       die "cant figure direction ".$direction;
   }
   my $partial = "00";
   if($row[8] =~ m/partial/){
      $partial = "10"; #not optimal, but will do for now.
   }
   my $call_type = 0;
   my $soft_ver = $row[1];
   my @soft_ver = split(/:/, $soft_ver);
   my $source = $soft_ver[0];
   my $version = $soft_ver[1];
    
   $last_gene_id++;
   print STDOUT $last_gene_id."\t".$contig_id."\t".$start."\t".$end."\t".$direction."\t".$partial."\t"."2"."\t".lc($source)."\t".$version."\t".""."\n";
   #print STDOUT $gene_id."\t".$contig_id."\t".$start."\t".$end."\t".$direction."\t".$partial."\t".$source."\t".$version."\n";

}
close(GFF);



