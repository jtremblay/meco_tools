#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
splitRnammer.pl

PURPOSE:

INPUT:
--infile_tab <string>   : Sequence file (output of hmmscan tblout against rrna models).
--infile_fasta <string> : Sequence file of contigs. rrna seqs will be extracted from this file using coords in --infile_tab
				
OUTPUT:
--S5 <string>   : ssu 5S
--S8 <string>   : ssu 5S
--S16 <string>  : ssu 16S
--S18 <string>  : ssu 18S
--S23 <string>  : lsu 23S
--S28 <string>  : lsu 28S

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_tab, $infile_fasta, $arc_5S, $arc_16S, $arc_23S, $bac_5S, $bac_16S, $bac_23S, $euk_5S, $euk_18S, $euk_28S);
my $verbose = 0;

GetOptions(
   'infile_tab=s'   => \$infile_tab,
   'infile_fasta=s' => \$infile_fasta,
   'arc_5S=s'       => \$arc_5S,
   'arc_16S=s'      => \$arc_16S,
   'arc_23S=s'      => \$arc_23S,
   'bac_5S=s'       => \$bac_5S,
   'bac_16S=s'      => \$bac_16S,
   'bac_23S=s'      => \$bac_23S,
   'euk_5S=s'       => \$euk_5S,
   'euk_18S=s'      => \$euk_18S,
   'euk_28S=s'      => \$euk_28S,
   'verbose' 	     => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT_arc_5S,  ">".$arc_5S)  or die "Can't open $arc_5S\n";
open(OUT_arc_16S, ">".$arc_16S) or die "Can't open $arc_16S\n";
open(OUT_arc_23S, ">".$arc_23S) or die "Can't open $arc_23S\n";
open(OUT_bac_5S,  ">".$bac_5S)  or die "Can't open $bac_5S\n";
open(OUT_bac_16S, ">".$bac_16S) or die "Can't open $bac_16S\n";
open(OUT_bac_23S, ">".$bac_23S) or die "Can't open $bac_23S\n";
open(OUT_euk_5S,  ">".$euk_5S)  or die "Can't open $euk_5S\n";
open(OUT_euk_18S, ">".$euk_18S) or die "Can't open $euk_18S\n";
open(OUT_euk_28S, ">".$euk_28S) or die "Can't open $euk_28S\n";

my %hash;
open(IN, "<".$infile_tab) or die "Can't open $infile_tab\n";
while(<IN>){
   chomp;
   if($_ =~ m/^#/){
      #$header .= $_."\n" if($k == 0);
      next;
   }

   my @row = split(/\s+/, $_);

   my $target_name = $row[0];
   my $query_name = $row[3];
   my $qlen = $row[5];
   my $evalue = $row[6];
   my $gene_start = $row[17];
   my $gene_end = $row[18];

   #next if($evalue >= $e);
   #next if($qlen < $ql);

   #if(exists $hash{$query_name} && $hash{$query_name}{contig_id} eq $query_name && $hash{$query_name}{target_name} eq $target_name && $evalue < $hash{$query_name}{evalue}){
   #if(exists $hash{$query_name} && $evalue < $hash{$query_name}{evalue}){
   if(exists $hash{$query_name}){
       #$hash{$query_name}{contig_id} = $query_name;
       #     $hash{$query_name}{target_name} = $target_name;
       #     $hash{$query_name}{start} = $gene_start;
       #     $hash{$query_name}{end} = $gene_start;
       #     $hash{$query_name}{evalue} = $evalue;
   }else{
      $hash{$query_name}{contig_id} = $query_name;
      $hash{$query_name}{target_name} = $target_name;
      $hash{$query_name}{start} = $gene_start;
      $hash{$query_name}{end} = $gene_end;
      $hash{$query_name}{evalue} = $evalue;
   }  
}
close(IN);

#print STDERR Dumper(\%hash);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my($contig_id) = $curr->header =~ m/^>(\S+) /;
   #print (STDERR "contig_id:$contig_id\n");
   if(exists $hash{$contig_id}){
      my $rrna_seq;
      my $curr_target_name = $hash{$contig_id}{target_name};
      my $gene_start = $hash{$contig_id}{start};
      my $gene_end = $hash{$contig_id}{end};

      if($gene_start < $gene_end){
         $rrna_seq = substr($curr->seq, $gene_start, abs($gene_end - $gene_start));
      }else{
         $rrna_seq = substr($curr->seq, $gene_end, abs($gene_start - $gene_end));
      }

      if($curr_target_name =~ m/arc_tsu/){
         print OUT_arc_5S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/bac_tsu/){
         print OUT_bac_5S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/euk_tsu/){
         print OUT_euk_5S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/arc_ssu/){
         print OUT_arc_16S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/bac_ssu/){
         print OUT_bac_16S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/euk_ssu/){
         print OUT_euk_18S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/arc_lsu/){
         print OUT_arc_23S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/bac_lsu/){
         print OUT_bac_23S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }elsif($curr_target_name =~ m/euk_lsu/){
         print OUT_euk_28S ">$contig_id $curr_target_name coords:$gene_start-$gene_end\n$rrna_seq\n";
      }

   }
}
close(OUT_arc_5S);
close(OUT_arc_16S);
close(OUT_arc_23S);
close(OUT_bac_5S); 
close(OUT_bac_16S);
close(OUT_bac_23S);
close(OUT_euk_5S); 
close(OUT_euk_18S);
close(OUT_euk_28S);

exit;
