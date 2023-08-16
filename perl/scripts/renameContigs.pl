#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Basename;
#use Cwd;
#use File::Spec;
use Path::Class qw(dir);

my $usage=<<'ENDHERE';
NAME:
renameContigs.pl

PURPOSE:

INPUT:
--infiles <string> : list of fasta files separated by a ,
				
OUTPUT:
Will write a new file called prefix_renamed.fasta in the same
directeory as the input file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles);
my $verbose = 0;

GetOptions(
   'infiles=s' => \$infiles,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @infiles = split(/,/, $infiles);

foreach(@infiles){
   my $dirname = dirname($_);
   my @levels = split(/\//, $_);
   my $dir = shift(@levels); $dir = shift(@levels); $dir = shift(@levels);
   my $currSampleName = $dir;
   #$currSampleName =~ s{.*/}{};      # removes path  
   #$currSampleName =~ s{\.[^.]+$}{}; # removes extension
   my $newFileName = $dirname."/Contigs_renamed.fasta";

   open(OUT, ">".$newFileName) or die "Can't open $newFileName\n";

   my $ref_fasta_db = Iterator::FastaDb->new($_) or die("Unable to open Fasta file, $_\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my $header = $curr->header;
      my $seq = $curr->seq;
      my $replacement = $currSampleName."_contig";
      
      $header =~ s/contig/$replacement/;
      my $formattedSeq = "";
      while (my $chunk = substr($seq, 0, 150, "")) {
         $formattedSeq .= "$chunk\n";
      }
      print OUT $header."\n".$formattedSeq."\n";
   }
   close(OUT);
}  

