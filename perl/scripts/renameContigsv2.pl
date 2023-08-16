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
renameContigsv2.pl

PURPOSE:

INPUT:
--infiles <string> : list of fasta files separated by a ,
				
OUTPUT:
STDOUT

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

   print STDERR "Processing $_\n";

   my $curr_file = $_;
   my $name = $curr_file;
   $name =~ s{.*/}{};      # removes path  
   $name =~ s{\.[^.]+$}{}; # removes extension

   my $ref_fasta_db = Iterator::FastaDb->new($curr_file) or die("Unable to open Fasta file, $curr_file\n");
   my $i = 1;
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my $header = $curr->header;
      my $seq = $curr->seq;
      my $replacement = ">".$name."_contig-".$i;
      
      my $formattedSeq = "";
      while (my $chunk = substr($seq, 0, 150, "")) {
         chomp($chunk);
         $formattedSeq .= "$chunk\n";
      }
      print STDOUT $replacement."\n".$formattedSeq;
      $i++;
   }
}  

