#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Encode;


my $usage=<<'ENDHERE';
NAME:
parseNR.pl

PURPOSE:

INPUT:
--infile <string> : NCBI's nr database.
				
OUTPUT:
STDOUT            : Spreadsheet with only parsed headers. No sequences.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header;

   my @header;
   if($header =~ m//){
      @header = split(/\cA/, $header);
   }else{
      $header[0] = $header;
   }
   my $parsed_header = $header[0];
   #print STDERR $parsed_header."\n";
   
   #my ($gi) = $parsed_header =~ />(.*)\| (.*)$/;
   #>gi|916623119|ref|WP_051230210.1| hypothetical protein [Haliea salexigens]
   #if($parsed_header =~ />(.*\|) (.*)$/){
   if($parsed_header =~ />(\S+) (.*)$/){
      my $gi = $1;
      my $desc = $2;
      $desc =~ s/MULTISPECIES\://;

      #print STDERR $gi."\n";
      if($gi =~ m/(gi\|\d+\|)/){
         print STDOUT "$1\t$desc\n";
      }else{
         print STDERR "No GI found...\n";
      }
   }
   else{
      print STDERR "$parsed_header\n";
   }
}
exit;
