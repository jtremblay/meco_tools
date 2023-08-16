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
--infile <string> : list of fasta files separated by a ,
--prefix <string> : prefix to add to contigs header. >Contigs.<prefix>.n		

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
my ($help, $infile, $prefix);
my $verbose = 0;

GetOptions(
   'infile=s'  => \$infile,
   'prefix=s'  => \$prefix,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
my $i = 1;
while( my $curr = $ref_fasta_db->next_seq() ) {
    #my $header = $curr->header;
   my $seq = $curr->seq;
   my $header = ">Contig$prefix.$i";
   
   my $formattedSeq = "";
   while (my $chunk = substr($seq, 0, 150, "")) {
      $formattedSeq .= "$chunk\n";
   }
   print STDOUT $header."\n".$formattedSeq."\n";
   $i++;
}

