#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--interleaved <string> : Interleaved fastq file.
OR
--reads1 <string>      : forward reads.
--reads2 <string>      : reverse reads.
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $reads1, $reads2);
my $verbose = 0;

GetOptions(
   'interleaved=s' => \$infile,
   'reads1=s' 	    => \$reads1,
   'reads2=s' 	    => \$reads2,
   'verbose' 	    => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

my $tmpdir = File::Temp->newdir(
    "tmpdir-mergePairs-XXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

## MAIN
if($infile){

   my $pipe = "$tmpdir/reads.pipe";
   system("mkfifo $pipe");
   system("gunzip -c ".$infile." > $pipe &");

   my $counter = 0;
   my $previous_pair;
   my $previous_entry;
   my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fastq file, $pipe\n");
   while( my $curr = $ref_fastq_db->next_seq() ) {

      if($counter == 0){
         $previous_pair = 1;
         $previous_entry = $curr->output;

     }else{
         if($curr->pair == 1){
            if($previous_pair != 2){
               die("Potential with these entries: ".$curr->output."\nAND\n".$previous_entry."\n");
            }
         }elsif($curr->pair == 2){
            if($previous_pair != 1){
               die("Potential with these entries: ".$curr->output."\nAND\n".$previous_entry."\n");
            }
         } 
         $previous_pair = $curr->pair;
         $previous_entry = $curr->output;
      }

      $counter++;
   }
}
