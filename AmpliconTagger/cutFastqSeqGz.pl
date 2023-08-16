#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Temp;
use Env qw(TMPDIR);

my $usage=<<'ENDHERE';
NAME:
cutFastqSeq.pl

PURPOSE:
To remove a specific number of nucleotides from either
the beginning or end of reads (or both). Useful if we 
have a library that consistently shows bad quality after 
nt 125..., or a library in which the last nucleotides 
are consistently of bad quality.
Also useful to trimm reads before assembly with FLASH.
INPUT:
--infile <string>  : fastq infile
--begin <int>      : Number of nucleotides to remove at the 5 prime region. 
                     Ex: if 3 is entered, the first 3 nucleotides will be removed
--end <int>        : Number of nucleotides to remove at the 3 prime region. 
                     Ex: if 3 is entered, the last 3 nucleotides will be removed
OUTPUT:
--outfile <string> : trimmed fastq file.
NOTES:

BUGS/LIMITATIONS:
TODO Parallelize this script.

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $begin, $end);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'infile=s'    => \$infile,
  'outfile=s'   => \$outfile,
  'begin=s'     => \$begin,
  'end=s'       => \$end,
  'verbose'     => \$verbose,
  'help'        => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fastq fastq file required\n")          unless $infile;
die("--begin or --end int value required\n")  unless(defined($begin) or defined($end));
die("--outfile outfile required\n")           unless $outfile;

## MAIN
# make pipe
## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-cutFastqByGz-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);
my $pipe = "$tmpdir/reads.pipe";

system("mkfifo $pipe");
# gunzip to pipes
print STDERR "[DEBUG] Reading $infile\n";
system("gunzip -c ".$infile." > $pipe &");

open(OUT, ">".$outfile) or die "Can't open file ".$!."\n";

my $in = new Iterator::FastqDb($pipe) or die("Unable to open Fastq file, $pipe\n");

if(!$begin){
  $begin = 0;
}
if(!$end){
  $end = 0;
}

while( my $curr = $in->next_seq() ){
   my $length = length($curr->seq()) - $end - $begin;

   #print STDERR "[DEBUG] header: ".$curr->header."\n";
   my $passed = "true"; 
   if($begin >= length($curr->seq())){
      print STDERR "[DEBUG] header: ".$curr->header."\n";
      print STDERR "[DEBUG] length curr seq: ".length($curr->seq)."\n";
      print STDERR "[DEBUG] begin: ".$begin."\n";
      print STDERR "[DEBUG] end: ".$end."\n";
      print STDERR "[DEBUG] length: ".$length."\n";
      print STDERR "[DEBUG] First nucleotides to cut must be lower than read length.\n";
      $passed = "false"; #next;
      die;
   }
   if($end >= length($curr->seq())){
      print STDERR "[DEBUG] header: ".$curr->header."\n";
      print STDERR "[DEBUG] length curr seq: ".length($curr->seq)."\n";
      print STDERR "[DEBUG] begin: ".$begin."\n";
      print STDERR "[DEBUG] end: ".$end."\n";
      print STDERR "[DEBUG] length: ".$length."\n";
      print STDERR "[DEBUG] Last nucleotides to cut must be lower than read length.\n";
      $passed = "false"; #next;
      die;
   }
   if($passed eq "true"){
      print OUT $curr->header."\n".substr($curr->seq, $begin, $length)."\n+\n".substr($curr->qual, $begin, $length)."\n";
   }else{
       #print OUT $curr->output."\n"; #else print unaletered fastq to not disrupt pairs.
       die;
   }
}
close(OUT);

exit;
