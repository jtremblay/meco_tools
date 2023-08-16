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
removeShortFastqsPEGz.pl

PURPOSE:
Script to filter weird fastqs...

INPUT:
--infile_R1 <string>  : fastq.gz infile
--infile_R2 <string>  : fastq.gz infile
--length <int>        : Remove reads with length lt <int>

OUTPUT:
--outfile_R1 <string> : trimmed fastq.gz file.
--outfile_R2 <string> : trimmed fastq.gz file.

NOTES:

BUGS/LIMITATIONS:
TODO Parallelize this script.

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2, $outfile_R1, $outfile_R2, $length);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'infile_R1=s'    => \$infile_R1,
  'infile_R2=s'    => \$infile_R2,
  'outfile_R1=s'   => \$outfile_R1,
  'outfile_R2=s'   => \$outfile_R2,
  'length=s'       => \$length,
  'verbose'        => \$verbose,
  'help'           => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile_R2 fastq.gz file required\n") unless $infile_R1;
die("--infile_R2 fastq.gz file required\n") unless $infile_R1;
die("--outfile_R1 outfile R1 required\n")   unless $outfile_R1;
die("--outfile_R2 outfile R2 required\n")   unless $outfile_R2;
die "--length missing"                      unless $length;

## MAIN
# make pipe
## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-cutFastqByGz-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);
my $pipe1 = "$tmpdir/reads1.pipe";
my $pipe2 = "$tmpdir/reads2.pipe";

system("mkfifo $pipe1");
system("mkfifo $pipe2");
# gunzip to pipes
print STDERR "[DEBUG] Reading $infile_R1 and $infile_R2\n";
system("gunzip -c ".$infile_R1." > $pipe1 &");
system("gunzip -c ".$infile_R2." > $pipe2 &");

open(OUT1, ">".$outfile_R1) or die "Can't open file ".$!."\n";
open(OUT2, ">".$outfile_R2) or die "Can't open file ".$!."\n";

my $in1 = new Iterator::FastqDb($pipe1) or die("Unable to open Fastq file, $pipe1\n");
my $in2 = new Iterator::FastqDb($pipe2) or die("Unable to open Fastq file, $pipe2\n");


while( my $curr1 = $in1->next_seq() ){
   my $curr2 = $in2->next_seq();

   if(length($curr1->seq) >= $length && length($curr2->seq) >= $length){
        print OUT1 $curr1->output;
        print OUT2 $curr2->output;
   }
}
close(OUT1);
close(OUT2);

exit;
