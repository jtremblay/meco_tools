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
createInterleavedFastq.pl

PURPOSE:

INPUT:
--reads1 <string> : Sequence file
--reads2 <string> : Sequence file
				
OUTPUT:
STDOUT            : Paired interleaved fastq file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $reads1, $reads2, $outfile);
my $verbose = 0;

GetOptions(
   'reads1=s'  => \$reads1,
   'reads2=s'  => \$reads2,
   'outfile=s' => \$outfile,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $tmpdir = File::Temp->newdir(
    "tmpdir-mergePairs-XXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

my $pipe1 = "$tmpdir/reads1.pipe";
my $pipe2 = "$tmpdir/reads2.pipe";
system("mkfifo $pipe1");
system("mkfifo $pipe2");
system("gunzip -c ".$reads1." > $pipe1 &");
system("gunzip -c ".$reads2." > $pipe2 &");

my $ref_fastq_db1 = Iterator::FastqDb->new($pipe1) or die("Unable to open Fastq file, $pipe1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($pipe2) or die("Unable to open Fastq file, $pipe2\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
   my $curr2 = $ref_fastq_db2->next_seq();
   
   my $base1 = $curr1->base;
   my $base2 = $curr2->base;

   if($base1 eq $base2){
      print STDOUT $curr1->output;
      print STDOUT $curr2->output;
   }else{
      print STDERR "Base 1 does not match Base 2 ".$base1." != ".$base2."\n";
      die();
   }
}
system("rm -f $pipe1");
system("rm -f $pipe2");
exit;
