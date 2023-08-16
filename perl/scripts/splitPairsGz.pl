#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
splitPairs.pl

PURPOSE:
Parallel wrapper for splitPairsGz.pl

INPUT:
--infile <string>    :  Sequence file
--outfile_1 <string> :  Sequence file
--outfile_2 <string> :  Sequence file

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $outfile_1, $outfile_2);

my $verbose = 0;

GetOptions(
  'infile=s'      => \$infile,
  'outfile_1=s'   => \$outfile_1,
  'outfile_2=s'   => \$outfile_2,
  'verbose'       => \$verbose,
  'help'          => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--outfile_1 arg required\n") unless($outfile_1);
die("--outfile_2 arg required\n") unless($outfile_2);

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-SplitPairsGz-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);


#################
## SUBROUTINES ##
#################
# Declare output;
my $z1 = new IO::Compress::Gzip $outfile_1 or die "gzip failed: $GzipError\n";
my $z2 = new IO::Compress::Gzip $outfile_2 or die "gzip failed: $GzipError\n";

my $pipe = "$tmpdir/reads.pipe";
system("mkfifo $pipe");
system("gunzip -c ".$infile." > $pipe &");

my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fastq file, $pipe\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
   my $pair = $curr->pair;
   
   if($pair == 1){
      $z1->print($curr->output);
   }elsif($pair == 2){
      $z2->print($curr->output);
   }else{
      die "Contains no pair information...\n";
   }
}
$z1->close();
$z2->close();

## REMOVE TEMP FILES
sub END{
  system("rm ".$tmpdir." -rf");
}
