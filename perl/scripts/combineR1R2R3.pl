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
combineR1R2R3.pl

PURPOSE:
Integrate barcode sequence into fastq file. All files shoulb be in .fastq.gz format.

INPUT:
--infile_R1 <string>    :  R1 Sequence file
--infile_R2 <string>    :  R2 Sequence file
--infile_R3 <string>    :  Index Sequence file
--outfile_R1 <string>   :  Sequence file
--outfile_R2 <string>   :  Sequence file

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2, $infile_R3, $outfile_R1, $outfile_R2);

my $verbose = 0;

GetOptions(
  'infile_R1=s'   => \$infile_R1,
  'infile_R2=s'   => \$infile_R2,
  'infile_R3=s'   => \$infile_R3,
  'outfile_R1=s'  => \$outfile_R1,
  'outfile_R2=s'  => \$outfile_R2,
  'verbose'       => \$verbose,
  'help'          => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile_R1 arg required\n") unless($infile_R1);
die("--infile_R2 arg required\n") unless($infile_R2);
die("--infile_R3 arg required\n") unless($infile_R3);
die("--outfile_R1 arg required\n") unless($outfile_R1);
die("--outfile_R2 arg required\n") unless($outfile_R2);

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
my $z1 = new IO::Compress::Gzip $outfile_R1 or die "gzip failed: $GzipError\n";
my $z2 = new IO::Compress::Gzip $outfile_R2 or die "gzip failed: $GzipError\n";

my $pipe1 = "$tmpdir/reads1.pipe";
my $pipe2 = "$tmpdir/reads2.pipe";
my $pipe3 = "$tmpdir/reads3.pipe";
system("mkfifo $pipe1");
system("mkfifo $pipe2");
system("mkfifo $pipe3");
system("gunzip -c ".$infile_R1." > $pipe1 &");
system("gunzip -c ".$infile_R2." > $pipe2 &");
system("gunzip -c ".$infile_R3." > $pipe3 &");

my $fastq_db1 = Iterator::FastqDb->new($pipe1) or die("Unable to open Fastq file, $pipe1\n");
my $fastq_db2 = Iterator::FastqDb->new($pipe2) or die("Unable to open Fastq file, $pipe2\n");
my $fastq_db3 = Iterator::FastqDb->new($pipe3) or die("Unable to open Fastq file, $pipe3\n");
while( my $curr1 = $fastq_db1->next_seq() ) {
    my $curr2 = $fastq_db2->next_seq();
    my $curr3 = $fastq_db3->next_seq();

    my $R3 = $curr3->seq();

    my $R1 = $curr1->base."#".$R3."/1\n".$curr1->seq."\n+\n".$curr1->qual()."\n";
    my $R2 = $curr2->base."#".$R3."/2\n".$curr2->seq."\n+\n".$curr2->qual()."\n";

    $z1->print($R1);
    $z2->print($R2);

    if(!($curr1->base() eq $curr2->base() && $curr3->base() eq $curr2->base())){
       die "Basenames do not match...\n";
    }
}
$z1->close();
$z2->close();

## REMOVE TEMP FILES
sub END{
  system("rm ".$tmpdir." -rf");
}
