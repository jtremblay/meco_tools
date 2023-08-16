#!/usr/bin/env perl

use strict;
use warnings;

use Env qw(SCRATCH TMPDIR PICARD_HOME);
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
bamToFastq.pl

PURPOSE:
To loop recursively through a directory and compress each .fastq
to .gz file. 

INPUT:
--infile <string>  : Directory where itags runs are to be backed-up.

OUTPUT:
--reads1 <string>
--reads2 <string>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $reads1, $reads2);
my $verbose = 0;

GetOptions(
    'infile=s'   => \$infile,
    'reads1=s'   => \$reads1,
    'reads2=s'   => \$reads2,
    'verbose'    => \$verbose,
    'help'       => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--infile arg missing\n" unless($infile);

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz

## MAIN
my $cmd = "";
$cmd .= " java -Djava.io.tmpdir=$SCRATCH -XX:ParallelGCThreads=3 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar";
$cmd .= " INPUT=".$infile;
$cmd .= " FASTQ=".$reads1;
#$cmd .= " SECOND_END_FASTQ=".$reads2;
print STDERR "Converting $infile bam to fastq\n";
print STDERR "============\n".$cmd."\n=============\n";
system($cmd);
$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly decompressed ".$infile." into .fastq archive...\n";

exit;
