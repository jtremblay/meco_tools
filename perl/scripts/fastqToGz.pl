#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
fastqToGz.pl

PURPOSE:
To loop recursively through a directory and compress each .fastq
to .gz file. 

INPUT:
--indir <string>    : Directory where itags runs are to be backed-up.
--num_threads <int> : Number of threads (for pigz).

OUTPUT:
No output options.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $num_threads);
my $verbose = 0;

GetOptions(
  'indir=s'       => \$indir,
  'num_threads=i' => \$num_threads,
  'verbose'       => \$verbose,
  'help'          => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--indir arg missing\n" unless($indir);
$num_threads = 1 unless($num_threads);

$indir = abs_path($indir);

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz


## MAIN
sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) { 
    
    if(substr($filename, -6) eq ".fastq"){
      print STDOUT "Compressing ".$fullpath." into .gz archive...\n";
      system("pigz -p ".$num_threads." ".$fullpath);
      $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
    }
  }
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

exit;
