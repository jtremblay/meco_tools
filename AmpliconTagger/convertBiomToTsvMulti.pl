#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;
use Data::Dumper;
use File::Find;
use Cwd 'abs_path';
use Parallel::ForkManager;

my $usage=<<'ENDHERE';
NAME:
convertBiomToTsvMulti.pl

PURPOSE:

INPUT:
--indir <string>    : Directory contianing OTU table (biom).
--num_threads <int> : Number of threads. default=1. 

OUTPUT:
STDOUT

NOTES:
Bug in single_rarefaction.py. lineages are on multiple lines...

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

## VALIDATE
die("--indir file required\n") unless $indir;
$indir = abs_path($indir);
$num_threads = 1 unless($num_threads);
## MAIN

my @array;

sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) { 
    
    if(substr($filename, -5) eq ".biom"){
      my $biom = $fullpath;
      my $tsv_filename = $filename;
      $tsv_filename =~ s/.biom/.tsv/g;
      my $tsv = $indir."/".$tsv_filename;
      push(@array, "biom convert -i $biom -o $tsv --table-type 'OTU table' --to-tsv --output-metadata-id='taxonomy' --header-key='taxonomy'")
    }
  }
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

my $pm = new Parallel::ForkManager($num_threads);
foreach my $command (@array){
    my $pid = $pm->start and next;
    print STDERR $command."\n";
    #system($command);
    #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly converted biom table into tsv table...\n";
    #$pm->finish;
    
	my $attempts_left = 3;
    LOOP: {
        if (!eval { system($command); 1 }) {
           warn $@;
           if (--$attempts_left) {
              warn "Retrying...\n";
              redo;
           } else {
              warn "Aborting.\n";
              $pm->finish(1);
           }
        }
    }
    $pm->finish(0);
}

exit;
