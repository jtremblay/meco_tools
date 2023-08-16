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
use POSIX;

my $usage=<<'ENDHERE';
NAME:
convertBiomToTsvMulti.pl

PURPOSE:

INPUT:
--indir <string>     : Directory contianing Feature table (biom).
--outdir <string>    : Directory contianing Feature table (biom).
--num_threads <int>  : Number of threads. default=1. 

OUTPUT:
STDOUT

NOTES:
Bug in single_rarefaction.py. lineages are on multiple lines...

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $indir, $outdir, $num_threads, $m);
my $verbose = 0;

GetOptions(
  'indir=s'       => \$indir,
  'outdir=s'      => \$outdir,
  'm=i'           => \$m,
  'num_threads=i' => \$num_threads,
  'verbose'       => \$verbose,
  'help'          => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--indir file required\n") unless $indir;
die("--outdir file required\n") unless $outdir;
$indir = abs_path($indir);
$num_threads = 1 unless($num_threads);
$m = 1000 unless($m);

## MAIN

my @array;
my @array2;
my @array3;

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
      push(@array, "biom convert -i $biom -o $tsv --table-type 'Feature table' --to-tsv --output-metadata-id='taxonomy' --header-key='taxonomy'");
      #my $tsv2 = $tsv;
      ##$tsv2 =~ s{.*/}{};      # removes path  
      #$tsv2 =~ s{\.[^.]+$}{}; # removes extension
      #$tsv2 = $tsv2."_d1000.tsv";
      
      # actually estimate the final number on the otu table labeling.
      my $tsv2 = $tsv;
      $tsv2 =~ s{.*/}{};
      $tsv2 =~ s{\.[^.]+$}{};
      #alpha_rarefaction_1_9_d1000.txt
      if($tsv2 =~ m/rarefaction_(\d+)_(\d+).*/){
        my $new_number = ceil($1 / 1000);
        $tsv2 = "rarefaction_".$new_number."_".$2.".tsv";
        
      }else{
        die "Problem in file parsing ... :".$tsv2."\n";
      }
      
      push(@array2, "divideFeatureTable.R -i $tsv -o $outdir/$tsv2 -m $m");
      my $biom2 = $tsv2;
      $biom2 =~ s{.*/}{};      # removes path  
      $biom2 =~ s{\.[^.]+$}{}; # removes extension
      $biom2 = $outdir."/".$biom2.".biom";
      

      push(@array3, "biom convert -i $outdir/$tsv2 -o $biom2 --table-type 'Feature table' --to-hdf5  --process-obs-metadata taxonomy");
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

$pm = new Parallel::ForkManager($num_threads);
foreach my $command (@array2){
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

$pm = new Parallel::ForkManager($num_threads);
foreach my $command (@array3){
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
