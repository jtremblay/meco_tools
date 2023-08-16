#!/usr/bin/env perl

use strict;
use warnings;

use POSIX;
use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
multipleRarefactions.pl

PURPOSE:
Feature tables rarefaction wrapper. Takes in input multiple Feature tables and
rarefy it to the lowest value of all the Feature tables. Implements the
single_rarefaction.py from Qiime.

INPUT:
--infile <string>          : Feature table in Biom format.
                             Can be multiple infile.
--infileTsv <string>       : Feature table in Qiime (tsv) format.
                             Can be multiple infile.
--threshold <float>        : Threshold to be used. For instance if 0.0001 means
                             that the samples will be rarefied at 0.01% of the total 
                             number of reads in the Feature table. For example, if we 
                             have a total 10,197,130 reads and we specify --threshold 0.0001
                             the script will rarefy at 1,020 reads per sample.
--n <int>                  : If table is to be rarefied at <int> reads. Optional - if not
                             specified, will rarefied at max value.
--m <int>                  : start rarefaction at <int> reads. 
--ratio <float>            : Will multiply the --m <int> and --n <int> by a ratio < 1.0
                             implemented because input Feature tables having decimal can
                             be problematic to rarefy with Qiime.
--step <int>               : Incrementation values in performing rarefactions. 
--permutations <int>       : Number of permutations. Default = 10.
--num_threads <int>        : Num threads.

OUTPUT:
--outdir <string>          : Filtered feature table in Qiime format
                             Must be the same number if outfile
                             as infiles.

NOTES:
Please make sure that the two header lines of the Feature table is 
properly formatted. For instance:

#Full Feature Counts
#FEATURE_ID 01  03  04  05  06  07  08  09  10  11  Consensus lineage  

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $infileTsv, $outdir, $step, $begin, $end, $threshold, $n, $perm, $num_threads, $m, $ratio, $number_of_points, $custom_max);
my $verbose = 0;

GetOptions(
  'infile=s'           => \$infile,
  'infileTsv=s'        => \$infileTsv,
  'outdir=s'           => \$outdir,
  'threshold=f'        => \$threshold,
  'n=i'                => \$n,
  'm=i'                => \$m,
  'ratio=f'            => \$ratio,
  'permutations=i'     => \$perm,
  'num_threads=i'      => \$num_threads,
  'number_of_points=i' => \$number_of_points,
  'step=i'             => \$step,
  'custom_max=i'       => \$custom_max,
  'verbose'            => \$verbose,
  'help'               => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile file required\n") unless $infile;
die("--outdir outdir required\n") unless $outdir;
#die("--step step required\n") unless $step;
#die("--n n required\n") unless $n;
$ratio = 1 unless($ratio);
#die("--ratio <float> has to be lower than 1\n") if($ratio > 1);
$perm = 10 unless($perm);
$num_threads = 1 unless($num_threads);
$m = 1 unless($m);
$number_of_points = 20 unless($number_of_points);
my $multipleRarefactions = `which multiple_rarefactions.py`; die if(!defined($multipleRarefactions)); chomp $multipleRarefactions;

$threshold = 0.05 unless($threshold);

## MAIN

# first pass through all infiles to find the lowest and highest abundance value.
my $max = 0;
my $min = 999999999;
my @final_sums;
open(IN, "<".$infileTsv) or die "Can't open infile ".$infileTsv."\n";
my @sums;
my $counter = 0;  

while(<IN>){
  chomp;
  if($_ =~ m/#/){
    next;
  }
  
  my @row = split(/\t/, $_);
  pop(@row); shift(@row);

  # For first row
  if($counter == 0){
    for(my $i=0; $i<@row; $i++){
      $sums[$i] = floor($row[$i]);
    }  
      
    $counter++;
    next;
  }

  for(my $i=0; $i<@row; $i++){
    if(defined($sums[$i])){
      $sums[$i] = floor($sums[$i] + floor($row[$i]));
    }else{
      print STDERR "Not defined...at indice ".$i." file: ".$infile." line ".$.."\n";
      #die ?
      #$sums[$i] = 0;
      #$sums[$i] = $sums[$i] + $row[$i];    
    }
  }

  $counter++;
}
#Push sums of these samples into the final sums array.
foreach my $sum (@sums){
#  print $sum."\t";
  push(@final_sums, $sum);
}
#print "========================\n";
close(IN);

#For Debug...
#foreach(@final_sums){
#  print STDERR "[DEBUG] ".$_."\n";
#}

# pick max value.
@final_sums = sort @final_sums;
#my $stat = Statistics::Descriptive::Sparse->new();
#foreach my $value (@final_sums){
#  $stat->add_data($value);
#}
#my $avg = $stat->mean();

# Then only consider a minimal value to be good if it is above at least 10% of the average sample abundance.
# This is to reject failed samples.
foreach my $value (@final_sums){
    $min = $value if($value < $min);
    $max = $value if($value > $max);
}
if($custom_max){
    $max = $custom_max;
}
print STDERR "Custom max was specified: ".$max."\n";
#$max = $max - $counter;
print STDERR "Will use Max value: ".$max."\n";
print STDERR "Will use Min value: ".$min."\n";
my $x = 0;
if($n){
    $x = $n;
}else{
    #my $number_of_points = 20;
    #my $number_of_points = floor($max / 50);
    $step = floor($max / $number_of_points);
    $x = floor($number_of_points * $step) - 10;
}
#print STDERR "Performing rarefaction at baseline <".$x."> for each tables.\n"; 
# Then, normalize all tables upon that value.
# Run single_rarefactions.py 
my $cmd = "parallel_multiple_rarefactions.py";
$cmd .= " -i ".$infile;
$cmd .= " -m ".$m;
$cmd .= " -x ".$x;
$cmd .= " -s ".$step;
$cmd .= " -n ".$perm;
$cmd .= " -o ".$outdir."/";
$cmd .= " -O ".$num_threads;
print STDERR $cmd."\n";
system($cmd);
$? != 0 ? die "command failed: $!\n" : print STDERR $cmd." done!\n"; 

