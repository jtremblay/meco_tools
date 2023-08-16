#!/usr/bin/env perl

use strict;
use warnings;

use POSIX;
use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
rtkSaturation.pl

PURPOSE:
Feature tables rarefaction wrapper for the RTK package. 

INPUT:
--infile <string>          : OTU table in tsv. (first column if not labeled (i.e. contains
                             only a tab.
--remove_last_col <string> : 'true' or 'false' if set to true, will remove last column before computing max, min values.
                             For instance, if using a feature table with a taxonomy last column,
                             this flag should be set. Default = true;

OUTPUT:
STDOUT                     : comma separated values.

NOTES:
Please make sure that the two header lines of the OTU table is 
properly formatted. For instance:

#contig_id sample1  sample2  sample3  sample4

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $number_of_points, $remove_last_col);
my $verbose = 0;

GetOptions(
  'infile=s'           => \$infile,
  'number_of_points=i' => \$number_of_points,
  'remove_last_col=s'  => \$remove_last_col,
  'verbose'            => \$verbose,
  'help'               => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile file required\n") unless $infile;
$number_of_points = 40 unless($number_of_points);
$remove_last_col = "true" unless($remove_last_col);
## MAIN

# first pass through infile (all columns/samples) to find the lowest and highest abundance value.
my $max = 0;
my $min = 999999999;
my @final_sums;
open(IN, "<".$infile) or die "Can't open infile ".$infile."\n";
my @sums;
my $counter = 0;  

while(<IN>){
  chomp;
  if($_ =~ m/#/ || $. == 1){
    next;
  }
  
  my @row = split(/\t/, $_);
  shift(@row);
  if($remove_last_col eq "true"){
    pop(@row); 
  }

  # For first row
  if($. == 2){ # actually, the first row is the header.
    for(my $i=0; $i<@row; $i++){
      $sums[$i] = floor($row[$i]);
    }  
      
    #$counter++;
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

  #$counter++;
}
#Push sums of these samples into the final sums array.
foreach my $sum (@sums){
#  print STDERR $sum."\t";
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

print STDERR "Will use Max value: ".$max."\n";
print STDERR "Will use Min value: ".$min."\n";
my $x = 0;
#my $number_of_points = 20;
#my $number_of_points = floor($max / 50);
my $step = floor($max / $number_of_points);
my @points_array = (100,200,400,800,1000);
my $start = 0;
#push(@points_array, $start);
for(my $j=1;$j<=$number_of_points;$j++){
    my $curr_value = ($j * $step) + $start;
   push(@points_array, $curr_value); 
}
$x = floor($number_of_points * $step) - 10;
$points_array[@points_array - 1] = $x;

print STDOUT join(",", @points_array);

#print STDERR "Performing rarefaction at baseline <".$x."> for each tables.\n"; 
# Then, normalize all tables upon that value.
# Run rtk
#rtk memory -i ./test.tsv -o ./output/ -r 10 -d 1000 
#my $cmd = "rtk ".$mode;
#$cmd .= " -i ".$infile;
#$cmd .= " -o ".$outdir."/";
#$cmd .= " -r ".$perm;
#$cmd .= " -d ".join(",", @points_array);
#print STDERR $cmd."\n";
#system($cmd);
#$? != 0 ? die "command failed: $!\n" : print STDERR $cmd." done!\n"; 

