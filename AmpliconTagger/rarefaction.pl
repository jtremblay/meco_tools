#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
rarefaction.pl

PURPOSE:
Feature tables rarefaction wrapper. Takes in input multiple Feature tables and
rarefy it to the lowest value of all the Feature tables. Implements the
single_rarefaction.py from Qiime.

INPUT:
--infile <string>          : Feature table in Qiime format.
                             Can be multiple infile.
--threshold <float>        : Fraction of average at which a sample
                             is considered failed. 
--n <int>                  : If table is to be rarefied at <int> reads. 

OUTPUT:
--outfile <string>         : Filtered otu table in Qiime format
                             Must be the same number if outfile
                             as infiles.

NOTES:
Please make sure that the two header lines of the Feature table is 
properly formatted. For instance:

#Full Feature Counts
#FEATURE_ID 01  03  04  05  06  07  08  09  10  11  Consensus lineage	

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjulien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @infile, @outfile, $begin, $end, $threshold, $n);
my $verbose = 0;

GetOptions(
  'infile=s'    => \@infile,
  'outfile=s'   => \@outfile,
  'threshold=f' => \$threshold,
  'begin=s'     => \$begin,
  'end=s'       => \$end,
  'n=i'         => \$n,
  'verbose'     => \$verbose,
  'help'        => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile file required\n") unless @infile;
die("--outfile outfile required\n") unless @outfile;
die("must be equal number of --infile and --outfile args...\n") if(@infile != @outfile);

$threshold = 0.05 unless($threshold);

## MAIN
# first pass through all infiles to find the lowest abundance value.

# Then if --n <int> arg is supplied, normalized to the sepcified value.
if($n){
	foreach my $infile (@infile){
		my $outfile = shift(@outfile);
		$outfile =~ s/\.tsv/_$n\.tsv/;
		# Run single_rarefactions.py 
		my $cmd = "single_rarefaction.py";
		$cmd .= " -i ".$infile;
		$cmd .= " -o ".$outfile;
		$cmd .= " -d ".$n;
      #print "$cmd\n" if $verbose;
      print "$cmd\n";
		system($cmd);
	}
}
exit;
