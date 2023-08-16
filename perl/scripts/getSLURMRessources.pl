#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
getSLURMRessources.pl

PURPOSE:

INPUT:
--infile <string> : GenPipes job submission list file.
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my @jobs;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   push(@jobs, $row[0]);

}
close(IN);
my $jobs = join(",", @jobs);

my $saccts = `sacct --format=jobid,jobname%50,maxvmsizenode,maxvmsize,avevmsize,avecpu,consumedenergy,MaxRSS,alloccpus,elapsed,exitcode -j $jobs`;
print STDOUT $saccts."\n";

