#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Parallel::ForkManager;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
sparCCWrapper.pl

PURPOSE:

INPUT:
--indir <string>    : Dir where bootsraps are located.
--n <int>           : Number of resampling. default = 20
--n_boot <int>      : Number of bootstrap. default = 100
--num_threads <int> : number of threads. default = 1


OUTPUT:
--outdir <string>   : outdirectory.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $n, $n_boot, $num_threads, $outdir);
my $verbose = 0;

GetOptions(
   'indir=s' 	    => \$indir,
   'outdir=s' 	    => \$outdir,
   'n=i'           => \$n,
   'n_boot=i'      => \$n_boot,
   'num_threads=i' => \$num_threads,
   'verbose' 	    => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

die "--indir missing\n" unless($indir);
die "--outdir missing\n" unless($outdir);
die "--n missing\n" unless($n);
die "--n_boot missing\n" unless($n_boot);
die "--num_threads missing\n" unless($num_threads);

$indir = abs_path($indir);
$outdir = abs_path($outdir);

## MAIN
my @commands;
for(my $i=0; $i<$n_boot; $i++){
   my $command = "SparCC.py $indir/boot_$i.txt -i $n --cor_file=$outdir/boot_$i.txt >> ./bootstrap_sparcc.log";
   push @commands, $command;
}

# Max 32 processes.
my $pm = new Parallel::ForkManager($num_threads); 

foreach my $command (@commands) {
   $pm->start and next; # do the fork

   #### Write code here:
   print STDERR "executing: $command\n";
   system($command);
   ####


   $pm->finish; # do the exit in the child process
}
$pm->wait_all_children;
exit;
