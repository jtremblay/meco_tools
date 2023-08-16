#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::File;
use Env qw/TMPDIR/;
use Parallel::ForkManager;
use File::Temp;
use Cache::FastMmap;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
generateAlphadivJobs.pl

PURPOSE:
Wrapper for single_rarefaction.pl because the parallel wrapper from Qiime
often hangs.

INPUT:
--infile <string>    : otu table in biom format.
--num_threads <int>  : Number of threads. Default=1.

OUTPUT:
--outdir <outdir>    :  outdir having all qscores in separate files

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $outdir, $num_threads, $n, $step, $permutations);
my $verbose = 0;
my $debug = 0;

## EXTERNAL TOOLS

GetOptions(
  'outdir=s'       => \$outdir,
  'infile=s'       => \$infile,
  'num_threads=i'  => \$num_threads,
  'n=i'            => \$n,
  'step=i'         => \$step,
  'permutations=i' => \$permutations,
  'verbose'        => \$verbose,
  'help'           => \$help
);
if ($help) { print $usage; exit; }

#VALIDATION
$num_threads = 1 unless($num_threads);

#my $basename = $infile;
#$basename =~ s{.*/}{};      # removes path  
#$basename =~ s{\.[^.]+$}{}; # removes extension

my @data;

# Generate biom file list
for(my $i=1;$i<1000;$i++){
    if($i %% $step){
        for(my $j=1;$j<=$permutations;$j++){
             push(@data, "rarefaction_".$i."_".$j.".biom");
        } 
    }
}
print STDERR @data;
exit;

##======== Parallel::ForkManager starts here ========##

my $pm = new Parallel::ForkManager($num_threads);

foreach my $data (@data) {
    # Forks and returns the pid for the child:
    my $pid = $pm->start and next; 

    #... do some work with $data in the child process ...
    #single_rarefaction.py -i ./gene_annotation/taxonomy_contigs/otu_table_normalized.biom -o ./alphadiv/all/rarefaction//RARIF_lqd_/rarefaction_1_0.biom   -d 1 ;
    #mv ./alphadiv/all/rarefaction//RARIF_lqd_/rarefaction_1_0.biom ./alphadiv/all/rarefaction///rarefaction_1_0.biom
    
    #single_rarefaction.py -i ./gene_annotation/taxonomy_contigs/otu_table_normalized.biom -o ./alphadiv/all/rarefaction//RARIF_lqd_/rarefaction_1_1.biom   -d 1 ; 
    #mv ./alphadiv/all/rarefaction//RARIF_lqd_/rarefaction_1_1.biom ./alphadiv/all/rarefaction///rarefaction_1_1.biom

    $pm->finish; # Terminates the child process
}
