#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Parallel::ForkManager;

my $usage=<<'ENDHERE';
NAME:
gunzipParallel.pl

PURPOSE:

INPUT:
--infiles <string>  : List of .fastq.gz files
--outfiles <string> : List of .fastq files
--num_threads <int> : Number of threads

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $outfiles, $num_threads);
my $verbose = 0;

GetOptions(
   'infiles=s'     => \$infiles,
   'outfiles=s'    => \$outfiles,
   'num_threads=i' => \$num_threads,
   'verbose' 	    => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

$num_threads = 1 unless($num_threads);

## MAIN
my @infiles = split(/,/, $infiles);
my @outfiles = split(/,/, $outfiles);

my $pm = new Parallel::ForkManager($num_threads);

foreach my $infile (@infiles){
   my $outfile = shift(@outfiles);
   $pm->start and next;
   print STDERR "$infile\t$outfile\n";
   system("gunzip -c ".$infile." > ".$outfile);
   $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly uncompressed ".$infile." into ".$outfile."...\n";
   $pm->finish;
}
$pm->wait_all_children;
exit;
