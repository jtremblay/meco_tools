#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;
use Parallel::ForkManager;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
fastqsGzToFastas.pl

PURPOSE

INPUT:
--infiles <string>    : Sequence Gz Fastqs file list sep by ,
--num_threads <int>   : number of threads.

OUTPUT:
--outfiles <string>   : Fasta files list sep by ,

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infiles, $outfiles, $num_threads);

my $verbose = 0;

GetOptions(
  'infiles=s'      => \$infiles,
  'outfiles=s'     => \$outfiles,
  'num_threads=i'  => \$num_threads,
  'verbose'        => \$verbose,
  'help'           => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infiles arg required\n") unless($infiles);
die("--outfiles arg required\n") unless($outfiles);
$num_threads = 1 unless($num_threads);

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-FastqsGzToFastas-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);


#################
## SUBROUTINES ##
#################
# Declare output;

my @fastqsGz = split(/,/, $infiles);
my @fastas = split(/,/, $outfiles);
my @commands;

for my $fastqGz (@fastqsGz){
    my $fasta = shift(@fastas);

    my $command = "gunzip -c $fastqGz | awk 'NR%4 == 1 {print \">\" substr(\$0, 2)} NR%4 == 2 {print}' > $fasta";
    push @commands, $command;
    #print STDERR "[DEBUG] $command\n";
    #system($command);
    #die "command failed: $!\n" if($? != 0);


}

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
## REMOVE TEMP FILES
sub END{
  system("rm ".$tmpdir." -rf");
}
