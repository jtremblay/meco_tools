#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Env qw(TMPDIR);
use File::Temp;
use Parallel::ForkManager;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
prodigalParallel.pl

PURPOSE

INPUT:
--infiles <string>       : list of fasta files sep by a ,
--num_threads <int>      : number of threads.

OUTPUT:
--outfiles <string>      : list of gff files list sep by ,
                           Should be in the same order as the list
                           provided in --infiles <list>
--outfiles_faa <string>  : list of faa files sep by ,
                           Should be in the same order as the list
                           provided in --infiles <list>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infiles, $outfiles, $outfiles_faa, $num_threads);

my $verbose = 0;

GetOptions(
  'infiles=s'      => \$infiles,
  'outfiles=s'     => \$outfiles,
  'outfiles_faa=s' => \$outfiles_faa,
  'num_threads=i'  => \$num_threads,
  'verbose'        => \$verbose,
  'help'           => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infiles arg required\n") unless($infiles);
die("--outfiles arg required\n") unless($outfiles);
die("--outfiles_faa arg required\n") unless($outfiles_faa);
$num_threads = 1 unless($num_threads);

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-parallelProdigal-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

# Declare output;

my @fastas = split(/,/, $infiles);
my @faas = split(/,/, $outfiles_faa);
my @gffs = split(/,/, $outfiles);
my @commands;

for my $fasta (@fastas){
    my $gff = shift(@gffs);
    my $faa = shift(@faas);

    my $command = "prodigal -p single -f gff -i $fasta -o $gff -a $faa";
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
