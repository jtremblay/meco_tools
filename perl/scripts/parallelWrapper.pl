#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Parallel::ForkManager;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $num_threads);
my $verbose = 0;

GetOptions(
   'infile=s'      => \$infile,
   'num_threads=i' => \$num_threads,
   'verbose'       => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @commands;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my $command = $_;
   push(@commands, $command);
}
close(IN);

my $pm = Parallel::ForkManager->new($num_threads);
 
#DATA_LOOP:
foreach my $command (@commands) {
# Forks and returns the pid for the child:
    my $pid = $pm->start and next;

    print STDERR $command."\n";
    system($command);
    $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly executed $command\n";

    $pm->finish; # Terminates the child process
}
$pm->wait_all_children;
exit;
