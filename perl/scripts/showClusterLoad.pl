#!/usr/bin/env perl

use strict;
#use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
showClusterLoad.pl

PURPOSE:
This script prints a summary of cluster load on 
guillimin compute cluster.

INPUT:
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $output = `qstat -a`;

my @output = split(/\n/, $output);

my %hash;
my %hash_q;

foreach my $line ( @output ) {

   next if($line !~ m/^\d+/);

   my @row = split(/\s+/, $line);

   if($row[9] eq "R"){
      my $value_5 = $row[5];
      my $value_6 = $row[6];
      $value_5 = 0 unless($value_5 =~ m/\d+/);
      $value_6 = 0 unless($value_6 =~ m/\d+/);
      $hash{$row[2]}{NODES_R} += $value_5;
      $hash{$row[2]}{CORES_R} += $value_6; 
   }elsif($row[9] eq "Q"){
      my $value_5 = $row[5];
      my $value_6 = $row[6];
      $value_5 = 0 unless($value_5 =~ m/\d+/);
      $value_6 = 0 unless($value_6 =~ m/\d+/);
      $hash{$row[2]}{NODES_Q} += $value_5;
      $hash{$row[2]}{CORES_Q} += $value_6; 
   }
}


printf("%-15s %-15s %-15s %-15s %-15s", "#QUEUE", "NODES-R", "CORES-R", "NODES-Q", "CORES-Q");
print "\n";
for my $k1 (keys %hash){
   printf("%-15s %-15s %-15s %-15s %-15s", $k1, $hash{$k1}{NODES_R}, $hash{$k1}{CORES_R}, $hash{$k1}{NODES_Q}, $hash{$k1}{CORES_Q});
   print "\n";
}

exit;

