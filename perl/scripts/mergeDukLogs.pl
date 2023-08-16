#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
mergeDukLogs.pl

PURPOSE:

INPUT:
--logs <string> : list of logs separated by comma
--ids <string>  : list of logs ids sepaeated by comma
				
OUTPUT:
STDOUT          : tsv of logs

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $logs, $ids);
my $verbose = 0;

GetOptions(
   'logs=s'  => \$logs,
   'ids=s'   => \$ids,
   'verbose' => \$verbose,
   'help'    => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @infiles1 = split(/,/, $logs); 
my @infiles2 = split(/,/, $ids);

my %hash1;
foreach my $file (@infiles1){
   my $prefix = shift(@infiles2);

   open(IN, "<".$file) or die "Can't open $file\n";
   while(<IN>){
      chomp;
      #print STDERR $_."\n";
      #my($matched) = $_ =~ m/Total number of matched reads:\s+(\d+)/g;
      # duk Total number of reads:    10778900
      # duk Total number of matched reads: 481670
      # bbduk:
      # #Total\t
      # #Matched\t
      #if($_ =~ m/Total number of matched reads:\s+(\d+)/){
      if($_ =~ m/^#Total\t(\d+)/){
          #print STDERR $1."\n";
         $hash1{$prefix}{TOTAL} = $1;
      }
      #if($_ =~ m/Total number of reads:\s+(\d+)/){
      if($_ =~ m/^#Matched\t(\d+)\t(\S+)/){
          #print STDERR $1."\n";
         $hash1{$prefix}{MATCH} = $1;
         $hash1{$prefix}{RATIO} = $2;
      }
      #if($_ =~ m/Match ratio:\s+(\d+\.\d+)/){
      #    #print STDERR $1."\n";
      #   $hash1{$prefix}{RATIO} = $1;
      #}
   }
   close(IN);
}
#print STDERR Dumper(\%hash1);

print STDOUT "sampleName\tTOTAL\tMATCH\tRATIO\n";
for my $key (keys %hash1){
   print STDOUT $key."\t";
   print STDOUT $hash1{$key}{TOTAL}."\t";
   print STDOUT $hash1{$key}{MATCH}."\t";
   print STDOUT $hash1{$key}{RATIO}."\n";
}

exit;
