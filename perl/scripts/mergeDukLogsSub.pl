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

   my $k = 0;
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
      if($_ =~ m/^Reads Used:\s+\t(\d+)\t/){
          #print STDERR $1."\n";
         $hash1{$prefix}{TOTAL} = $1;
      }

      #if($_ =~ m/Total number of reads:\s+(\d+)/){
      if($_ =~ m/^Read 1 data:\s+\t/){
          $k = 1;
      }elsif($_ =~ m/^Read 2 data:\s+\t/){
          $k = 2;
      }
      
      if($_ =~ m/mapped:\s+\t/){
          my @row = split(/\t/, $_);
          $hash1{$prefix}{"R".$k}{RATIO} = $row[1];
          $hash1{$prefix}{"R".$k}{MATCH} = $row[2];
      }
   }
   close(IN);
}
#print STDERR Dumper(\%hash1);

print STDOUT "sampleName\tTOTAL\tMATCH_R1\tRATIO_R1\tMATCH_R2\tRATIO_R2\n";
for my $key (keys %hash1){
    print STDOUT $key."\t";
    print STDOUT $hash1{$key}{TOTAL}."\t";
    
    print STDOUT $hash1{$key}{R1}{MATCH}."\t";
    print STDOUT $hash1{$key}{R1}{RATIO}."\t";
    print STDOUT $hash1{$key}{R2}{MATCH}."\t";
    print STDOUT $hash1{$key}{R2}{RATIO}."\n";
    
}

exit;
