#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
keepUblastBestHits.pl

PURPOSE:
Takes keeps best n hits from ublast table.

INPUT:
--infiles <string> : Sequence file
--n <int>          : Number of best hits to keep (default = 1)	
--e <float>        : 1e-02 default
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $n, $e);
my $verbose = 0;

GetOptions(
  'infiles=s' => \$infiles,
  'n=i'       => \$n,
  'e=f'       => \$e,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

die("--infiles missing...\n") unless($infiles);
$n = 1 unless($n);
$e = 1e-02 unless($e);
## MAIN

my @infiles = split(/,/, $infiles);
my %hash;
foreach my $infile(@infiles){
   open(IN, "<".$infile) or die "Can't open $infile\n";
   print STDERR "Processing $infile\n";
   while(<IN>){
     chomp;
     #if($_ !~ m/^\d+;#/){
       #$header .= $_."\n";
       #next;
     #}
   
     my @row = split(/\t/, $_);
     #my $id = (split /\t/, $_)[0];
     my $id = $row[0];
     my $target = $row[1];
     my $evalue = $row[10];
  
     # if exists, check if evalue is lower. If so, replace. 
     if(exists $hash{$id}){
       if($evalue < $hash{$id}{EVALUE}){
          $hash{$id}{EVALUE} = $evalue; 
          $hash{$id}{LINE} = $_; 
       }
   
     # if does not exists, populate hash with query-line-evalue
     }else{
       $hash{$id}{EVALUE} = $evalue;
       $hash{$id}{LINE} = $_;
     }  
   }
   close(IN);
}

#print STDOUT $header;

for my $k1 (sort keys %hash) {
    if($hash{$k1}{EVALUE} < $e){
        print STDOUT $hash{$k1}{LINE}."\n";
    }
  #for my $k2 (sort keys %{ $hash{$k1} }){
  #  print STDOUT $hash{$k1}{$k2}."\n";

  #}
}
exit;
