#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
keepHmmerBestHit.pl

PURPOSE:
Takes keeps best n hits from blast table.

INPUT:
--infile <string> : .domtblout file
--n <int>         : Number of hits to keep (default:1).
--e <int>         : Evalue cutoff. Default = 1e-05
--q <int>         : minimum query length. Default = 100
--a <int>         : minimum alignment length. Default = 100

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $n, $e, $q, $a);
my $verbose = 0;

GetOptions(
  'infile=s'  => \$infile,
  'n=i'       => \$n,
  'e=f'       => \$e,
  'q=i'       => \$q,
  'a=i'       => \$a,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $header = "";
my %hash;
#my $i=1;
#my $j=0;

$a = 100 unless($a);
$q = 100 unless($q);
$e = 1e-05 unless($e);
$n = 1 unless($n);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
  chomp;
  if($_ =~ m/^#/){
    $header .= $_."\n";
    next;
  }

  my @row = split(/\s+/, $_);

  #print STDERR join("\n", @row);
  #exit;

  my $id = $row[3];
  my $evalue = $row[6];
  my $qlen = $row[5];
  my $start = $row[17];
  my $end = $row[18];
  my $alen = abs($end - $start);

  next if($evalue >= $e);
  next if($alen < $a);
  next if($qlen < $q);

  if(exists $hash{$id}){
    next;

  }else{
    $hash{$id} = $_;

  }  
}

print STDOUT $header;

for my $k1 (sort keys %hash) { 
  print STDOUT $hash{$k1}."\n";
  #for my $k2 (sort keys %{ $hash{$k1} }){
  #  print STDOUT $hash{$k1}{$k2}."\n";

  #}
}
exit;
