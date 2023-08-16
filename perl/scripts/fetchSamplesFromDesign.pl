#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--design <string>       : design file
--projectSheet <string> : nanuq project sheet
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $design, $projectSheet);
my $verbose = 0;

GetOptions(
    'design=s' 	      => \$design,
    'projectSheet=s' 	=> \$projectSheet,
    'verbose' 	      => \$verbose,
    'help' 		        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(DESIGN, "<".$design) or die "Can't open $design\n";

my %hash;
while(<DESIGN>){
  chomp;
  my @row = split(/\t/, $_);
  my $sampleName = $row[0];  

  if(exists $hash{$sampleName}){
    die "$sampleName is present at least twice in your design file. Please correct your design file.\n";
  }else{
    $hash{$sampleName} = 1;
  }
}
close(DESIGN);

open(PRJ, "<".$projectSheet) or die "Can't open $projectSheet\n";
while(<PRJ>){
  chomp;
  my @row = split(/\,/, $_);
  my $sampleName = $row[0];
  $sampleName =~ s/"//g;
  
  if(exists $hash{$sampleName}){
    print STDOUT $_."\n";
    delete($hash{$sampleName}); 
  } 
}
close(PRJ);

print STDERR "These are the following sample names in your design file that weren't found in the project sample sheet. Cheers.\n";
for my $key (keys %hash){
  print STDERR "$key\n";
}


