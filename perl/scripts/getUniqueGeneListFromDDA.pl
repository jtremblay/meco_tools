#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
getUniqueGeneListFromDDA.pl

PURPOSE:

INPUT:
--indir <string> : indirectory
				
OUTPUT:
STDOUT           : unique gene list

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
   'indir=s' 	=> \$indir,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
$indir = abs_path($indir);

my %hash;

sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) { 
    
    if(substr($filename, -10) eq "_edger.tsv"){
      
      open(IN, "<".$fullpath) or die "Can't open $fullpath\n";
      while(<IN>){
         chomp;
         next if($. == 1);
         my @row = split(/\t/, $_);
         $hash{$row[0]} = $row[0];
      }
      close(IN);
    }
  }
}

find (\&eachFile, $indir);

foreach my $key (keys %hash){
   print STDOUT $key."\n";
}
