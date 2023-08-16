#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
hpss_itags.pl

PURPOSE:
To loop recursively through a directory and compress each .fastq
to .gz file. 

INPUT:
--indir <string>  : Directory where itags runs are to be backed-up.

OUTPUT:
No output options.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
    'indir=s' 	=> \$indir,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--indir arg missing\n" unless($indir);

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz


my $numberOfFiles = 0;
my %hash;
my %hash_lineage;
my %hash_filename;

sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if(substr($filename, -4) eq ".tsv"){
          $hash_filename{$fullpath} = $fullpath;
          $numberOfFiles++;
		}
	}
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

print STDERR "number of files: ".$numberOfFiles."\n";

my $i=0;
for my $key (keys %hash_filename){
   open(IN, "<".$key) or die "Can't open $key\n";
   print STDERR "Processing file: ".$key."\n";
   
   while(<IN>){
      chomp; 
      if($_ =~ m/^#/){
         print STDOUT $_."\n" if($i == 0);
         next;
      }else{
         my @row = split(/\t/, $_);
         my $otuID = shift(@row);
         my $lineage = pop(@row);
         $hash_lineage{$otuID} = $lineage;

         my $j = 1;
         foreach my $value (@row){
            if(exists $hash{$otuID}{$j}){
               $hash{$otuID}{$j} += $value; 
            }else{
               $hash{$otuID}{$j} = $value; 
            }
            $j++;
         }
      }
   }
   close(IN);
   $i++;
}

foreach my $key (keys %hash) {
   print STDOUT $key;
   foreach my $key2 (sort {$a <=> $b} keys %{ $hash{$key} }) {
      my $currValue = $hash{$key}{$key2};
      my $normalizedValue = $currValue / $numberOfFiles;
      $normalizedValue = sprintf "%.0f", $normalizedValue;
      print STDOUT "\t".$normalizedValue;
   }
   print STDOUT "\t".$hash_lineage{$key}."\n";
}

exit;
