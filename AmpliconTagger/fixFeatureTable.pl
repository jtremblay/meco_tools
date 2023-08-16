#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
rarefaction.pl

PURPOSE:

INPUT:
--infile <string> : Feature table (tsv).


OUTPUT:
STDOUT            : Filtered feature table in Qiime format
                    Must be the same number if outfile
                    as infiles.

NOTES:
Bug in single_rarefaction.py. lineages are on multiple lines...

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjulien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
  'infile=s'    => \$infile,
  'verbose'     => \$verbose,
  'help'        => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile file required\n") unless $infile;

## MAIN
my %hash;
my $id = "";
my $lineage = "";
my $abundance = "";
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   if($_ =~ m/^#/){
      print $_."\n";
      next;
   }

   

   if($_ =~ m/^\s\'\S__/){
      #print STDERR "found $_\n";
      #Means we are in broken tax lineage.
      $lineage = $_;
      $lineage =~ s/\['//g;
      $lineage =~ s/'\]//g;
      $lineage =~ s/'/;/g;
      $lineage =~ s/^ '/;/g;
      #$lineage =~ s/ /;/g;
   
      #my @lineage = split(/ /, $lineage);
      
   }else{
      #print STDERR "not found\n";

      my @row = split(/\t/, $_);
      $id = shift(@row);
      $lineage = pop(@row);
      $abundance = join("\t", @row);
      #print STDERR "id: ".$id."\n";
      #print STDERR "lineage: $lineage\n";
 
      $lineage =~ s/\['//g;
      $lineage =~ s/'\]//g;
      $lineage =~ s/' '/;/g;
      $lineage =~ s/'//g;
      #$lineage =~ s/ /;/g;
   
      #my @lineage = split(/ /, $lineage);
   }
   
   if(exists $hash{$id}){
      $hash{$id}{lineage} .= $lineage;
   }else{
      $hash{$id}{abundance} = $abundance;
      $hash{$id}{lineage} = $lineage;
   }
   
}
#print STDERR Dumper(\%hash);

for my $key (keys %hash){
   print STDOUT $key."\t";
   print STDOUT $abundance = $hash{$key}{abundance}."\t";
   my $lineage = $hash{$key}{lineage};
   $lineage =~ s/;;/;/g;
   $lineage =~ s/; ;/;/g;
   $lineage =~ s/; /;/g;
   $lineage =~ s/ ;/;/g;
   print STDOUT $lineage."\n";
}


exit;
