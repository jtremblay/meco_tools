#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile_otu_table <string>    : Sequence file
--infile_map <string>          : mapping file
				
OUTPUT:
STDOUT (outfile_otu_table)

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_otu_table, $infile_map);
my $verbose = 0;

GetOptions(
   'infile_map=s' 	    => \$infile_map,
   'infile_otu_table=s'  => \$infile_otu_table,
   'verbose'             => \$verbose,
   'help'                => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

open(IN, "<".$infile_map) or die "Can't open $infile_map\n";
while(<IN>){
   next if($_ =~ m/^#/);
   chomp;
   my @row = split(/\t/, $_);
   my $old_sample_name = $row[0];
   my $new_sample_name = $row[1];

   $hash{$old_sample_name} = $new_sample_name; 
}
close(IN);

open(IN, "<".$infile_otu_table) or die "Can't open $infile_otu_table\n";
while(<IN>){
   chomp;
   if($_ =~ m/^#OTU ID/){
      my @row = split(/\t/, $_);
      my $last = pop(@row);
      my $first = shift(@row);
      print STDOUT "#OTU_ID\t";

      foreach my $el(@row){
         if(exists $hash{$el}){
            print STDOUT $hash{$el}."\t";
            delete($hash{$el});
         }else{
            print STDERR "Did not find $el\n";
         }
      }
      print STDOUT "$last\n";
   }else{
      print $_."\n";
   }
}
close(IN);

print STDERR Dumper(\%hash);
