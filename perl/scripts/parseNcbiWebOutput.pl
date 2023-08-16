#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
parseNcbiWebOutput.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
--outdir <string> : outdir

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outdir);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'outdir=s'  => \$outdir,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   my @row = split(/\t/, $_);
   $hash{$row[0]."==".$row[1]} = $_;
}
close(IN);

for my $key (keys %hash){
   my @key = split(/==/, $key);
   my $query = $key[0];
   open(OUT, ">", $outdir."/".$query.".tsv") or die "Can't open ".$outdir."/".$query.".tsv";
   close(OUT);
}

for my $key (keys %hash){
   my @key = split(/==/, $key);
   my $query = $key[0];
   open(OUT, ">>", $outdir."/".$query.".tsv") or die "Can't open ".$outdir."/".$query.".tsv";
   print OUT $hash{$key}."\n";
   close(OUT);
}

