#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
getTophatAlignRate.pl

PURPOSE:

INPUT:
--indir <string> : Sequence file
				
OUTPUT:

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
    'help' 		  => \$help
);
if ($help) { print $usage; exit; }

my %hash;

sub commify {
  my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text
}

sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if($filename =~ m/align_summary.txt/){
			#print STDOUT "Found ".$fullpath."...\n";

      my @path = split(/\//, $fullpath);
      my $currFile = pop(@path);
      my $run = pop(@path);
      my $sampleName = pop(@path);

      open(IN, "<".$fullpath) or die "Can't open $fullpath\n";
      while(<IN>){
        chomp;
        #print STDERR $_."\n";
        if($_ =~ m/(.*\%) concordant pair alignment rate/){
          #40.6% concordant pair alignment rate
          $hash{$sampleName}{$run}{perc} = $1;
        }if($_ =~ m/Input:\s+(\d+)/){
          $hash{$sampleName}{$run}{input} = commify($1);
          #die "Did not found % alignment value...\n";
        }if($_ =~ m/Aligned pairs:\s+(\d+)/){
          $hash{$sampleName}{$run}{ap} = commify($1);
          #die "Did not found % alignment value...\n";
        }

      }

      #print STDERR $fullpath."\n";
      #print STDERR $currFile."\n";
      #print STDERR $run."\n";
      #print STDERR $sampleName."\n";
		}
	}
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

print STDOUT "#sampleName\trunId\tleftReads\trightReads\talignedPairs\talignedPerc\n";
for my $k1 (sort keys %hash) {
    for my $k2 (sort keys %{ $hash{$k1} }){
      print STDOUT "$k1\t$k2\t$hash{$k1}{$k2}{input}\t$hash{$k1}{$k2}{input}\t$hash{$k1}{$k2}{ap}\t$hash{$k1}{$k2}{perc}\n";
    }
}

exit;

