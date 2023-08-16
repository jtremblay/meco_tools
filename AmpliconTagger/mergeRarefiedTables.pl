#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;
use Data::Dumper;
use File::Find;
use Cwd 'abs_path';
use List::Util qw(sum);

my $usage=<<'ENDHERE';
NAME:
mergeRarefiedTables.pl

PURPOSE:

INPUT:
--indir <string>    : Directory contianing Feature table (biom).
--infile <string>   : Feature table containing all FEATURE_IDs with their taxonomic lineage.

OUTPUT:
STDOUT              : feature table computed with the mean of each Feature

NOTES:
Takes all tsv Feature tables in a dir and will output one
single Feature table having the median value of each case.
Assumes all Feature tables in the dir have been rarefied to
the same depth. In the indir, all tsv have to be Feature tables.

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $indir, $infile);
my $verbose = 0;

GetOptions(
  'indir=s'          => \$indir,
  'infile=s'         => \$infile,
  'verbose'          => \$verbose,
  'help'             => \$help
);
if ($help) { print $usage; exit; }

sub median{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2){ #odd?
        return $vals[int($len/2)];
    }else{#even
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub mean {
    return sum(@_)/@_;
}

## VALIDATE

die("--infile file required\n") unless $infile;
die("--indir file required\n") unless $indir;
$indir = abs_path($indir);
#open(OUT_MEDIAN, ">".$out_median) or die "Can't open $out_median\n";

## MAIN
my @array;

sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) { 
    if($filename =~ m/.*_\d+.tsv/ && substr($filename, -4) eq ".tsv"){
      my $tsv_filename = $filename;
      push(@array, $fullpath)
    }
  }
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

my %hash;
my $number_of_tables = 0;
foreach my $table (@array){
    print STDERR "Processing $table\n";
    open(IN, "<".$table) or die "Cant open $table\n";
    my %hash_samples;
    while(<IN>){
        chomp;

        # There is a bug in rtk 0.93.2: when specifying --ns arg, there are two \t\t char introduced
        # in output file. Replace these \t\t by \t
        $_ =~ s/\t\t/\t/g;

        if($_ =~ m/^Rarefied/){
            my @row = split(/\t/, $_);
            #pop(@row); #remove taxonomy column
            my $i = 0;
            shift(@row);
            foreach my $el (@row){
                $hash_samples{$i} = $el;
                $i++;   
            }
            print STDERR "Dumper(hash_samples:)\n";
            print STDERR Dumper(\%hash_samples);
            #exit;
        }elsif($_ =~ m/^#/){
            next;
        }else{
            my @row = split(/\t/, $_);
            my $feature_id = $row[0];
            my $i = 0;
            shift(@row);
            #pop(@row);
            my $last_index = scalar @row;
            #exit;
            foreach my $el (@row){
                my $sample_id = $hash_samples{$i};
                #push @{$hash{$feature_id}{$sample_id}}, $el;
                if($i != $last_index){
                    $hash{$feature_id}{$sample_id} .= int($el).":";
                }else{
                    $hash{$feature_id}{$sample_id} = $el; #el of last index is supposed to be taxonomic lineage.
                } 
                $i++;
            }
        }
    }
    close(IN);
    $number_of_tables++;
}

#print STDERR Dumper(\%hash);
my $i = 0;
print STDOUT "#FEATURE_ID";
#print OUT_MEDIAN "#FEATURE_ID";
foreach my $feature_id (sort keys %hash){
    foreach my $sample_id (sort keys %{ $hash{$feature_id} }) {
		print STDOUT "\t".$sample_id;
        #print OUT_MEDIAN "\t".$sample_id;
	}
	if($i == 0){
		last;
	}
}
print STDOUT "\ttaxonomy\n";

# Finally put taxonmic values
open(IN, "<".$infile) or die "Cant open $infile\n";
while(<IN>){
    chomp;
    if($_ =~ m/^#FEATURE_ID/){ next;}
    my @row = split(/\t/, $_);
    my $feature_id = shift(@row);
    my $taxonomy = pop(@row);
    if(exists $hash{$feature_id}){
        $hash{$feature_id}{taxonomy} = $taxonomy;
    }
}
close(IN);

foreach my $feature_id (sort keys %hash){
    my $taxonomy = "";
	print STDOUT $feature_id;
    foreach my $sample_id (sort keys %{ $hash{$feature_id} }) {
        if($sample_id eq "taxonomy"){
            $taxonomy = $hash{$feature_id}{$sample_id};
        }else{
            my @values = split(/:/, $hash{$feature_id}{$sample_id});
            #print OUT_MEDIAN "\t".median(@values);
            #print STDOUT "\t".mean(@values);
            print STDOUT "\t".sum(@values)/$number_of_tables;
        }
    }
	print STDOUT "\t".$taxonomy."\n";
}

exit;
