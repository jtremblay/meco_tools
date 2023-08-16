#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
validateMiRNATarget.pl

PURPOSE:

INPUT:
--infiles <string>         : parsed mirnatarget output files separated by a ,
--infile_reference <string>  : psrnatarget output file.
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $infile_reference);
my $verbose = 0;

GetOptions(
   'infiles=s'          => \$infiles,
   'infile_reference=s' => \$infile_reference,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @infiles = split(/,/, $infiles);
for my $infile (@infiles){
    my $name = $infile;
    $name =~ s{^.*/}{};     # remove the leading path
    $name =~ s{\.[^.]+$}{}; # remove the extension
    my %hash;

    open(IN, "<".$infile_reference) or die "Can't open $infile_reference\n";
    while(<IN>){
        chomp;
        if($_ =~ m/#/){
            next;
        }
        my @row = split(/\t/, $_);
        $hash{$row[0]."_".$row[1]."_".$row[5]}{ref}++;
        $hash{$row[0]."_".$row[1]."_".$row[5]}{start_pos} = $row[5];
    }
    close(IN);

    open(IN, "<".$infile) or die "Can't open $infile\n";
    while(<IN>){
        chomp;
        if($_ =~ m/#/){
            next;
        }
        my @row = split(/\t/, $_);
        my $value =  $row[0]."_".$row[1]."_".$row[7];
        if(exists $hash{$value}){
            $hash{$value}{test}++;
        }else{
            $hash{$value}{test} = 1;
            $hash{$value}{start_pos} = $row[7];
        }
    }
    close(IN);

    my $j = 0;
    print STDERR "------------------------------------------------------------------\n";
    foreach my $key (sort {$hash{$a}{start_pos} <=> $hash{$b}{start_pos}} keys %hash) {
    #foreach my $name (sort {lc $a cmp lc $b} keys %planets) {
        print STDERR $key."\t".$hash{$key}{start_pos}."\t";
        if(defined($hash{$key}{ref})){
            print STDERR $hash{$key}{ref}."\t";
        }else{
            print STDERR "0\t";
        }
        
        if(defined($hash{$key}{test})){
            print STDERR $hash{$key}{test}."\n";
        }else{
            print STDERR "0\n";
        }
        $j++;
    }
    print $infile."\t".$j."\n";
}

#print STDERR Dumper(\%hash);

