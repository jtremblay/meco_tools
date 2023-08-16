#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
generateSubCPMMatrixFromContrasts.pl

PURPOSE:

INPUT:
--infile <string>    : CPM matrix (or other abundance matrix).
--indir <string>     : Directory containing edgeR contrast tables.
				
OUTPUT:
--outdir <string>    : Directory where CPM matrices of each contrast file 
                       will be written to disk.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $indir, $outdir);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'indir=s'  => \$indir,
   'outdir=s' => \$outdir,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }


$indir = abs_path($indir);

my %hash;
my %hash_filenames;

sub eachFile{
    my $filename = $_;
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_

    if (-e $filename) { 
    
        if(substr($filename, -10) eq "_edger.tsv"){
            print STDERR "Found ".$fullpath." file...\n";
            $hash_filenames{$filename} = $filename;
            open(IN, "<".$filename);
            while(<IN>){
                chomp;
                my @row = split(/\t/, $_);
                $hash{$row[0]}{$filename} = "";
            }
            close(IN);
        }
    }
}


## MAIN
find (\&eachFile, $indir);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = $row[0];

    if($. == 1){
        foreach my $key (keys %hash_filenames){
            open(OUT, ">".$hash_filenames{$key});
            print OUT $_."\n";
            close(OUT);
        }
    }else{ # if exists in hash, loop through geneid to find corresponding file for each geneid. open, print and close.
        if(exists $hash{$gene_id}){
            foreach my $filename (keys %{ $hash{$gene_id} }) {
                open(OUT, ">>".$filename);
                print OUT $_."\n";
                close(OUT);
            }
        }
   }
}
close(IN);
