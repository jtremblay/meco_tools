#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateMAGFileForAnvio.pl

PURPOSE:

INPUT:
--indir <string> : Directory where are located MAGs .fa files. Non-recursive.
				
OUTPUT:
STDOUT           : MAG file compatible with anvio. each line: <Contig Id>\t<MAG/Bin Id>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
   'indir=s'  => \$indir,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my @files = glob "$indir/*.fa";
for (0..$#files){
    my $curr_file = $files[$_];
    #print STDERR $curr_file."\n";
    my $name = $curr_file;
    $name =~ s{^.*/}{};     # remove the leading path  
    $name =~ s{\.[^.]+$}{}; # remove the extension
    $name =~ s/\./_/g;
    print STDERR "name: ".$name."\n";
    #$files[$_] =~ s/\.fa$//;
    open(IN, "<".$curr_file) or die "Can't open $curr_file\n";
    while (my $line = <IN>){
        chomp $line;
        if($line =~ m/^>(\S+)/){
            print STDOUT $1."\t".$name."\n";
        }
    }
    close(IN);
}

