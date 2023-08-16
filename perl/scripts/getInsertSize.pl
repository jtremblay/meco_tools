#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;
use Cwd 'abs_path';
use File::Temp;
use Env qw/TMPDIR/;
use File::Basename;
use File::stat;

my $usage=<<'ENDHERE';
NAME:
getInsertSize.pl

PURPOSE:

INPUT:
--indir <string> : dir where are located your bam files.
				
OUTPUT:
STDOUT           : Spreadsheet with insert size per lib

NOTES:
nrc/samtools/1.9 nrc/nrc_tools/1.0

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

## Validate
die "--indir arg missing\n" unless($indir);

$indir = abs_path($indir);

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz
my $tmpDir = File::Temp->newdir(
   "prepareReadsForGamngsXXXXXXX",
   DIR => $TMPDIR."/",
   CLEANUP => 1
);


## MAIN
sub eachFile{
   my $filename = $_;
   my $fullpath = $File::Find::name;
   #remember that File::Find changes your CWD, 
   #so you can call open with just $_

   if (-e $filename) { 
    
      if($filename =~ m/(\S+\.bam$)/){
         my $bam = $1;
         next if($bam =~ m/sort/g);
         my $prefix = $1;
         $prefix =~ s{.*/}{};      # removes path  
         $prefix =~ s{\.[^.]+$}{}; # removes extension
         my $bam_fullpath = $fullpath;
         print STDERR "Found ".$fullpath." bam file... will print output in $tmpDir/$prefix.txt\n";
         #system("bamtools stats -in $fullpath -insert > $tmpDir/$prefix");
         my $filesize = stat($fullpath)->size;
         #next if($filesize < 1000);
         system("samtools view $fullpath | head -n 1000000 | getinsertsize.py - > $tmpDir/$prefix.txt");
         $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly computed insert size stats for ".$fullpath."\n";


         #Then open $tmpDir/$prefix.txt and only keep median insert size.
         open(IN, "<".$tmpDir."/".$prefix.".txt");
         while(<IN>){
            chomp;
            #Read span: mean 185.026863278, STD=64.2423040882
            if($_ =~ m/Read span: mean (\d+\.\d+), STD=(\d+.\d+)/){
               my $average = $1;
               my $std = $2;
               my $min = $average - $std;
               $min = 1 if($min <= 0);
               my $max = $average + $std;
               if($max <= $min){
                  $max = $min + 2;
               }
              
               print STDOUT "$prefix bwa ../raw_reads/$prefix\_R1.fastq.gz ../raw_reads/$prefix\_R2.fastq.gz $average 0.35 FR\n";
               #print STDOUT $bam_fullpath."\n".$min." ".$max."\n";
               last;
            }
         }
         close(IN);
      }
   }
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

exit;


