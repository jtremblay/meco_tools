#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;
use Cwd 'abs_path';
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
--indir <string>  : indir
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $indir);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'indir=s' 	=> \$indir,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

$indir = abs_path($indir);

## MAIN

# First take hash
my %hashR1;
my %hashR2;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   $hashR1{$row[0]} = $row[1]; 
   $hashR2{$row[0]} = $row[1]; 
}
close(IN);

#print STDERR Dumper(\%hash);


sub eachFile{
   my $filename = $_;
   my $fullpath = $File::Find::name;
   #remember that File::Find changes your CWD, 
   #so you can call open with just $_

   if (-e $filename) { 
    
      if(substr($filename, -9) eq ".fastq.gz"){
          #print STDOUT "Compressing ".$fullpath." into .gz archive...\n";
         for my $key (keys %hashR1){
            if($filename =~ m/($key)(.*)(_R1_)(.*)(\.fastq.gz)/){
               print STDOUT "ln -s $fullpath $hashR1{$key}".$2."_R1_".$4.$5."\n";
               #print STDOUT "ln -s $fullpath $hash{$key}".$2."_R2_".$4.$5."\n";

               #print STDERR "filename: $filename\n";
               #print STDERR "fullpath: $fullpath\n";
               #print STDERR $hash{$key}.$2."_R1_".$4.$5."\n";
               #print STDERR $hash{$key}.$2."_R2_".$4.$5."\n";

               delete($hashR1{$key});
               next; 
            }
         }
         for my $key (keys %hashR2){
            if($filename =~ m/($key)(.*)(_R2_)(.*)(\.fastq.gz)/){
               print STDOUT "ln -s $fullpath $hashR2{$key}".$2."_R2_".$4.$5."\n";
               #print STDOUT "ln -s $fullpath $hash{$key}".$2."_R2_".$4.$5."\n";

               #print STDERR "filename: $filename\n";
               #print STDERR "fullpath: $fullpath\n";
               #print STDERR $hash{$key}.$2."_R1_".$4.$5."\n";
               #print STDERR $hash{$key}.$2."_R2_".$4.$5."\n";

               delete($hashR2{$key});
               next; 
            }
         }
         #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
      }
   }
} 


## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

exit;
