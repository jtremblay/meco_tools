#!/usr/bin/env perl
#
use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use File::Find;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateReadsetSheet.pl

PURPOSE:
Takes an indir containing fastqs and generate a readset sheet.

INPUT:
--indir <string>   : sheet having only sample prefix in them.
--type <string>    : Default = PAIRED_END

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $type);
my $verbose = 0;

GetOptions(
    'indir=s'  => \$indir,
    'type=s'   => \$type,
    'verbose'  => \$verbose,
    'help'     => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "missing --indir\n" unless($indir);

## MAIN
my %hash1;
my %hash2;
my %hash;
my $header;
my $indir_relpath = $indir;
$indir = abs_path($indir);
my $prefix;
$type = "PAIRED_END" unless($type);

#system("mkdir -p ".$indir."/raw_reads");

sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_
  #print STDERR "Pattern: $pattern\n";

  if (-e $filename) { 
    
    if($filename =~ m/\.fastq\.gz/){
      my $R;
      if($filename =~ m/\S+_L\d\d\d_(R\d)_\d\d\d.fastq\.gz/ || $filename =~ m/_(R\d)\.fastq\.gz/ || $filename =~ m/_(R\d)_\d\d\d\.fastq\.gz/){
        $R = $1;
        #print STDERR $R."\n"; 
      }
      my $pattern;
      my $pattern2;
      my $pattern3;
      #ANT2022-0037_S29_L001_R1_001.fastq.gz
      if($filename =~ m/(\S+)(_L\d\d\d)_R\d(_\d\d\d\.fastq\.gz)/){
        $pattern = $1;
        $pattern2 = $2;
        $pattern3 = $3;
        #print STDERR "1\n";
      }elsif($filename =~ m/(\S+)_R\d(\.fastq\.gz)/){
        $pattern = $1;
        $pattern2 = "";
        $pattern3 = $2;
        #print STDERR "2\n";
      }elsif($filename =~ m/(\S+)_R\d(_\d\d\d\.fastq\.gz)/){
        $pattern = $1;
        $pattern2 = "";
        $pattern3 = $2;
        #print STDERR "3\n";
      }
      #print STDERR $pattern."\n"; 
      my $found_R1 = 0;
      my $found_R2 = 0;
      if(defined $R){
         if($R eq "R1"){
            $found_R1 = 1;
         }elsif($R eq "R2"){
            $found_R2 = 1;
         }else{
            next;
         }
         $hash{$pattern}{prefix}         = $pattern;
         $hash{$pattern}{R1}             = "".$pattern.$pattern2."_R1".$pattern3;
         $hash{$pattern}{R2}             = "".$pattern.$pattern2."_R2".$pattern3;
         #$hash{$pattern}{R1}             = " ".$pattern."_R1.fastq.gz";
         #$hash{$pattern}{R2}             = " ".$pattern."_R2.fastq.gz";
         $hash{$pattern}{libraryBarcode} = "AAAAAAAA";
         $hash{$pattern}{library}        = "MPS000000";
         $hash{$pattern}{run}            = "9999"; 
         $hash{$pattern}{lane}           = "1";
         $hash{$pattern}{adaptor1}       = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
         $hash{$pattern}{adaptor2}       = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG";
         $hash{$pattern}{qualityOffset}  = "33";
         $hash{$pattern}{readset}        = $pattern;
       }
     }
  }
}

find (\&eachFile, $indir);

#print STDERR Dumper(\%hash);

# Then loop through hash to print to STDOUT
print STDOUT "Sample\tReadset\tRunType\tBAM\tFASTQ1\tFASTQ2\tLibrary\tRun\tLane\tAdaptor1\tAdaptor2\tQualityOffset\tBED\n";

for my $key (keys %hash){
   print STDOUT $key."\t".$key."\t$type\t";     
   #print STDOUT $key.".bam\t./raw_reads/".$key."_R1.fastq.gz\t./raw_reads/".$key."_R2.fastq.gz\t".$hash{$key}{libraryBarcode}."\t".$hash{$key}{run}."\t".$hash{$key}{lane}."\t".$hash{$key}{adaptor1}."\t".$hash{$key}{adaptor2}."\t".$hash{$key}{qualityOffset}."\t".$key.".bed\n";
   print STDOUT $key.".bam\t$indir_relpath/".$hash{$key}{R1}."\t$indir_relpath/".$hash{$key}{R2}."\t".$hash{$key}{libraryBarcode}."\t".$hash{$key}{run}."\t".$hash{$key}{lane}."\t".$hash{$key}{adaptor1}."\t".$hash{$key}{adaptor2}."\t".$hash{$key}{qualityOffset}."\t".$key.".bed\n";
}
close(STDOUT);

exit;
