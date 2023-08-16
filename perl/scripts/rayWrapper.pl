#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
rayWrapper.pl

PURPOSE:

INPUT:
--indir <string>  : directory where fastq.gz files are located
--range <string>  : range of canopies to assemble.

OUTPUT:
--outdir <string> : directory where assemblies will be located.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $outdir, $num_threads, $kmer, $minContigLength, $range);
my $verbose = 0;

GetOptions(
   'indir=s' 	         => \$indir,
   'outdir=s' 	         => \$outdir,
   'num_threads=s' 	   => \$num_threads,
   'kmer=s' 	         => \$kmer,
   'range=s'            => \$range,
   'minContigLength=s' 	=> \$minContigLength,
   'verbose' 	         => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my @range = split(/-/, $range);
my $start = $range[0];
my $end = $range[1];

die "start has to be higher than end...\n" if($end <= $start);

my %hash;

sub eachFile{
   my $filename = $_;
   my $fullpath = $File::Find::name;
   #remember that File::Find changes your CWD, 
   #so you can call open with just $_

   if (-e $filename) { 
      
      if(substr($filename, -9) eq ".fastq.gz"){
         if($filename =~ m/canopy(\d{5})\.R1\.fastq\.gz/){
            $hash{$1}{R1paired} = $fullpath;

         }elsif($filename =~ m/canopy(\d{5})\.R2\.fastq\.gz/){
            $hash{$1}{R2paired} = $fullpath;
         
         }elsif($filename =~ m/canopy(\d{5})\.R1\.unpaired\.fastq\.gz/){
            $hash{$1}{R1unpaired} = $fullpath;

         }elsif($filename =~ m/canopy(\d{5})\.R2\.unpaired\.fastq\.gz/){
            $hash{$1}{R2unpaired} = $fullpath;

         }
      }
   }
}

# build hash structure for ray input.
find (\&eachFile, $indir);

#for my $key (sort{$a cmp $b} keys %hash){
#   my $R1unpaired = $hash{$key}{R1unpaired};
#   my $R2unpaired = $hash{$key}{R2unpaired};
#
#   if(-e "$outdir/$key.singlets.fastq.gz" && -s "$outdir/$key.singlets.fastq.gz"){
#      # Do nothin... file exists and is not empty.
#   }else{
#      print STDERR "cat $R1unpaired $R2unpaired > $outdir/$key.singlets.fastq.gz\n";
#      system("cat $R1unpaired $R2unpaired > $outdir/$key.singlets.fastq.gz");
#      $hash{$key}{singlets} = "$outdir/$key.singlets.fastq.gz";
#   }
#}

#Then build commands.
for my $key (sort{$a cmp $b} keys %hash){

   if($key >= $start && $key <= $end){  

      my $currOutdir = $outdir."/".$key;
   
      system("rm $currOutdir -rf");

      unless(-e "$currOutdir/$key.done"){

         my $R1paired = $hash{$key}{R1paired};
         my $R2paired = $hash{$key}{R2paired};
         my $R1unpaired = $hash{$key}{R1unpaired};
         my $R2unpaired = $hash{$key}{R2unpaired};

         #system("mkdir -p $currOutdir");

         my $cmd .= "mpiexec -n $num_threads Ray -k $kmer -minumum-contig-length $minContigLength -p $R1paired $R2paired -s $R1unpaired -s $R2unpaired -o $currOutdir/";
         print STDERR $cmd."\n";
         system($cmd);
         if($? == 0){
            print STDERR "Successfully assembled $key...\n";
            system("touch $currOutdir/$key.done");
         }else{
            die "Failed to assemble $key... command failed: $!\n";
         }
      }
   }
}


