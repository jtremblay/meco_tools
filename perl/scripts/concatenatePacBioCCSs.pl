#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use File::Find;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
concatenatePacBioCCSs.pl

PURPOSE:
Writes a cat script to obtain one fastq per file. Generates a readset as well
(to use with a pipeline). Assumes that raw reads are all located in the directory
./raw_data. In the --cat output, output reads will be concatenated in ./raw_reads dir.

INPUT:
--infile <string>  : Text file where each line represent a unique sample name. 
                     Also make sure that reads downloaded from nanuq are located in 
                     ./raw_data/ directory in your current location. (Optional).
                     If no --infile arg is specified, will determine automatically 
                     sample names.
--indir <string>   : Directory were are located pacbio zip fastqs to be concatenated.

OUTPUT:
--readset <string> : File into which will be written read sets
                     (mandatory for pipeline usage).
--makeSL <string>  : Optional. Directory in which to make symlinks.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $readset, $makeSL, $indir);
my $verbose = 0;

GetOptions(
    'infile=s'  => \$infile,
    'indir=s'   => \$indir,
    'readset=s' => \$readset,
    'makeSL=s'  => \$makeSL,
    'verbose'   => \$verbose,
    'help'      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
#die "missing --infile\n" unless($infile);
die "missing --readset\n" unless($readset);
#die "missing --makeSL\n" unless($makeSL);

## MAIN
my %hash1;
my %hash2;
my %hash;
my $header;
#my $indir = abs_path("./raw_data/");
$indir = abs_path($indir);
my $pattern;
my $prefix;
my %hash_samples;

system("mkdir -p ./raw_reads");

sub eachFile{
  my $filename = $_;
  #my $pattern = shift @_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_
  #print STDERR "Pattern: $pattern\n";

  if (-e $filename) { 
    
    if($filename =~ m/$pattern/){
      #print STDERR "Found ".$fullpath." ...\n";
      my($R) = $filename =~ /_(R\d)\.fastq\.gz/;
      #print STDERR "R: $R\n";
      my $found_R1 = 0;
      my $found_R2 = 0;
      if($R eq "R1"){
         $hash1{$pattern}{$filename} = $fullpath;
         $found_R1 = 1;
      }elsif($R eq "R2"){
         $hash2{$pattern}{$filename} = $fullpath;
         $found_R2 = 1;
      }
      #if($found_R1 == 0){
      #   print STDERR "R1 file missing for ".$prefix."\n";
      #}
      #if($found_R2 == 0){
      #   print STDERR "R2 file missing for ".$prefix."\n";
      #}
      $hash{$prefix}{prefix}         = $prefix;
      $hash{$prefix}{R1}             = " ".$filename."_R1.fastq.gz";
      $hash{$prefix}{R2}             = " ".$filename."_R2.fastq.gz";
      $hash{$prefix}{libraryBarcode} = "AAAAAAAA";
      $hash{$prefix}{library}        = "MPS000000";
      $hash{$prefix}{run}            = "9999"; 
      $hash{$prefix}{lane}           = "1";
      $hash{$prefix}{adaptor1}       = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
      $hash{$prefix}{adaptor2}       = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG";
      $hash{$prefix}{qualityOffset}  = "33";
      $hash{$prefix}{readset}        = $prefix;
    }
  }
}

sub eachFileName{
    my $filename = $_;
    #my $pattern = shift @_;
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_
    #print STDERR "Pattern: $pattern\n";

    if (-e $filename) { 
    
        if($filename =~ m/\.ccs.fastq.zip/){
            my $R;
            #if($filename =~ m/_(R\d)\.fastq\.gz/){
                #print STDERR "Found ".$fullpath." ...\n";
                $R = "R1";
                #print STDERR "R: $R\n";
                my $found_R1 = 0;
                my $found_R2 = 0;
                #.Mal_39.ccs.fastq.zip
                my $pattern;
                if($filename =~ m/_forward\.(\S+)\.ccs\.fastq\.zip/){
                    print STDERR "match: $1\n";
                    $pattern = $1;
                    if($R eq "R1"){
                        $hash1{$pattern}{$filename} = $fullpath;
                        $found_R1 = 1;
                        $prefix = $pattern;
                        $pattern =~ s/_R1//g;
                    }elsif($R eq "R2"){
                        $hash2{$pattern}{$filename} = $fullpath;
                        $found_R2 = 1;
                        $prefix = $pattern;
                        $pattern =~ s/_R2//g;
                    }
                    $hash{$prefix}{prefix}         = $prefix;
                    $hash{$prefix}{R1}             = " ".$filename."_R1.fastq.gz";
                    $hash{$prefix}{R2}             = " ".$filename."_R2.fastq.gz";
                    $hash{$prefix}{libraryBarcode} = "AAAAAAAA";
                    $hash{$prefix}{library}        = "MPS000000";
                    $hash{$prefix}{run}            = "9999"; 
                    $hash{$prefix}{lane}           = "1";
                    $hash{$prefix}{adaptor1}       = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
                    $hash{$prefix}{adaptor2}       = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG";
                    $hash{$prefix}{qualityOffset}  = "33";
                    $hash{$prefix}{readset}        = $prefix;
                }
            #}
        }elsif($filename =~ m/_R\d\.fastq\.gz/){
            #print STDERR "Found ".$fullpath." ...\n";
            my $pattern;
            my $R;
            if($filename =~ /(.*)_(R\d)\.fastq\.gz/){
                $pattern = $1;
                $R = $2;
            }
            #print STDERR "R: $R\n";
            my $found_R1 = 0;
            my $found_R2 = 0;
            if($R eq "R1"){
                $hash1{$pattern}{$filename} = $fullpath;
                $found_R1 = 1;
            }elsif($R eq "R2"){
                $hash2{$pattern}{$filename} = $fullpath;
                $found_R2 = 1;
            }
            $hash{$pattern}{pattern}         = $pattern;
            $hash{$pattern}{R1}             = " ".$filename."_R1.fastq.gz";
            $hash{$pattern}{R2}             = " ".$filename."_R2.fastq.gz";
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

# one sample name per row.
if($infile){
    open(IN, "<".$infile);
    while(<IN>){
       chomp;
       $pattern = $_."_";
       $prefix = $_;
       #print STDERR "Pattern ".$_."\n";
       $hash_samples{$pattern} = $pattern;
       find (\&eachFile, $indir);
    }
    close(IN);
}else{
    find (\&eachFileName, $indir);
}
#exit;
#close(IN);
#print STDERR Dumper(%hash_samples);
#exit;


#R1
my $cat;
for my $k1 (sort keys %hash1) {
   $cat = "zcat ";
   for my $k2 (sort keys %{ $hash1{$k1} }){
       #print STDERR $hash1{$k1}{$k2}."\n";
      $cat .= $hash1{$k1}{$k2} ." ";
   }
   my $prefix = substr($k1, 0);
   $cat .= "> ./raw_reads/".$prefix."_R1.fastq && gzip ./raw_reads/".$prefix."_R1.fastq";
   print STDOUT $cat."\n";
}

#R2
undef $cat;
for my $k1 (sort keys %hash2) {
   $cat = "cat ";
   for my $k2 (sort keys %{ $hash2{$k1} }){
       #print STDERR $hash2{$k1}{$k2}."\n";
      $cat .= $hash2{$k1}{$k2} ." ";
   }
   my $prefix = substr($k1, 0, -1);
   $cat .= "> ./raw_reads/".$prefix."_R2.fastq.gz";
   print STDOUT $cat."\n";
}
#print STDERR Dumper(\%hash1);

# Then loop through hash to print to STDOUT
open(READSET, ">".$readset) or die "Can't open $readset\n";
print READSET "Sample\tReadset\tRunType\tBAM\tFASTQ1\tFASTQ2\tLibrary\tRun\tLane\tAdaptor1\tAdaptor2\tQualityOffset\tBED\n";

system("mkdir -p ".$makeSL) if($makeSL);

for my $key (keys %hash){
   print READSET $key."\t".$key."\tSINGLE_END\t"; 
    
   print READSET $key.".bam\t./raw_reads/".$key."_R1.fastq.gz\t./raw_reads/".$key."_R2.fastq.gz\t".$hash{$key}{libraryBarcode}."\t".$hash{$key}{run}."\t".$hash{$key}{lane}."\t".$hash{$key}{adaptor1}."\t".$hash{$key}{adaptor2}."\t".$hash{$key}{qualityOffset}."\t".$key.".bed\n";

   if($makeSL){
      symlink(abs_path("./raw_reads/".$key."_R1.fastq.gz"), $makeSL."/".$key."_R1.fastq.gz"); 
      symlink(abs_path("./raw_reads/".$key."_R2.fastq.gz"), $makeSL."/".$key."_R2.fastq.gz"); 
   }
   
}
close(READSET);

for my $key (keys %hash_samples){
    if(exists $hash1{$key}){
        #print STDERR "$key Exists!\n";
    }else{
        print STDERR "$key does not exists! - R1\n";
    }
    if(exists $hash2{$key}){
        #print STDERR "$key Exists!\n";
    }else{
        print STDERR "$key does not exists! - R2\n";
    }
}

exit;
