#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
mergeFlagStatsAndQC.pl

PURPOSE:

INPUT:
--infilesFlagstat <string> : list of flagstat files separated by a comma.
--infilesQC <string>       : list of trimmomatic qc (summary) files separated by a comma.
--infilesBBduk <string>    : list of bbduk files separated by a comma.
--infilesBBdukSub <string> : list of bbduk subtraction files separated by a comma. (optional)
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
Assumes both infilesFlagstat, infilesBBduk, infilesBBdukSub and infilesQC are in the exact same order.

AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infilesFlagstat, $infilesQC, $infilesBBduk, $infilesBBdukSub);
my $verbose = 0;

GetOptions(
   'infilesFlagstat=s' => \$infilesFlagstat,
   'infilesQC=s'       => \$infilesQC,
   'infilesBBduk=s'    => \$infilesBBduk,
   'infilesBBdukSub=s' => \$infilesBBdukSub,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## Validate
die("--infilesFlagstat arg needed") unless($infilesFlagstat);
die("--infilesQC arg needed") unless($infilesQC);
die("--infilesBBduk arg needed")    unless($infilesBBduk);    
#die("--infilesBBdukSub arg needed") unless($infilesBBdukSub);

## MAIN
my %hash;
my @infilesFlagstat = split(/,/, $infilesFlagstat);
my @infilesQC = split(/,/, $infilesQC);
my @infilesBBduk = split(/,/, $infilesBBduk);
my @infilesBBdukSub = split(/,/, $infilesBBdukSub);

# Flagstats
foreach my $file (@infilesFlagstat){
   open(IN, "<".$file) or die "Can't open $file\n";
   my $name = $file;
   $name =~ s{.*/}{};      # removes path  
   $name =~ s{\.[^.]+$}{};    # removes extension

   print STDERR $name."\n";

   while(<IN>){
      #1,5,9
      if($. == 1){
         my($totalReads) = $_ =~ m/^(\d+)/;
         $hash{$name}{total_reads} = $totalReads;
      }if($. == 5){
         my($mapped) = $_ =~ m/^(\d+)/;
         if($mapped == 0){
            $hash{$name}{mapped} = 0;
            $hash{$name}{mapped_perc} = 0;
        }else{
            $hash{$name}{mapped} = $mapped;
            $hash{$name}{mapped_perc} = int(($mapped/$hash{$name}{total_reads}) * 100);
         }
      }if($. == 9){
         my($properlyPaired) = $_ =~ m/^(\d+)/;
         if($properlyPaired == 0){
            $hash{$name}{properly_paired} = 0;
            $hash{$name}{properly_paired_perc} = 0;
         }else{
            $hash{$name}{properly_paired} = $properlyPaired;
            $hash{$name}{properly_paired_perc} = int(($properlyPaired/$hash{$name}{total_reads}) * 100);
         }
      }
   }
   close(IN);
}

# --infilesFlagstat 
# CLA17-cyano/contigs_abundance/megahit/      CLA17-cyano_S52_L001.flagstats,
# CLA17-cyano/contigs_abundance/megahit/      CLA17-cyano_S53_L001.flagstats --infilesQC i
# CLA17-cyano/qced_reads/CLA17-cyano_S52_L001/CLA17-cyano_S52_L001_trim_stats.csv,
# CLA17-cyano/qced_reads/CLA17-cyano_S53_L001/CLA17-cyano_S53_L001_trim_stats.csv

# QC
foreach my $file (@infilesQC){
   open(IN, "<".$file) or die "Can't open $file\n";
   my $name = $file;
   $name =~ s{.*/}{};         # removes path  
   $name =~ s{\.[^.]+$}{};    # removes extension
   $name =~ s/\.trim\.stats//;
   $name =~ s/_trim_stats//;
   print STDERR "Name QC file: ".$name."\n";

   while(<IN>){
      my $line = $_;
      my @row = split(/,/, $_);
      my $value = $row[1];
      #1,2,3
      #print STDERR $value."\n";
      if($. == 1){
         my($raw_frag) = $value =~ m/(\d+)/;
         $hash{$name}{raw_frag} = $raw_frag;
      }if($. == 2){
         my($frag_surv) = $value =~ m/(\d+)/;
         $hash{$name}{frag_surv} = $frag_surv;
         $hash{$name}{frag_surv_perc} = int(($frag_surv/$hash{$name}{raw_frag}) * 100);
      }if($. == 3){
         my($single_surv) = $value =~ m/(\d+)/;
         $hash{$name}{single_surv} = $single_surv;
         $hash{$name}{single_surv_perc} = int(($single_surv/$hash{$name}{raw_frag}) * 100);
      }
   }
   close(IN);
}

## Infiles BBduk
foreach my $file (@infilesBBduk){
   my $prefix = $file;
   $prefix =~ s{.*/}{};         # removes path  
   $prefix =~ s{\.[^.]+$}{};    # removes extension
   $prefix =~ s/\.duk_contam_interleaved_log\.txt//;
   $prefix =~ s/\.duk_contam_interleaved_log//;
   print STDERR "Name BBDUK contam file: ".$prefix."\n";

   open(IN, "<".$file) or die "Can't open $file\n";
   while(<IN>){
      chomp;
      #print STDERR $_."\n";
      #my($matched) = $_ =~ m/Total number of matched reads:\s+(\d+)/g;
      # duk Total number of reads:    10778900
      # duk Total number of matched reads: 481670
      # bbduk:
      # #Total\t
      # #Matched\t
      #if($_ =~ m/Total number of matched reads:\s+(\d+)/){
      if($_ =~ m/^#Total\t(\d+)/){
          #print STDERR $1."\n";
         $hash{$prefix}{TOTAL} = $1;
      }
      #if($_ =~ m/Total number of reads:\s+(\d+)/){
      if($_ =~ m/^#Matched\t(\d+)\t(\S+)/){
          #print STDERR $1."\n";
         $hash{$prefix}{MATCH} = $1;
         $hash{$prefix}{RATIO} = $2;
      }
      #if($_ =~ m/Match ratio:\s+(\d+\.\d+)/){
      #    #print STDERR $1."\n";
      #   $hash1{$prefix}{RATIO} = $1;
      #}
   }
   close(IN);
}

## Infiles BBduk sub
if($infilesBBdukSub){
    foreach my $file (@infilesBBdukSub){
        my $prefix = $file;
        $prefix =~ s{.*/}{};         # removes path  
        $prefix =~ s{\.[^.]+$}{};    # removes extension
        $prefix =~ s/\.ncontam_paired_sub_log\.txt//;
        $prefix =~ s/\.ncontam_paired_sub_log//;
        print STDERR "Name BBDUK subtract file: ".$prefix."\n";

       my $k = 0;
       open(IN, "<".$file) or die "Can't open $file\n";
       while(<IN>){
          chomp;
          #print STDERR $_."\n";
          #my($matched) = $_ =~ m/Total number of matched reads:\s+(\d+)/g;
          # duk Total number of reads:    10778900
          # duk Total number of matched reads: 481670
          # bbduk:
          # #Total\t
          # #Matched\t
          #if($_ =~ m/Total number of matched reads:\s+(\d+)/){
          if($_ =~ m/^Reads Used:\s+\t(\d+)\t/){
              #print STDERR $1."\n";
             $hash{$prefix}{TOTAL} = $1;
          }

          #if($_ =~ m/Total number of reads:\s+(\d+)/){
          if($_ =~ m/^Read 1 data:\s+\t/){
              $k = 1;
          }elsif($_ =~ m/^Read 2 data:\s+\t/){
              $k = 2;
          }
          
          if($_ =~ m/mapped:\s+\t/){
              my @row = split(/\t/, $_);
              $hash{$prefix}{"R".$k}{RATIO} = $row[1];
              $hash{$prefix}{"R".$k}{MATCH} = $row[2];
          }
       }
       close(IN);
    }
}

print STDOUT "sampleName\trawFragments\tsurvivingFragments\tsurvivingFragments%\tsurvivingSingle\t";
print STDOUT "contamQCedReads\tsubtractedReads\ttotalQCedReads\tratioQCedReads\tmapped\tmapped%\tproperlyPaired\tproperlyPaired%\n";
for my $key (sort keys %hash) {
    print STDOUT $key."\t";
    print STDOUT commify($hash{$key}{raw_frag})."\t";
    print STDOUT commify($hash{$key}{frag_surv})."\t";
    print STDOUT $hash{$key}{frag_surv_perc}."%\t";
    print STDOUT commify($hash{$key}{single_surv})."\t";
    #print STDOUT $hash{$key}{single_surv_perc}."\t";
    print STDOUT commify($hash{$key}{MATCH})."\t";
    if(exists $hash{$key}{R1}{MATCH}){
        my $R1 = $hash{$key}{R1}{MATCH};
        my $R2 = $hash{$key}{R2}{MATCH};
        my $subtractedReads = $R1 + $R2;
        print STDOUT commify($subtractedReads)."\t";
    }else{
        print STDOUT "NA\t";
    }
    print STDOUT commify($hash{$key}{total_reads})."\t"; # reads that enter the mapping against ref step.
    # compute ratio QCed passed reads
    my $ratio = $hash{$key}{total_reads} / ($hash{$key}{frag_surv} * 2) * 100;
    $ratio = int($ratio);
    print STDOUT $ratio."%\t";
    print STDOUT commify($hash{$key}{mapped})."\t";
    print STDOUT $hash{$key}{mapped_perc}."%\t";
    print STDOUT commify($hash{$key}{properly_paired})."\t";
    print STDOUT $hash{$key}{properly_paired_perc}."%\n";
}

sub commify {
    scalar reverse join ',',
    unpack '(A3)*',
    reverse scalar shift;
}
