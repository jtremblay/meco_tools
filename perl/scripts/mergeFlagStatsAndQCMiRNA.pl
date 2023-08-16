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
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
Assumes both infilesFlagstat and infilesQC are in the exact same order.

AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infilesFlagstat, $infilesQC);
my $verbose = 0;

GetOptions(
   'infilesFlagstat=s' => \$infilesFlagstat,
   'infilesQC=s'       => \$infilesQC,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my @infilesFlagstat = split(/,/, $infilesFlagstat);
my @infilesQC = split(/,/, $infilesQC);

# Flagstats
foreach my $file (@infilesFlagstat){
   open(IN, "<".$file) or die "Can't open $file\n";
   my $name = $file;
   $name =~ s{.*/}{};      # removes path  
   $name =~ s{\.[^.]+$}{};    # removes extension
   $name =~ s{\.[^.]+$}{};    # removes extension

   #print STDERR $name."\n";

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
# #Input Reads: 30941049 Surviving: 27145431 (87.73%) Dropped: 3795618 (12.27%)
foreach my $file (@infilesQC){
   open(IN, "<".$file) or die "Can't open $file\n";
   my $name = $file;
   $name =~ s{.*/}{};         # removes path  
   $name =~ s{\.[^.]+$}{};    # removes extension
   $name =~ s/\.trim\.stats//;
   $name =~ s/_trim_stats//;
   #print STDERR "Name QC file: ".$name."\n";

   while(<IN>){
      my $line = $_;

      if($line =~ m/Input Reads\: (\d+) Surviving\: (\d+) \((\S+)\%\) Dropped: (\d+) \((\S+)\%\)$/){
        $hash{$name}{raw_frag} = $1;
        $hash{$name}{frag_surv} = $2;
        $hash{$name}{frag_surv_perc} = $3;
        #print STDERR "1: $1\n";
        #print STDERR "2: $2\n";
        #print STDERR "3: $3\n";
      }
   }
   close(IN);
}

print STDOUT "sampleName\trawFragments\tsurvivingFragments\tsurvivingFragments%\tpassedSizeSelect\tmapped\tmapped%\n";
for my $key (sort keys %hash) {
    print STDOUT $key."\t";
    print STDOUT commify($hash{$key}{raw_frag})."\t";
    print STDOUT commify($hash{$key}{frag_surv})."\t";
    print STDOUT $hash{$key}{frag_surv_perc}."%\t";
    #print STDOUT commify($hash{$key}{single_surv})."\t";
    #print STDOUT $hash{$key}{single_surv_perc}."\t";
    print STDOUT commify($hash{$key}{total_reads})."\t";
    print STDOUT commify($hash{$key}{mapped})."\t";
    print STDOUT $hash{$key}{mapped_perc}."%\n";
    #print STDOUT commify($hash{$key}{properly_paired})."\t";
    #print STDOUT $hash{$key}{properly_paired_perc}."%\n";
}

sub commify {
    scalar reverse join ',',
    unpack '(A3)*',
    reverse scalar shift;
}
