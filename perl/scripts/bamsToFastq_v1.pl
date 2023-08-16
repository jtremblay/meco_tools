#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Which;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
bamToFastq.pl

PURPOSE:

INPUT:
--infiles <string>    : list of bam files separated by a ,
				
OUTPUT:
--outfile1 <string> : fastq R1  
--outfile2 <string> : fastq R2

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $outfile1, $outfile2);
my $verbose = 0;

GetOptions(
   'infiles=s'    => \$infiles,
   'outfile1=s'   => \$outfile1,
   'outfile2=s'   => \$outfile2,
   'verbose'      => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $bamToFastq = which('bamToFastq'); chomp $bamToFastq;
print STDERR "path of bamToFastq:\t".$bamToFastq."\n";;
die "Can't find bamToFastq on path\n" if(!defined($bamToFastq));

print STDERR $outfile1."\n".$outfile2."\n";
my $tmpR1 = $outfile1.".tmp";
my $tmpR2 = $outfile2.".tmp";
print STDERR $tmpR1."\n".$tmpR2."\n";
system("echo '' > $outfile1");
system("echo '' > $outfile2");

my @infiles = split(/,/, $infiles);

foreach my $file (@infiles){
    my $cmd = "bamToFastq -bam ".$file." -fq1 ".$tmpR1." -fq2 ".$tmpR2." -useTags\n";
    print STDERR "Executing: ".$cmd."\n";
    system($cmd);
    die "command failed: $!\n" if($? != 0);

    # Then add results to final fastq file.
    my $cmd_cat_R1 = "cat $tmpR1 >> $outfile1";
    print STDERR $cmd_cat_R1."\n";
    system($cmd_cat_R1);
    die "command failed: $!\n" if($? != 0);
    
    my $cmd_cat_R2 = "cat $tmpR2 >> $outfile2";
    print STDERR $cmd_cat_R2."\n";
    system($cmd_cat_R2);
    die "command failed: $!\n" if($? != 0);

    system("rm $tmpR1 $tmpR2 -rf");
    die "command failed: $!\n" if($? != 0);

}
exit;
