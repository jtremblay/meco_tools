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
--infile_paired <string>    : bam file
--infile1 <string>          : bam file
--infile2 <string>          : bam file
				
OUTPUT:
--outfile_paired <string>   : bam file
--outfile1 <string>         : fastq R1  
--outfile2 <string>         : fastq R2

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_paired, $infile1, $infile2, $outfile_paired1, $outfile_paired2, $outfile1, $outfile2);
my $verbose = 0;

GetOptions(
   'infile_paired=s'   => \$infile_paired,
   'infile1=s'         => \$infile1,
   'infile2=s'         => \$infile2,
   'outfile_paired1=s' => \$outfile_paired1,
   'outfile_paired2=s' => \$outfile_paired2,
   'outfile1=s'        => \$outfile1,
   'outfile2=s'        => \$outfile2,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $bamToFastq = which('bamToFastq'); chomp $bamToFastq;
print STDERR "path of bamToFastq:\t".$bamToFastq."\n";;
die "Can't find bamToFastq on path\n" if(!defined($bamToFastq));

my $tmpR1 = $outfile1.".tmp";
my $tmpR2 = $outfile2.".tmp";
system("echo '' > $outfile1");
system("echo '' > $outfile2");

my $cmd = "bamToFastq -bam ".$infile_paired." -fq1 ".$outfile_paired1." -fq2 ".$outfile_paired2;
print STDERR "Executing: ".$cmd."\n";
system($cmd);
die "command failed: $!\n" if($? != 0);

$cmd = "";
$cmd = "bamToFastq -bam ".$infile1." -fq1 ".$outfile1." -fq2 /dev/null";
print STDERR "Executing: ".$cmd."\n";
system($cmd);
die "command failed: $!\n" if($? != 0);

$cmd = "";
$cmd = "bamToFastq -bam ".$infile2." -fq1 /dev/null -fq2 ".$outfile2;
print STDERR "Executing: ".$cmd."\n";
system($cmd);
die "command failed: $!\n" if($? != 0);

exit;
