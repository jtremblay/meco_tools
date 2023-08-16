#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
addSymlinkForQCedReads.pl

PURPOSE:

INPUT:
--indir <string>  : directory containing qced reads
--outdir <string> : outdir where symlinks should be written
                    Should be typically ./project_root/qced_reads/   
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $outdir);
my $verbose = 0;

GetOptions(
   'indir=s' 	=> \$indir,
   'outdir=s'  => \$outdir,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

$indir = abs_path($indir);
$outdir = abs_path($outdir);

system("mkdir -p $outdir\n");

## MAIN
sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if($filename =~ m/(\S+)\.ncontam_paired_R1\.fastq\.gz/){
          #print STDERR "Found fullpath ".$fullpath."\n";
          #print STDERR "Found filename ".$filename."\n";

         my $subdir = $1;
         #print STDERR "Found filename :".$outdir."/$subdir/".$filename."\n";
         system("mkdir -p $outdir/$subdir\n");

         symlink($fullpath, $outdir."/$subdir/".$filename);
         #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
		}
		if($filename =~ m/(\S+)\.ncontam_paired_R2\.fastq\.gz/){
          #print STDERR "Found fullpath ".$fullpath."\n";
          #print STDERR "Found filename ".$filename."\n";

         my $subdir = $1;
         #print STDERR "Found filename :".$outdir."/$subdir/".$filename."\n";
         system("mkdir -p $outdir/$subdir\n");

         symlink($fullpath, $outdir."/$subdir/".$filename);
         #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
		}
		if($filename =~ m/(\S+)\.trim\.stats\.csv/){
          #print STDERR "Found fullpath ".$fullpath."\n";
          #print STDERR "Found filename ".$filename."\n";

         my $subdir = $1;
         #print STDERR "Found filename :".$outdir."/$subdir/".$filename."\n";
         system("mkdir -p $outdir/$subdir\n");

         symlink($fullpath, $outdir."/$subdir/".$filename);
         #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
		}
	}
}

## MAIN

## Compress .fastq into .gz
find (\&eachFile, $indir);

exit;
