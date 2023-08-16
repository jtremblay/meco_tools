#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use Cwd;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
generateSymLinksMG.pl

PURPOSE:
Creates symlinks for the qced_reads directory when
merging two MG projects for a single analysis.

INPUT
--indir <string> : MG project directory
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
   'indir=s' => \$indir,
   'verbose' => \$verbose,
   'help'    => \$help
);
if ($help) { print $usage; exit; }

die "--indir arg missing\n" unless($indir);

## SUB
my $cwd;

sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) { 
   
    #if($filename =~ m/\.ncontam_paired_R1\.fastq\.gz/){
    # print STDERR $fullpath."\n";
    #}

    if($fullpath =~ m/\/qced_reads\/(\S+)\/(\S+\.ncontam_paired_R[1|2]\.fastq\.gz)/){
      print STDERR "Found ".$1."/".$2."\n";
      system("mkdir -p ".$cwd."/".$1);
      #symlink(TARGET, NAME);
      symlink($fullpath, $cwd."/".$1."/".$2);
      $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly created symbolic symlink ".$fullpath." into ".$cwd."/".$1."/".$2."\n";
    }
    
    if($fullpath =~ m/\/qced_reads\/(\S+)\/(\S+\.trim\.stats\.csv)/){
      print STDERR "Found ".$1."/".$2."\n";
      system("mkdir -p ".$cwd."/".$1);
      #symlink(TARGET, NAME);
      symlink($fullpath, $cwd."/".$1."/".$2);
      $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly created symbolic symlink ".$fullpath." into ".$cwd."/".$1."/".$2."\n";
    }
  }
}

## MAIN

$indir = abs_path($indir);
$indir = $indir."/qced_reads/";

$cwd = getcwd();
$cwd = abs_path($cwd);
system("mkdir -p ".$cwd."/qced_reads");
$cwd = $cwd."/qced_reads";

# Compress .fastq into .gz
find (\&eachFile, $indir);
