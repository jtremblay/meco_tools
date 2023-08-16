#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Cwd 'abs_path';
use File::Find;
use Cwd qw(getcwd);
use File::Spec;
use Data::Dumper;
use Cwd;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;


my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
   'indir=s'  => \$indir,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }
my $cwd = getcwd;

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-FastqsGzToFastas-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

sub eachFile{
    my $filename = $_; 
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_
     
    if (-e $filename) { 
        #if(substr($filename, -(length($extension)+1)) eq $extension ){
        if($filename =~ m/.*R1_001\.fastq.gz/){
            #print STDERR "Filename: ".$filename."\n";
            #print STDERR "Fullpath: ".$fullpath."\n";
            my $pipe = "$tmpdir/reads.pipe";
            system("mkfifo $pipe");
            system("gunzip -c ".$fullpath." > $pipe &");

            my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fastq file, $pipe\n");
            #my $i = 0;
            #while( my $curr = $ref_fastq_db->next_seq()) {
            #print STDERR $curr->header;
            my $curr = $ref_fastq_db->next_seq();
            if($curr->header =~ m/^\@M00270:228/){
                print STDOUT $fullpath."\n";
            }                    
            #    $i++;
            #}
            system("rm $pipe");
        }
    }
}

## MAIN

find(\&eachFile, $indir);

