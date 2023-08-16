#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
pacBioGenerateFofns.pl

PURPOSE:

INPUT:
--indir <string>       : Indir where are located the fofns

OUTPUT:
--outdir <string>      : outdir where will be written the subreads for
                         that specific movie in the fofn.

NOTES:


BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $outdir);
my $verbose = 0;

GetOptions(
   'indir=s'   => \$indir,
   'outdir=s'  => \$outdir,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN


#mkdir $outdir unless(-d $outdir);
opendir(DIR, $indir) or die "Could not open $indir\n";
print STDERR $indir."\n";
while (my $filename = readdir(DIR)) {
    print STDERR "$filename\n";
    $filename = $indir."/".$filename;
    print STDERR "$filename\n" if -f $filename;
    if(-f $filename){
        if($filename =~ m/.*\.txt$/){
            my $cmd = "bax2bam -f ".$filename." -o ".$outdir." --output-xml=".$outdir.".xml --subread --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag";
            print STDERR "Executing: $cmd\n";
            system($cmd);
            die "command failed: $!\n" if($? != 0);   
        }
    }
}

exit;

