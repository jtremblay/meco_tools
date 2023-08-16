#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Parallel::ForkManager;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
bwaMemSamtoolsMga.pl

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
my ($help, $infiles1, $infiles2, $num_threads, $reference, $outfile_mapped, $ram, $min_id);
my $verbose = 0;

GetOptions(
   'infiles1=s'        => \$infiles1,
   'infiles2=s'        => \$infiles2,
   #'tmp_bams=s'        => \$tmp_bams,
   'ram=s'             => \$ram,
   'min_id=f'          => \$min_id,
   'outfile_mapped=s'  => \$outfile_mapped,
   'num_threads=i'     => \$num_threads,
   'reference=s'       => \$reference,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my @infiles1 = split(/,/, $infiles1);
my @infiles2 = split(/,/, $infiles2);
#my @tmp_bams = split(/,/, $tmp_bams);
my $internal_num_threads = 1;
#my $mem_per_thread = "3500M";
#my @fastqs_pairedmapped;
#my @bams_R1mapped;
#my @bams_R2mapped;
my @fastqs_pairedmapped;
$min_id = 0.98 unless($min_id);

#my $reference_basename = $reference;
#$reference_basename =~ s{^.*/}{};  
#$reference_basename =~ s{\.[^.]+$}{};

my ($ref_filename, $ref_dir) = fileparse($reference);
$ref_filename =~ s{^.*/}{};  
$ref_filename =~ s{\.[^.]+$}{};
print STDERR "ref_filename: ".$ref_filename."\n";

my $cmd = "bbmap.sh ref=$reference path=$ref_dir/$ref_filename/";
system($cmd);
die "command failed: $!\n" if($? != 0);


my @cmds;
foreach my $infile1 (@infiles1){
    my $infile2 = shift(@infiles2);

    my $name = $infile1;
    $name =~ s{^.*/}{};  
    $name =~ s{\.[^.]+$}{};
    $name =~ s/_R1.fastq//;

    my ($filename, $dir) = fileparse($outfile_mapped);
    
    #print STDERR "Name: ".$name."\n";
    #print STDERR "Dir: ".$dir."\n";
    #exit(1);
 
    print STDERR "directories: ".$dir."\n";
    print STDERR "fastqs basename : ".$name."\n";

    #bbmap.sh -threads=1 -Xmx50g in=./TN14GW.ncontam_paired_R1.fastq.gz in2=./TN14GW.ncontam_paired_R2.fastq.gz out=mapped.sam path=./assembly/  bamscript=bs.sh; sh bs.sh
    my $cmd = "bbmap.sh";
    $cmd .= " -threads=$internal_num_threads -Xmx$ram";
    $cmd .= " -in=$infile1 -in2=$infile2 outm=$dir/$name.fq path=$ref_dir/$ref_filename minid=$min_id";

    push(@cmds, $cmd);
    push(@fastqs_pairedmapped, $dir."/".$name.".fq");
}

my $pm = new Parallel::ForkManager($num_threads); 

foreach my $command (@cmds) {
    $pm->start and next; # do the fork
    print STDERR "executing: $command\n";
    system($command);
    $pm->finish; # do the exit in the child process
}
$pm->wait_all_children;

#Once done merge all bams into one large bam and delete the rest
#
#my $cmd_merge = "samtools merge -f ".$outfile_mapped." ".join(" ", @fastqs_pairedmapped);
my $cmd_merge = "cat ".join(" ", @fastqs_pairedmapped)." > ".$outfile_mapped;
print STDERR "Attemtping to execute following cmd: $cmd_merge\n";
system($cmd_merge);
die "command failed: $!\n" if($? != 0);

my $cmd_del = "rm ".join(" ", @fastqs_pairedmapped);
print STDERR "Attemtping to execute following cmd: $cmd_del\n";
system($cmd_del);
die "command failed: $!\n" if($? != 0);

exit(0);
