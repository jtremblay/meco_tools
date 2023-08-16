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
my ($help, $infiles1, $infiles2, $num_threads, $reference, $tmp_bams, $outfile_paired, $outfile_R1, $outfile_R2);
my $verbose = 0;

GetOptions(
   'infiles1=s'        => \$infiles1,
   'infiles2=s'        => \$infiles2,
   'tmp_bams=s'        => \$tmp_bams,
   'outfile_paired=s'  => \$outfile_paired,
   'outfile_R1=s'      => \$outfile_R1,
   'outfile_R2=s'      => \$outfile_R2,
   'num_threads=i'     => \$num_threads,
   'reference=s'       => \$reference,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my @infiles1 = split(/,/, $infiles1);
my @infiles2 = split(/,/, $infiles2);
my @tmp_bams = split(/,/, $tmp_bams);
my $internal_num_threads = 1;
my $mem_per_thread = "3500M";
my @bams_pairedmapped;
my @bams_R1mapped;
my @bams_R2mapped;

my @cmds;
foreach my $infile1 (@infiles1){
    my $infile2 = shift(@infiles2);
    my $tmp_bam = shift(@tmp_bams);

    my $name = $tmp_bam;
    $name =~ s{^.*/}{};  
    $name =~ s{\.[^.]+$}{};

    my($filename, $dir) = fileparse($tmp_bam);

    print STDERR "directories: ".$dir."\n";

    #-q 1 -F 4 - > $tmp_bam.tmp &&";
    my $cmd = "bwa mem -M";
    $cmd .= " -t $internal_num_threads";
    $cmd .= " $reference";
    $cmd .= " $infile1 $infile2";
    $cmd .= " | samtools view -Sbh - > $tmp_bam.tmp &&";
    $cmd .= " samtools view -bh -F 524 -f 2 $tmp_bam.tmp > $tmp_bam.tmp_pairedmapped &&";
    $cmd .= " samtools view -bh -F 260 -f 8 $tmp_bam.tmp > $tmp_bam.tmp_R1mapped &&";
    $cmd .= " samtools view -bh -F 264 -f 4 $tmp_bam.tmp > $tmp_bam.tmp_R2mapped &&";
    $cmd .= " samtools sort -@ $num_threads -m $mem_per_thread $tmp_bam.tmp_pairedmapped -o $tmp_bam.tmp_pairedmapped.sorted.bam &&";
    $cmd .= " samtools sort -@ $num_threads -m $mem_per_thread $tmp_bam.tmp_R1mapped -o $tmp_bam.tmp_R1mapped.sorted.bam &&";
    $cmd .= " samtools sort -@ $num_threads -m $mem_per_thread $tmp_bam.tmp_R2mapped -o $tmp_bam.tmp_R2mapped.sorted.bam &&";
    $cmd .= " mv $tmp_bam.tmp_pairedmapped.sorted.bam ".$dir."/".$name."_pairedmapped.bam"." &&";
    $cmd .= " mv $tmp_bam.tmp_R1mapped.sorted.bam ".$dir."/".$name."_R1mapped.bam &&";
    $cmd .= " mv $tmp_bam.tmp_R2mapped.sorted.bam ".$dir."/".$name."_R2mapped.bam &&";
    $cmd .= " rm $tmp_bam.tmp && ";
    $cmd .= " rm $tmp_bam.tmp_pairedmapped &&";
    $cmd .= " rm $tmp_bam.tmp_R1mapped &&";
    $cmd .= " rm $tmp_bam.tmp_R2mapped";

    push(@cmds, $cmd);
    push(@bams_pairedmapped, $dir."/".$name."_pairedmapped.bam");
    push(@bams_R1mapped, $dir."/".$name."_R1mapped.bam");
    push(@bams_R2mapped, $dir."/".$name."_R2mapped.bam");
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
my $cmd_merge = "samtools merge -f ".$outfile_paired." ".join(" ", @bams_pairedmapped);
print STDERR "Attemtping to execute following cmd: $cmd_merge\n";
system($cmd_merge);
die "command failed: $!\n" if($? != 0);

$cmd_merge = "";
$cmd_merge = "samtools merge -f ".$outfile_R1." ".join(" ", @bams_R1mapped);
print STDERR "Attemtping to execute following cmd: $cmd_merge\n";
system($cmd_merge);
die "command failed: $!\n" if($? != 0);

$cmd_merge = "";
$cmd_merge = "samtools merge -f ".$outfile_R2." ".join(" ", @bams_R2mapped);
print STDERR "Attemtping to execute following cmd: $cmd_merge\n";
system($cmd_merge);
die "command failed: $!\n" if($? != 0);

#Then delete all bams
my $cmd_del = "rm ".join(" ", @bams_pairedmapped)." ".join(" ", @bams_R1mapped)." ".join("", @bams_R2mapped);
print STDERR "Attemtping to execute following cmd: $cmd_del\n";
system($cmd_del);
die "command failed: $!\n" if($? != 0);



exit;
