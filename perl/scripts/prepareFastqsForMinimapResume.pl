#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
prepareFastqsForMinimapResume.pl

PURPOSE:

INPUT:
--infile_paf <string>    : Sequence file
--infile_fastq <string>  : Fastq file
				
OUTPUT:
STDOUT : Fastq file to use as query sequence file in next run of minimap 

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_paf, $infile_fastq);
my $verbose = 0;

GetOptions(
   'infile_paf=s'   => \$infile_paf,
   'infile_fastq=s' => \$infile_fastq,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-prepareFastqsPaf-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);
print STDERR "$tmpdir\n";

# MAIN
my %hash;
my $last_query_seq = ""; # important to remove from hash as the alignment process of the last query seq is most probably incomplete
open(IN, "zcat $infile_paf |") or die "zcat $infile_paf: $!";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    $hash{$row[0]} = "";
    $last_query_seq = $row[0];
}
close(IN);
delete($hash{$last_query_seq});

## MAIN

my $pipe = "$tmpdir/reads.pipe";
system("mkfifo $pipe");
system("gunzip -c ".$infile_fastq." > $pipe &");

print STDERR Dumper(\%hash);

my $total_seqs = 0;
my $total_aligned_seqs = 0;
my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fastq file, $pipe\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
    $total_seqs++;

    #print STDERR "[DEBUG]: ".$curr->header."\n";
    my $header = $curr->header;
    $header =~ s/^@//;

    if(exists $hash{$header}){
        $total_aligned_seqs++;
    }else{
        print STDOUT $curr->output;
    }
}
print STDERR "Aligned seqs: $total_aligned_seqs\n";
print STDERR "Total seqs: $total_seqs\n";
exit;
