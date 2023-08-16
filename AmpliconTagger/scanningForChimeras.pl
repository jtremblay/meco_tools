#!/usr/bin/env perl

use warnings;
use strict;

use Env qw/PATH/;
use Getopt::Long;
use File::Copy;
use Cwd;
use Iterator::FastaDb;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use File::Which;

$| = 1;    # UNBUFFER OUTPUT

my $usage = <<'ENDHERE';
NAME:
scanningForChimeras.pl

PURPOSE:
Scan for chimeras

INPUT:
--infile_fasta <string>  : Sequence file in fasta format.
--infile_tsv <string>    : Feature/SNV abundance matrix.
--ref_db <string>        : Reference DB for chimera UCHIME reference
--num_threads <int>      : Number of threads. Default = 1.
--ref_only               : Set this flag if Uchime-reference only is to be performed

OUTPUT:
--outdir <string>        : Outfile string.

OPTIONS:
--verbose                : print status messages to stdout

AUTHOR/SUPPORT:
Genomics and Microbiomes

Julien Tremblay - jtremblay514@gmail.com

NOTE:

ENDHERE

my ($help, $infile_tsv, $infile_fasta, $ref_db, $start_at, $outdir, $num_threads, $ref_only);
GetOptions(
  'infile_fasta=s'  => \$infile_fasta,
  'infile_tsv=s'    => \$infile_tsv,
  'ref_db=s'        => \$ref_db,
  'outdir=s'        => \$outdir,
  'num_threads=i'   => \$num_threads,
  'ref_only'        => \$ref_only,
  'start_at=i'      => \$start_at,
  'help'            => \$help
);

if ($help){ print $usage; exit;}

## Validate
die "--infile_fasta arg missing.\n" if(!defined($infile_fasta));
die "--ref_db arg is missing (fasta reference database).\n" unless($ref_db);
die "--outdir arg missing.\n" unless($outdir);
$num_threads = 1 unless($num_threads);

## Check paths.
$start_at = 0 unless($start_at);

## Check paths.
my $vsearch= which('vsearch'); chomp $vsearch;
print "path of vsearch:\t".$vsearch."\n";
die "Can't find vsearch on path\n" if(!defined($vsearch));

$start_at = 0 unless($start_at);

#Chimera detection (UCHIME de novo mode):
#  usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
#     [-uchimeout results.uch] [-uchimealns results.alns]
#  Input is estimated amplicons with integer abundances specified using ";size=N".
#usearch -uchime_denovo

my $basename = $infile_fasta;
$basename =~ s{.*/}{};      # removes path  
$basename =~ s{\.[^.]+$}{}; # removes extension

if($ref_only){
    print STDERR "Performing uchime-reference only...\n";
    my $cmd = $vsearch;
    $cmd .= " --uchime_ref ".$infile_fasta;
    $cmd .= " --db ".$ref_db;
    $cmd .= " --nonchimeras ".$outdir."/".$basename."_nonChimRef.fasta";
    $cmd .= " --threads ".$num_threads;
    print STDERR $cmd."\n";
    system($cmd) if($start_at <= 0);
    die "command failed: $!\n" if($? != 0);

}else{

    my $cmd = $vsearch;
    $cmd .= " --uchime_denovo ".$infile_fasta;
    $cmd .= " --nonchimeras ".$outdir."/".$basename."_nonChimDeNovo.fasta";
    print STDERR $cmd."\n";
    system($cmd) if($start_at <= 0);
    die "command failed: $!\n" if($? != 0);
    
    #usearch -uchime_ref
    #Chimera detection (UCHIME ref. db. mode):
    #  usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
    #    [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]
    $cmd = $vsearch;
    $cmd .= " --uchime_ref ".$outdir."/".$basename."_nonChimDeNovo.fasta";
    $cmd .= " --db ".$ref_db;
    $cmd .= " --nonchimeras ".$outdir."/".$basename."_nonChimDeNovoRef.fasta";
    $cmd .= " --threads ".$num_threads;
    print STDERR $cmd."\n";
    system($cmd) if($start_at <= 1);
    die "command failed: $!\n" if($? != 0);
}

# Then re-generate obs table free of chimeras. And remove ;size=<int> pattern.
open(OUT, ">".$outdir."/obs.fasta") or die "Can't open ".$outdir."/obs.fasta"."\n";
my %hash;

my $ref_fasta;
if($ref_only){
    $ref_fasta = $outdir."/".$basename."_nonChimRef.fasta";
}else{
    $ref_fasta = $outdir."/".$basename."_nonChimDeNovoRef.fasta";
}

my $ref_fasta_db = Iterator::FastaDb->new($ref_fasta) or die("Unable to open Fasta file ".$ref_fasta."\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/>//;
    $header =~ s/;size=\d+//;
    $hash{$header} = $curr->seq;
    #print OUT ">".$header."\n".$curr->seq."\n";
}
for my $header ( sort {$a<=>$b} keys %hash) {
    print OUT ">".$header."\n".$hash{$header}."\n";
}
close(OUT);

open(OUT, ">".$outdir."/obs.tsv") or die "Can't open ".$outdir."/obs.tsv"."\n";
open(IN, "<".$infile_tsv) or die "Can't open $infile_tsv\n";
while(<IN>){
    chomp;
    if($. == 1){ print OUT $_."\n";next;}
    my @row = split(/\t/, $_);
    my $cluster_id = $row[0];
    if(exists $hash{$cluster_id}){
        print OUT $_."\n";
    }
}
close(IN);
close(OUT);

 exit;
