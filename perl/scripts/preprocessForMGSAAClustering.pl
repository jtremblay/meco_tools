#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw/ TMPDIR SCRATCH /;
use File::Temp;
use Cwd 'abs_path';

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
preprocessForMGSAAClustering.pl

PURPOSE:
Takes an indirectory containing .fastq.gz. and generates one interleaved fastq file for all .fastq.gz
files in your directory. This is intended for fastq libs that have already be multiplexed but for which
barcode information/sequence is missing. The resulting files (barcodes and fastq) can then be used
as input for bioniformatics pipelines (rRNA tags, Ray, etc.).

INPUT:
--infiles <string>          : List of infiles in fastq.gz separated by a comma
				
OUTPUT:
--outfile_barcodes <string> : barcodes sequence file. 
--outfile_fastq <string>    : outfile where fastq (gz format) are written.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $outfile_barcodes, $outfile_fasta);
my $verbose = 0;

GetOptions(
   'infiles=s' 	      => \$infiles,
   'outfile_barcodes=s' => \$outfile_barcodes,
   'outfile_fasta=s'    => \$outfile_fasta,
   'verbose' 	         => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $tmpdir = File::Temp->newdir(
    "tmpdir-assignBarcodes-XXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

my %hash;
my @infiles = split(/,/, $infiles);
foreach(@infiles){
   chomp;
   $hash{$_}{R1} = $_;
}

# Get barcodes
my $i=0;
my $j=0;
my @barcodes;
for(glob '{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}'){
	if($i % 30 == 0){ 
		$barcodes[$j] = $_; 
		$j++;
	}else{

	}   
	$i++;
} 

# Declare output;
#my $z = new IO::Compress::Gzip $outfile_fastq or die "gzip failed: $GzipError\n";

open(OUTBAR, ">".$outfile_barcodes) or die "Can't open $outfile_barcodes\n";
open(OUT, ">".$outfile_fasta) or die "Can't open $outfile_fasta\n";

my $x = 0;
for my $key (keys %hash){

   # define barcode sequence
   my $currBarcode = $barcodes[$x];

   # print to barcodes file
   print OUTBAR ">".$key."\n".$currBarcode."\n";

   my $ref_fasta_db1 = Iterator::FastaDb->new($hash{$key}{R1}) or die("Unable to open Fastq file, ".$hash{$key}{R1}."\n");
   while( my $curr1 = $ref_fasta_db1->next_seq() ) {
       #$z->print(">".$curr1->header."#".$currBarcode."/1\n".$curr1->seq."\n");
       print OUT $curr1->header."#".$currBarcode."\n".$curr1->seq."\n";
   }
   $x++;
}
close(OUTBAR);
close(OUT);
#$z->close();


## REMOVE TEMP FILES
sub END{
  local $?;
  system("rm ".$tmpdir." -rf");
}
