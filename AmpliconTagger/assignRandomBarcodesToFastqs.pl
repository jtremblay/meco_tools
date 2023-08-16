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
use Env qw(TMPDIR);
use File::Temp;
use Cwd 'abs_path';
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
assignRandomBarcodesToFastqs.pl

PURPOSE:
Takes an indirectory containing .fastq.gz. and generates one interleaved fastq file for all .fastq.gz
files in your directory. This is intended for fastq libs that have already be multiplexed but for which
barcode information/sequence is missing. The resulting files (barcodes and fastq) can then be used
as input for bioniformatics pipelines (rRNA tags, metagenome assembly, etc.).

INPUT:
--indir <string>            : Directory in which are your sequence files
--se                        : Set this flag if you have ion torrent data instead of miseq. 
                              Will assume reads are single ended.
--external                  : Set this flag to avoid searching for S\d+_L\d+
--flexible                  : Set this flag if samples have various labelling nomenclature.

OUTPUT:
--outfile_barcodes <string> : barcodes sequence file. 
--outfile_fastq <string>    : outfile where fastq (gz format) are written.
--min_len <int>             : minimum length filter for a read to be accepted.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $outfile_barcodes, $outfile_fastq, $min_len, $se, $external, $flexible);
my $verbose = 0;

GetOptions(
   'indir=s'             => \$indir,
   'outfile_barcodes=s'  => \$outfile_barcodes,
   'outfile_fastq=s'     => \$outfile_fastq,
   'min_len=i'           => \$min_len,
   'se'                  => \$se,
   'external'            => \$external,
   'flexible'            => \$flexible,
   'verbose'             => \$verbose,
   'help'                => \$help
);
if ($help) { print $usage; exit; }
my %hash;

die "--min_len <int> missing...\n" unless($min_len);

# SUBS
sub eachFile{
    my $filename = $_;
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_

    if(-e $filename) { 
        if($se){
            if(substr($filename, -9) eq ".fastq.gz"){
                print STDERR "Found $fullpath...\n";
                if($filename =~ m/(.*)\.fastq\.gz/){
                    my $prefix = $1;
                    my $R = "R1";
                    $hash{$prefix}{$R} = $fullpath;
                }else{
                    print STDERR "$filename do not match expected nomenclature and will be left out.\n";
                }
            }
        }elsif($external){
            if(substr($filename, -9) eq ".fastq.gz"){
                print STDERR "Found $fullpath...\n";
                if($filename =~ m/(.*)_(R\d)\.fastq.gz/){
                    my $prefix = $1;
                    my $R = $2;
                    $hash{$prefix}{$R} = $fullpath;
                }else{
                    print STDERR "$filename do not match expected nomenclature and will be left out.\n";
                }
            }
        }elsif($flexible){
            if(substr($filename, -9) eq ".fastq.gz"){
                print STDERR "Found $fullpath...\n";
                if($filename =~ m/(.*)_(R\d).*\.fastq.gz/){
                    my $prefix = $1;
                    my $R = $2;
                    $hash{$prefix}{$R} = $fullpath;
                }else{
                    print STDERR "$filename do not match expected nomenclature and will be left out.\n";
                }
            }
        }else{
            if(substr($filename, -9) eq ".fastq.gz"){
                print STDERR "Found $fullpath...\n";
                if($filename =~ m/(.*)_S\d+_L\d+_(R\d)_.*/){
                    my $prefix = $1;
                    my $R = $2;
                    $hash{$prefix}{$R} = $fullpath;
                }
                else{
                    print STDERR "$filename do not match expected nomenclature and will be left out.\n";
                }
            }
        }
    }
}


## MAIN

my $tmpdir = File::Temp->newdir(
    "tmpdir-assignBarcodes-XXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

## REMOVE TEMP FILES
sub END{
  local $?;
  system("rm ".$tmpdir." -rf");
}

# Populate hash
$indir = abs_path($indir)."/";

find (\&eachFile, $indir);

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

print STDERR Dumper(\%hash);

for my $key (keys %hash){
    my $num = 0;
    foreach my $key2 (keys %{ $hash{$key} }) {
        $num++;
    }
    if($num == 1){
        if($se){
            print STDERR $key." is single-end fastq.\n";
        }else{
            print STDERR $key." is single-end fastq and will be left out.\n";
            delete $hash{$key};
        }
    }
}


# Declare output;
my $z = new IO::Compress::Gzip $outfile_fastq or die "gzip failed: $GzipError\n";

open(OUTBAR, ">".$outfile_barcodes) or die "Can't open $outfile_barcodes\n";

my $x = 0;
for my $key (keys %hash){
   # make pipe
   my $pipe1 = "$tmpdir/reads1.pipe";
   my $pipe2 = "$tmpdir/reads2.pipe" if(!$se);
   system("mkfifo $pipe1");
   system("mkfifo $pipe2") if(!$se);
   # gunzip to pipes
   my $gz1 = $hash{$key}{R1};
   my $gz2 = $hash{$key}{R2} if(!$se);
   print STDERR "[DEBUG] Integrating $gz1\n";
   system("gunzip -c ".$gz1." > $pipe1 &");
   print STDERR "[DEBUG] Integrating $gz2\n";
   system("gunzip -c ".$gz2." > $pipe2 &") if(!$se);

   # define barcode sequence and remove all potentially problematic characters.
   my $currBarcode = $barcodes[$x];
   my $currKey = $key;
   #$currKey =~ s/_/\./;
   #$currKey =~ s/-/\./;
   $currKey =~ s/\//\./;
   $currKey =~ s/\\/\./;
   $currKey =~ s/\(/\./;
   $currKey =~ s/\)/\./;
   $currKey =~ s/\s+/\./;

   # Print to barcodes file
   print OUTBAR ">".$currKey."\n".$currBarcode."\n";

   my $ref_fastq_db1 = Iterator::FastqDb->new($pipe1) or die("Unable to open Fastq file, $pipe1\n");
   my $ref_fastq_db2 = Iterator::FastqDb->new($pipe2) or die("Unable to open Fastq file, $pipe2\n") if(!$se);
   while( my $curr1 = $ref_fastq_db1->next_seq() ) {
      my $curr2 = $ref_fastq_db2->next_seq() if(!$se);
      
      my $length_R1 = length($curr1->seq);
      my $length_R2;
      if(!$se){
         $length_R2 = length($curr2->seq) 
      }else{
         $length_R2 = $min_len + 1;
      }

      if($length_R1 < $min_len || $length_R2 < $min_len){
          #print to failed file archive?
      }else{
         $z->print("@".$curr1->base."#".$currBarcode."/1\n".$curr1->seq."\n+\n".$curr1->qual."\n");
         $z->print("@".$curr2->base."#".$currBarcode."/2\n".$curr2->seq."\n+\n".$curr2->qual."\n") if(!$se);
      }
   }
   system("rm $pipe1");
   system("rm $pipe2") if(!$se);
   $x++;
}
close(OUTBAR);
$z->close();

