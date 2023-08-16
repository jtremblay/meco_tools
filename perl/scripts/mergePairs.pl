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

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
mergePairs.pl

PURPOSE:

INPUT:
--reads1 <string>     : Sequence file gz
--reads2 <string>     : Sequence file gz
				
OUTPUT:
--reads1_out <string> : paired reads1 fastq.gz
--reads2_out <string> : paired reads2 fastq.gz
STDOUT                : stdout in gzip. Interleaved fastq.gz

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $reads1, $reads2, $reads1_out, $reads2_out);
my $verbose = 0;

GetOptions(
   'reads1=s' 	  => \$reads1,
   'reads2=s'     => \$reads2,
   'reads1_out=s' => \$reads1_out,
   'reads2_out=s' => \$reads2_out,
   'verbose' 	  => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

# VALIDATE
die "--reads1 missing\n" unless($reads1);
die "--reads2 missing\n" unless($reads2);
die "--reads1_out missing\n" unless($reads1_out);
die "--reads2_out missing\n" unless($reads2_out);

# MAIN

my $tmpdir = File::Temp->newdir(
    "tmpdir-mergePairs-XXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

my $pipe1 = "$tmpdir/reads1.pipe";
my $pipe2 = "$tmpdir/reads2.pipe";
system("mkfifo $pipe1");
system("mkfifo $pipe2");
system("gunzip -c ".$reads1." > $pipe1 &");
system("gunzip -c ".$reads2." > $pipe2 &");

# Loop through R1
my %hash1;
my %hash2;
my $ref_fastq_db = Iterator::FastqDb->new($pipe1) or die("Unable to open Fastq file, $pipe1\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
   $hash1{$curr->base} = "";
   print STDERR "R1 ".$curr->base."\n" if($verbose);
}
system("rm $pipe1");

# Loop through R2
undef $ref_fastq_db;
$ref_fastq_db = Iterator::FastqDb->new($pipe2) or die("Unable to open Fastq file, $pipe2\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
   if(exists $hash1{$curr->base}){
      $hash2{$curr->base} = "";
      print STDERR "R2 ".$curr->base."\n" if($verbose);
   }else{
      delete $hash1{$curr->base};
      #print STDERR "R2 ".$curr->base." does not exists... in reads 1\n";
   }
}

# close pipes
system("rm $pipe2");

# open new pipes to read R1 and R2 and write them in a single interleaved file and R1 and R2 paired.
# could probably use same pipes as above, but create new ones just to make sure.
# Declare output;
my $output = "-";
my $z = new IO::Compress::Gzip $output or die "gzip failed: $GzipError\n";
my $z1 = new IO::Compress::Gzip $reads1_out or die "gzip failed: $GzipError\n";
my $z2 = new IO::Compress::Gzip $reads2_out or die "gzip failed: $GzipError\n";
my $pipe1b = "$tmpdir/reads1b.pipe";
my $pipe2b = "$tmpdir/reads2b.pipe";
system("mkfifo $pipe1b");
system("mkfifo $pipe2b");
system("gunzip -c ".$reads1." > $pipe1b &");
system("gunzip -c ".$reads2." > $pipe2b &");

# Dump hash for debugging
#print STDERR Dumper(\%hash1);
#print STDERR Dumper(\%hash2);

# Print paired 1 and 2 separately.
my $ref_fastq_db1 = Iterator::FastqDb->new($pipe1b) or die("Unable to open Fastq file, $pipe1b\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
   
    if( exists $hash1{$curr1->base} && exists $hash2{$curr1->base} ){
      $z1->print($curr1->output);
   }
}
$z1->close();
system("rm $pipe1b");

my $ref_fastq_db2 = Iterator::FastqDb->new($pipe2b) or die("Unable to open Fastq file, $pipe2b\n");
while( my $curr2 = $ref_fastq_db2->next_seq() ) {
   
    if( exists $hash1{$curr2->base} && exists $hash2{$curr2->base} ){
      $z2->print($curr2->output);
   }
}
$z2->close();
system("rm $pipe2b");

# Then from created files write interleaved file.
my $pipe1c = "$tmpdir/reads1c.pipe";
my $pipe2c = "$tmpdir/reads2c.pipe";
system("mkfifo $pipe1c");
system("mkfifo $pipe2c");
system("gunzip -c ".$reads1_out." > $pipe1c &");
system("gunzip -c ".$reads2_out." > $pipe2c &");
undef $ref_fastq_db1;
undef $ref_fastq_db2;
$ref_fastq_db1 = Iterator::FastqDb->new($pipe1c) or die("Unable to open Fastq file, $pipe1c\n");
$ref_fastq_db2 = Iterator::FastqDb->new($pipe2c) or die("Unable to open Fastq file, $pipe2c\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
   my $curr2 = $ref_fastq_db2->next_seq();
   $z->print($curr1->output."".$curr2->output);
}
$z->close();
system("rm $pipe1c");
system("rm $pipe2c");

## REMOVE TEMP FILES
sub END{
  local $?;
  system("rm ".$tmpdir." -rf");
}
