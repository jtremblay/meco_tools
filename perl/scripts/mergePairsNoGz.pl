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
mergePairsNoGz.pl

PURPOSE:

INPUT:
--reads1 <string>     : Sequence
--reads2 <string>     : Sequence
				
OUTPUT:
STDOUT                : stdout Interleaved fastq

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

# MAIN
my $in_1 = new Iterator::FastqDb($reads1,{trimN=>0, trim3=>0});
my $in_2 = new Iterator::FastqDb($reads2,{trimN=>0, trim3=>0});
while( my $curr_1 = $in_1->next_seq() ){
  my $curr_2 = $in_2->next_seq();
  print STDOUT "@".$curr_1->base()."#".$curr_1->barcode()."/".$curr_1->pair()."\n".$curr_1->seq()."\n+\n".$curr_1->qual."\n";
  print STDOUT "@".$curr_2->base()."#".$curr_2->barcode()."/".$curr_2->pair()."\n".$curr_2->seq()."\n+\n".$curr_2->qual."\n";
  
  die "Problem in file seq order...\n" if($curr_1->base() ne $curr_2->base());
}
