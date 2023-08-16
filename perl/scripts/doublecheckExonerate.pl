#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
doublecheckExonerate.pl

PURPOSE:
Sometimes, exonerate does not output all files. This script will create (touch) empty files.

INPUT:
--indir <string>  : Directory where exonerate output is
--no_chunks <int> : number of chunks
--prefix <string> : prefix to files

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $no_chunks, $prefix);
my $verbose = 0;

GetOptions(
   'indir=s' 	  => \$indir,
   'no_chunks=i' => \$no_chunks,
   'prefix=s'    => \$prefix,
   'verbose' 	  => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN

$indir = abs_path($indir);

for(my $i=0; $i<$no_chunks; $i++){
   if(-e $indir."/".$prefix.sprintf("%07d", $i)){
      #do nothing..
   }else{
       system("touch ".$indir."/".$prefix.sprintf("%07d", $i));
       print STDERR ("touch ".$indir."/".$prefix.sprintf("%07d", $i))." not found... empty file created...\n";
   }
}
exit;
