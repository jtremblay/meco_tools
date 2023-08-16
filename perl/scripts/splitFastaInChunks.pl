#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
splitFastaInChunks.pl
PURPOSE:

INPUT:
--fasta <string>  : fasta file
--n <int>         : how many sequences per file.

OUTPUT:
--outdir <string> : outdir

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $n, $fasta, $outdir);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'fasta=s'   => \$fasta,
  'n=i'       => \$n,
  'outdir=s'  => \$outdir,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fasta fasta required\n") unless $fasta;
die("--outdir outdir required\n") unless $outdir;
die("--n required\n") unless $n;

## MAIN
my $i = 1;
my $j = 1;
$outdir = abs_path($outdir);
my $db = new Iterator::FastaDb($fasta);
   
my $outfile = $outdir."/nr_".$j.".fasta";
open(OUT, ">".$outfile) or die "Can't open $outfile\n";
while(my $seq=$db->next_seq) {
   
    if($i % $n == 0){
      close(OUT);  
      $j++;
      $outfile = $outdir."/nr_".$j.".fasta";
      open(OUT, ">".$outfile) or die "Can't open $outfile\n";
      print OUT $seq->output;
      print STDERR "j: ".$j."\n"; 
    }else{
      print OUT $seq->output; 
    }
    $i++;
}

close(OUT);

exit;
