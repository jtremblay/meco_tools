#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
    convert_fastq_to_fasta.pl
PURPOSE:
The name says it all.
INPUT:
--fastq <fasta_infile>

OUTPUT:
STDOUT <fasta_outfile> 

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $fastq, $fasta);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'fastq=s'   => \$fastq,
  'fasta=s'   => \$fasta,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fastq fastq file required\n") unless $fastq;
#die("--fasta fasta required\n") unless $fasta;

## MAIN
#open(OUT, ">".$fasta) or die "Can't open file ".$!."\n";

my $db = new Iterator::FastqDb($fastq);
while(my $seq=$db->next_seq) {
  my $header = $seq->header;
  #$header =~ s/#/-/g;
  $header =~ s/@// if(substr($header,0,1) eq "@");
  print STDOUT ">".$header."\n".$seq->seq()."\n";
}

#close(OUT);


exit;
