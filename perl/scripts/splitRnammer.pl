#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
splitRnammer.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file (output of rnammmer).
				
OUTPUT:
--S5 <string>   : ssu 5S
--S8 <string>   : ssu 5S
--S16 <string>  : ssu 16S
--S18 <string>  : ssu 18S
--S23 <string>  : lsu 23S
--S28 <string>  : lsu 28S

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $S5, $S8, $S16, $S18, $S23, $S28);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'S5=s'      => \$S5,
   'S8=s'      => \$S8,
   'S16=s'     => \$S16,
   'S18=s'     => \$S18,
   'S23=s'     => \$S23,
   'S28=s'     => \$S28,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT5, ">".$S5) or die "Can't open $S5\n";
open(OUT8, ">".$S8) or die "Can't open $S8\n";
open(OUT16, ">".$S16) or die "Can't open $S16\n";
open(OUT18, ">".$S18) or die "Can't open $S18\n";
open(OUT23, ">".$S23) or die "Can't open $S23\n";
open(OUT28, ">".$S28) or die "Can't open $S28\n";

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   if($curr->header =~ m/molecule=5s_rRNA/){
      print OUT5 $curr->output;
   }elsif($curr->header =~ m/molecule=8s_rRNA/){
      print OUT8 $curr->output;
   }elsif($curr->header =~ m/molecule=16s_rRNA/){
      print OUT16 $curr->output;
   }elsif($curr->header =~ m/molecule=18s_rRNA/){
      print OUT18 $curr->output;
   }elsif($curr->header =~ m/molecule=23s_rRNA/){
      print OUT23 $curr->output;
   }elsif($curr->header =~ m/molecule=28s_rRNA/){
      print OUT28 $curr->output;
   }
}
close(OUT5);
close(OUT8);
close(OUT16);
close(OUT18);
close(OUT23);
close(OUT28);

exit;
