#!/usr/bin/env perl

use strict;
use warnings;

#use Env qw/TMPDIR/;
use Bio::Seq;
use Bio::Tools::IUPAC;
#use String::Approx 'aslice';
#use List::Util qw(sum);
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
generateAllPossibleCombinations.pl

PURPOSE:
Take fasta with one entry and generate multi fasta file with all possible combinations of sequences.

INPUT:
--infile <string>          : Sequence file.

OUTPUT:
STDOUT
--comp <string> 
--rev <string> 
--revcomp <string>

OPTIONAL:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $insert_size, $rev, $revcomp, $comp);
my $verbose = 0;

GetOptions(
  'infile=s'                   => \$infile,
  'rev=s'                      => \$rev,
  'revcomp=s'                  => \$revcomp,
  'comp=s'                     => \$comp,
  'insert_size=i'              => \$insert_size,
  'verbose'                    => \$verbose,
  'help'                       => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);

## GENERATE PRIMER SEQS
my @left_primer;
my @right_primer;
my $ref_fasta_db = new Iterator::FastaDb($infile) or die("Unable to open Fasta file, $infile\n");
my $i=1;
while( my $ref_seq = $ref_fasta_db->next_seq() ) {
    my $primer_seq = uc($ref_seq->seq());
    if($i == 1){
        @left_primer = generate_primer_combination($primer_seq);
    }
    #elsif($i == 2){
    #    @right_primer = generate_primer_combination($primer_seq);
    #}elsif($i > 2){
    #    die "Please provide only one sequence in fasta file\n";
    #}
    $i++;
}
if($rev){
    open(OUT_REV, ">",$rev) or die "Can't open $rev\n";
}
if($revcomp){
    open(OUT_REVCOMP, ">",$revcomp) or die "Can't open $revcomp\n";
}
if($comp){
    open(OUT_COMP, ">",$comp) or die "Can't open $comp\n";
}

foreach my $left (@left_primer){
    #foreach my $right (@right_primer){
        print STDOUT $left."\n";#"\t".$right."\t".$insert_size."\n";
        if($rev){
            print OUT_REV reverse($left)."\n";
        }
        if($revcomp){
            print OUT_REVCOMP reverse_complement_IUPAC($left)."\n";
        }
        if($comp){
            print OUT_COMP complement_IUPAC($left)."\n";
        }
    #}
}


## GENERATE ALL POSSIBLE PRIMER COMBINATION
## input: one dna sequence that can contain ambiguous bases.
## output: an array of dna sequences with no ambiguous bases.
sub generate_primer_combination{
    my $primer_seq = $_[0];
    my @array = ();

    #For some obscur reason, Bio::Seq would always accept sequence as provided...
    $primer_seq =~ m/(\S+)/;
    $primer_seq = $1;

    my $ambiseq = Bio::Seq->new(-seq => $primer_seq, -alphabet => 'dna');
    my $stream  = Bio::Tools::IUPAC->new(-seq => $ambiseq);

    while (my $uniqueseq = $stream->next_seq()) {
        push(@array, $uniqueseq->seq());
    }
    return @array;
}


sub reverse_complement_IUPAC {
   my $dna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $revcomp;
}

sub complement_IUPAC {
   my $dna = shift;


   # complement the reversed DNA sequence
   my $comp = $dna;
   $comp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $comp;
}

exit;
