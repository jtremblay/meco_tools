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
--fasta <string>          : fasta file

OUTPUT:
--outdir <string>         : outdir

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@nrc-cnrc.gc.ca
ENDHERE

## OPTIONS
my ($help, $genomes_list, $fasta, $outdir);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'fasta=s'        => \$fasta,
  'genomes_list=i' => \$genomes_list,
  'outdir=s'       => \$outdir,
  'verbose'        => \$verbose,
  'help'           => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fasta fasta required\n") unless $fasta;
die("--outdir outdir required\n") unless $outdir;
#die("--genomes_list required\n") unless $genomes_list;

## MAIN
#my $i = 1;
#my $j = 1;
$outdir = abs_path($outdir);
#my %hash;

#open(IN, ">".$infile) or die "Can't open $infile\n";  
#while(<IN>){
#    chomp;
#    my @row = split(/\n/, $_);
#    $hash{$_} = $_;
#}
#close(IN);

my $db = new Iterator::FastaDb($fasta);
while(my $seq=$db->next_seq) {
    my $header = $seq->header;
    if($header =~ m/^>(\d+)_contig*/){
        my $genome_id = $1;
        open(OUT, ">>".$outdir."/".$genome_id.".fna") or die "Can't open ".$outdir."/".$genome_id.".fna\n";
        print OUT $seq->output;
        close(OUT); 
    }else{
        die("Something went wrong in the parsing of sequence headers...\n")
    }
}
exit;
