#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
contigsToCanopy.pl

PURPOSE:

INPUT:
--infile_canopy <string> : output of canopy algorithm.
--contigs <string>       : string of contigs
--n <int>                : number of canopies to process. Default=1000

OUTPUT:
--outdir <string>        : outdir

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_canopy, $contigs, $outdir, $n);
my $verbose = 0;

GetOptions(
   'infile_canopy=s' => \$infile_canopy,
   'contigs=s' 	   => \$contigs,
   'outdir=s' 	      => \$outdir,
   'n=i'             => \$n,
   'verbose' 	      => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

$n = 1000 unless($n);

## MAIN
my $outfile = $outdir."/canopy_contigs.fasta";
open(OUT, ">".$outfile) or die "Can't open $outfile\n";

my %hash;
my %seen;
my %hashOutfiles;
my $counter = 0;
open(IN, "<".$infile_canopy) or die "Can't open $infile_canopy\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $canopy = $row[0];
   my $gene = $row[1];
   my @gene = split(/=/, $gene);
   my $contig = $gene[0];
   $contig =~ s/\./ /g;

   $hash{$contig} = $canopy;
   $seen{$canopy} = 1;

   # Open files. One file per canopy
   if(exists $seen{$canopy}){
   
   }else{
      #open (my $fh, ">", "$outdir/$canopy.fasta") or die "Can't open $outdir/$canopy.fasta\n";
      #$hashOutfiles{$canopy} = $fh;
      #open $hashOutfiles{$canopy}, ">", "$outdir/$canopy.fasta" or die "Can't open $outdir/$canopy.fasta\n";
      $counter++;
   }
   last if($counter == $n);
}
close(IN);
#print STDERR Dumper(\%hash);

my @contigs = split(/,/, $contigs);
foreach my $currContigs (@contigs){

   print STDERR "Processing contig: $currContigs\n";

   my $ref_fasta_db = Iterator::FastaDb->new($currContigs) or die("Unable to open Fasta file, $currContigs\n"); 
   while( my $curr = $ref_fasta_db->next_seq() ) {
      my $seq = $curr->seq;
      my $header = $curr->header;
      $header =~ s/^>//;
         
      if(exists $hash{$header}){
         my $key = $hash{$header};
         #print {$hashOutfiles{$key}} ">".$header."\n".$seq."\n";
         print OUT ">".$hash{$header}."=".$header."\n".$seq."\n";
      }else{
       
      }
   }
}
close(OUT);
# Close all file handles
foreach my $key (keys %hashOutfiles){
   close $hashOutfiles{$key};
}

   
my $cmd = "bwa index $outfile";
system($cmd);
$? != 0 ? die "command failed: $!\n" : print STDERR "Succesfuly made index...\n";

exit;

