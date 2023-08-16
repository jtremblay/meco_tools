#!/usr/bin/env perl

use warnings;
use strict;

use Env qw/PATH/;
use Getopt::Long;
use File::Copy;
use Cwd;
use Iterator::FastaDb;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use File::Which;
use Data::Dumper;

$| = 1;    # UNBUFFER OUTPUT

my $usage = <<'ENDHERE';
NAME:
clustering4.pl

PURPOSE:
To count occurrences of all sequences. Intended for amplicon libraries. 
This script will cluster 16S rRNA sequences in a fashion similar to R. 
Edgar's OTU pipe pipeline, but using dnaclust (i.e. open source) 
instead of USEARCH. Briefly, there is: 1) 100% dereplication. 2) 99% 
clustering. 3)removal of clusters < 3. 4) removal of chimeras. 5) 
Clustering at 97% id.

INPUT:
--infile_fastq <fastq>     : Sequence file in fastq format.
--infile_fasta <fasta>     : Sequence file in fasta format.
--barcodes <string>        : Barcodes file in fasta format. Optional: used
                             for puting names in final cluster tables
                             instead of barcodes sequences.
--num_threads              : Number of threads. Only used in the 100%
                             dereplication step.
--clustering_perc_id <int> : clustering percentage similarity

OUTPUT:
--outdir <string>          : Output directory where files will be generated.

OPTIONS:
--verbose                  : print status messages to stdout

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

NOTE:

ENDHERE

my ($help, $infile_fastq, $infile_fasta, $start_at, $verbose, $outdir, $barcodes, $num_threads, $clustering_perc_id);
GetOptions(
  'infile_fasta=s'  => \$infile_fasta,
  'infile_fastq=s'  => \$infile_fastq,
  'barcodes=s'      => \$barcodes,
  'outdir=s'        => \$outdir,
  'clustering_perc_id=f' => \$clustering_perc_id,
  'num_threads=i'   => \$num_threads,
  'verbose'         => \$verbose,
  'start_at=i'      => \$start_at,
  'help'            => \$help
);

if ($help){ print $usage; exit;}

## Validate
die "--infile_fastq or --infile_fasta arg missing.\n" if(!defined($infile_fastq) and !defined($infile_fasta));
die "--outdir arg missing.\n" unless($outdir);
$num_threads = 1 unless($num_threads);

## Check paths.
my $cdhit = which('cd-hit-est'); chomp $cdhit;
print "path of cdhit:\t".$cdhit."\n" if($verbose);
die "Can't find cd-hit-est on path\n" if(!defined($cdhit));

$start_at = 0 unless($start_at);

##############
## BEG SUBS ##
##############

## Write clusters table.
## Input: clusters fasta sequence rep.
## Output: cluster table outfile and cluster fasta file.
sub writeClusters{
  
   my ($clstr, $obsTsv) = @_;

   my %HoH = ();
   my $i = 0;
   my %barcode_hash = ();
  
   ## Loop to get all barcodes count.
   open(IN, "<".$clstr) or die "Can't open ".$clstr."\n";
   my $seed;
   while(<IN>){
      chomp;
      if($_ =~  m/^>Cluster (\d+)/){
         $seed = "cluster_$1";
         $HoH{$seed}{size} = 0;
      }elsif($_ =~ m/>gene_id_\d+#([ACGTN]+)\/1/){
         if(exists $HoH{$seed}{$1}){
            $HoH{$seed}{$1}++;
         }else{
            $HoH{$seed}{$1} = 1;
         }
         $HoH{$seed}{size}++;
      }
      $i++;
      #last if($i > 100);
   }
  # Open outfile for writing.
  open(OUT, ">".$obsTsv) or die "Can't open ".$obsTsv."\n";
  
  # Loop through barcodes
  my %names = ();
  if($barcodes){
    my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file ".$barcodes."\n");
    while( my $curr = $ref_fasta_db->next_seq() ) {
      my $header = $curr->header;
      $header =~ s/>//;
      $header =~ s/\t//g;
      $names{$curr->seq} = $header;
    }    
  }
  
  print OUT "#CLUSTER";
  foreach my $ref_barcode (sort keys %names) {
     print OUT "\t".$names{$ref_barcode};
  }
  print OUT "\n";
  
  # Then print table.
  foreach my $seqid (sort {$HoH{$b}{size} <=> $HoH{$a}{size} } keys %HoH){
    print OUT $seqid;  
  
    # Print barcodesin order.
    foreach my $ref_barcode (sort keys %names) {
      if(exists $HoH{$seqid}{$ref_barcode} ){
        print OUT "\t".$HoH{$seqid}{$ref_barcode};
      }else{
        print OUT "\t0";
      }
    }
    print OUT "\n";        
  }
  close(OUT);
}
##############
## END SUBS ##
##############

## MAIN

## Cluster at 99%
my $cmd = $cdhit;
$cmd .= " -T ".$num_threads;
$cmd .= " -c 0.95";
$cmd .= " -d 100";
$cmd .= " -n 11";
$cmd .= " -M 0";
$cmd .= " -i ".$infile_fasta; 
$cmd .= " -o ".$outdir."/cdhit_".$clustering_perc_id.".fasta";
system($cmd) if($start_at <= 1);
die "command failed: $!\n" if($? != 0);

## Generate Obs table
writeClusters($outdir."/cdhit_".$clustering_perc_id.".fasta.clstr", $outdir."/obs.tsv") if($start_at >= 2);

exit;
