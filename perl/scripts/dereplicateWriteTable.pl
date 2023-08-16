#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use File::Temp;
use POSIX qw(mkfifo);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
dereplicateWriteTable.pl

PURPOSE:

INPUT:
--fasta <string>           : Sequence file (fasta format).
        
OUTPUT:
--outdir <string>          : Sequence file 100% clustered.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Microbiomes and Genomics
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $fasta, $outdir);
my $verbose = 0;

GetOptions(
  'fasta=s'              => \$fasta,
  'outdir=s'             => \$outdir,
  'verbose'              => \$verbose,
  'help'                 => \$help
);
if ($help) { print $usage; exit; }
 
## Validate
die "--fasta OR --fastq required.\n" if(!defined($fasta));

my $tmpdir = File::Temp->newdir(
    "tmpDirDereplicateGzXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0
);


####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################

## Write clusters table.
## Input: clusters fasta sequence rep.
## Output: cluster table outfile
sub writeClusters{
  
  my ($clusters_fasta, $clusters_table, $new_cluster_fasta) = @_;

  my %HoH = ();
  my $i = 0;
  my %barcode_hash = ();
  
  ## Loop to get all barcodes count.
  my $ref_fasta_db = Iterator::FastaDb->new($clusters_fasta) or die("Unable to open Fasta file ".$clusters_fasta);
  while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    my @row = split(/;/, $header);
    
    my $seqid = shift(@row);
    $seqid =~ s/>//;
    my $abundance = pop(@row);
  
    $abundance =~ m/(\d+)/;
    $HoH{$seqid}{abundance} = $1;
    
    foreach(@row){
      $_ =~ m/\#(\S+)=(\d+)/;
      $HoH{$seqid}{$1} = $2;
      $barcode_hash{$1} = 1;
    }
  }
  
  # Open outfile for writing.
  open(OUT, ">".$clusters_table) or die "Can't open ".$clusters_table."\n";
  
  # Print header
  my %names = ();
  #if($barcodes){
  #  $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file ".$barcodes."\n");
  #  while( my $curr = $ref_fasta_db->next_seq() ) {
  #    my $header = $curr->header;
  #    $header =~ s/>//;
  #    $names{$curr->seq} = $header;
  #  }    
  #}
  
  print OUT "#CLUSTER";
  foreach my $ref_barcode (sort keys %barcode_hash) {
      #if($barcodes){
      #print OUT "\t".$names{$ref_barcode};
      #}else{
      $ref_barcode =~ s/_unmapped\.fastq\.gz//;
      print OUT "\t".$ref_barcode;
      #}
  }
  print OUT "\n";
  
  # Then print table.
  foreach my $seqid (sort {$HoH{$b}{abundance} <=> $HoH{$a}{abundance} } keys %HoH){
    print OUT $seqid;  
  
    # Print barcodesin order.
    foreach my $ref_barcode (sort keys %barcode_hash) {
      if(exists $HoH{$seqid}{$ref_barcode} ){
        print OUT "\t".$HoH{$seqid}{$ref_barcode};
      }else{
        print OUT "\t0";
      }
    }
    print OUT "\n";        
  }
  close(OUT);

  # Print final clusters in seqobs format.  
  open(OUT, ">".$new_cluster_fasta) or die "Can't open ".$new_cluster_fasta." for writing\n";
  $ref_fasta_db = Iterator::FastaDb->new($clusters_fasta) or die("Unable to open Fasta file ".$clusters_fasta);
  while(my $curr = $ref_fasta_db->next_seq() ){
    $curr->header =~ m/(\d+)/;
    my $header = $1;
    print OUT ">".$header."\n".$curr->seq."\n";
  }
  close(OUT);

}
writeClusters($fasta, $outdir."/obs.tsv", $outdir."/obs.fasta");

## REMOVE TEMP FILES
sub END{
  local $?;
  system("rm ".$tmpdir." -rf");
}

exit;
