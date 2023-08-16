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

$| = 1;    # UNBUFFER OUTPUT

my $usage = <<'ENDHERE';
NAME:
clusteringLongReadsUsearch.pl

PURPOSE:
To count occurrences of all sequences and optionally bin by abundance.
Intended for amplicon libraries. This script will cluster sequences at
97% id and remove chimeras wit UCHIME (Vsearch) (both denovo and ref).

INPUT:
--infile_fastq <fastq>  : Sequence file in fastq format.
--infile_fasta <fasta>  : Sequence file in fasta format.
--ref_db <string>       : Reference fasta database for chimera hunting.
--barcodes <string>     : Barcodes file in fasta format. Optional: used
                          for puting names in final cluster tables
                          instead of barcodes sequences.
--num_threads           : Number of threads. Only used in the 100%
                          dereplication step.
--first_round <float>   : First round of clustering id (default = 0.99)
--lowAbunCutoff <int>   : Only sequences >= lowAbunCutOff <int> will be
                          kept. So sequences < lowAbunCutOff are discarded.

OUTPUT:
--outdir <string>       : Output directory where files will be generated.

OPTIONS:
--verbose               : print status messages to stdout

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

NOTE:

ENDHERE

my ($help, $infile_fastq, $infile_fasta, $ref_db, $start_at, $verbose, $outdir, $barcodes, $num_threads, $first_round, $low_abun_cutoff);
GetOptions(
  'infile_fasta=s'  => \$infile_fasta,
  'infile_fastq=s'  => \$infile_fastq,
  'barcodes=s'      => \$barcodes,
  'ref_db=s'        => \$ref_db,
  'outdir=s'        => \$outdir,
  'num_threads=i'   => \$num_threads,
  'first_round=f'   => \$first_round,
  'lowAbunCutoff=i' => \$low_abun_cutoff,
  'verbose'         => \$verbose,
  'start_at=i'      => \$start_at,
  'help'            => \$help
);

if ($help){ print $usage; exit;}

## Validate
die "--infile_fastq or --infile_fasta arg missing.\n" if(!defined($infile_fastq) and !defined($infile_fasta));
die "--ref_db arg is missing (fasta reference database).\n" unless($ref_db);
die "--outdir arg missing.\n" unless($outdir);
$num_threads = 1 unless($num_threads);

## Check paths.
my $dereplicate = which('dereplicate.pl');
print "path of dereplicate script :\t".$dereplicate."\n" if($verbose);
my $vsearch = which('vsearch');
print "path of vsearch:\t".$vsearch."\n" if($verbose);
$start_at = 0 unless($start_at);

##############
## BEG SUBS ##
##############
## Parse uclust table to generate a fasta file.
## Input uclust file
## Ref fasta file
## Fasta outfile
sub uclustToFasta{
  my ($uclust_file, $ref_fasta, $fasta_outfile) = @_;
  
  my %HoH = ();
  my %hash_seq = ();

  # Loop through uclust table.
  open(IN, $uclust_file) or die "Can't open ".$uclust_file."\n";
  while(<IN>){
    chomp;
    next if( substr($_,0,1) eq "#");
  
    my @row = split(/\t/, $_);
    my $status = $row[0];
    my $cluster_id = $row[1];
    my $query = $row[8];

    my @query = split(/;/, $query);
    my $seqid = shift(@query);
    my $size_field = pop(@query);
    $size_field =~ m/size=(\d+)/;      
    my $size = $1;
  
    next if($status eq "C");  
    $HoH{$cluster_id}{size} += $size;
    $HoH{$cluster_id}{seqid} = $seqid if($status eq "S");
    
    # Parse barcodes
    foreach my $barcode_field (@query){
      $barcode_field =~ m/#([ACGT]*)/;
      my $barcode_seq = $1;
      $barcode_field =~ m/[ACGT]*=(\d+)/;
      my $barcode_count = $1;
      $HoH{$cluster_id}{$barcode_seq} += $barcode_count;
    }
  }
  close(IN);
  
  # Loop through fasta reference to grab actual sequence.
  my $ref_fasta_db = Iterator::FastaDb->new($ref_fasta) or die("Unable to open Fasta file ".$ref_fasta."\n");
  while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/>//;
    my @row = split(/;/, $header);
    my $seqid = shift(@row);
    $hash_seq{$seqid} = $curr->seq;    
  }
  
  # Open fasta outfile for writing.
  open(OUT, ">".$fasta_outfile) or die "Can't open ".$fasta_outfile."\n";

  # Loop through hash and print.
  foreach my $cluster_id (sort {$HoH{$b}{size} <=> $HoH{$a}{size} } keys %HoH){
    print OUT ">".$cluster_id.";";
    
    for my $barcode_seq (keys %{ $HoH{$cluster_id} }){
      next if( ($barcode_seq eq "size") or ($barcode_seq eq "seqid") ); #if size or seqid
      print OUT "#".$barcode_seq."=".$HoH{$cluster_id}{$barcode_seq}.";";
    }
    print OUT "size=".$HoH{$cluster_id}{size}."\n";
    print OUT $hash_seq{ $HoH{$cluster_id}{seqid} }."\n";
  }
  close(OUT);  
}

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
      $_ =~ m/\#([ACGT]*)=(\d+)/;
      $HoH{$seqid}{$1} = $2;
      $barcode_hash{$1} = 1;
    }
  }
  
  # Open outfile for writing.
  open(OUT, ">".$clusters_table) or die "Can't open ".$clusters_table."\n";
  
  # Print header
  my %names = ();
  if($barcodes){
    $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file ".$barcodes."\n");
    while( my $curr = $ref_fasta_db->next_seq() ) {
      my $header = $curr->header;
      $header =~ s/>//;
      $names{$curr->seq} = $header;
    }    
  }
  
  print OUT "#CLUSTER";
  foreach my $ref_barcode (sort keys %barcode_hash) {
    if($barcodes){
      print OUT "\t".$names{$ref_barcode};
    }else{
      print OUT "\t".$ref_barcode;
    }
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
##############
## END SUBS ##
##############

## MAIN
my $fr = $first_round;
$fr =~ s/\.//g;
#my $sr = $second_round;
#$sr =~ s/\.//g;

## dereplicate reads at 100%
my $cmd = $dereplicate;
$cmd .= " --fasta ".$infile_fasta if($infile_fasta);
$cmd .= " --fastq ".$infile_fastq if($infile_fastq);
$cmd .= " --outfile ".$outdir."/derep1.fasta";
$cmd .= " --num_threads ".$num_threads;
if($start_at <= 1){
    print STDERR $cmd."\n";
    system($cmd);
    die "command failed: $!\n" if($? != 0);
}

## Cluster at 97%
$cmd = $vsearch." --cluster_smallmem ";
$cmd .= " ".$outdir."/derep1.fasta ";
$cmd .= " --id ".$first_round;
$cmd .= " --uc ".$outdir."/derep1_".$fr.".uc";
$cmd .= " --log ".$outdir."/derep1_".$fr.".log";
$cmd .= " --usersort";
if($start_at <= 2){
    print STDERR $cmd."\n";
    system($cmd);
    die "command failed: $!\n" if($? != 0);
}

## Parse Uclust file and write clusters
uclustToFasta($outdir."/derep1_".$fr.".uc", $outdir."/derep1.fasta", $outdir."/derep1_".$fr.".fasta");

## Sort by size
$cmd = $dereplicate;
$cmd .= " -fasta ".$outdir."/derep1_".$fr.".fasta";
$cmd .= " -outfile ".$outdir."/derep1_".$fr."_derep2.fasta";
$cmd .= " -minsize ".$low_abun_cutoff;
$cmd .= " -sort ";
if($start_at <= 3){
    print STDERR $cmd."\n";
    system($cmd);
    die "command failed: $!\n" if($? != 0);
}
#Chimera detection (UCHIME de novo mode):
#  usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
#     [-uchimeout results.uch] [-uchimealns results.alns]
#  Input is estimated amplicons with integer abundances specified using ";size=N".
#usearch -uchime_denovo
#$cmd = $usearch;
#$cmd .= " -uchime_denovo ".$outdir."/derep1_".$fr."_derep2.fasta";
#$cmd .= " -nonchimeras ".$outdir."/derep1_".$fr."_derep2_nonChimDeNovo.fasta";
##$cmd .= " -threads ".$num_threads;
#system($cmd) if($start_at <= 4);
#die "command failed: $!\n" if($? != 0);

#usearch -uchime_ref
#Chimera detection (UCHIME ref. db. mode):
#  usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
#    [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]
$cmd = $vsearch;
$cmd .= " --uchime_ref ".$outdir."/derep1_".$fr."_derep2.fasta";
$cmd .= " --db ".$ref_db;
$cmd .= " --strand both";
$cmd .= " --nonchimeras ".$outdir."/derep1_".$fr."_derep2_nonChimRef.fasta";
$cmd .= " --threads ".$num_threads;
if($start_at <= 4){
    print STDERR $cmd."\n";
    system($cmd);
    die "command failed: $!\n" if($? != 0);
}

## Sort by size
$cmd = $dereplicate;
$cmd .= " --fasta ".$outdir."/derep1_".$fr."_derep2_nonChimRef.fasta";
$cmd .= " --outfile ".$outdir."/derep1_".$fr."_derep2_nonChimRef_sorted.fasta";
$cmd .= " --sort";
if($start_at <= 5){
    print STDERR $cmd."\n";
    system($cmd);
    die "command failed: $!\n" if($? != 0);
}

## Generate Obs table
writeClusters($outdir."/derep1_".$fr."_derep2_nonChimRef_sorted.fasta", $outdir."/obs.tsv", $outdir."/obs.fasta");

exit;
