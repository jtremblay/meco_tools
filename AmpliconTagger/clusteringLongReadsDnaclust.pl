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
clusteringLongReadsDnaclust.pl

PURPOSE:
Clustering wrapper for PacBio.
To count occurrences of all sequences. Intended for amplicon libraries. 
This script will cluster 16S rRNA sequences in a fashion similar to R. 
Edgar's OTU pipe pipeline, but using dnaclust (i.e. open source) 
instead of USEARCH. Briefly, there is: 1) 100% dereplication. 2) 97% 
clustering. 3) removal of chimeras with UCHIME ref.

INPUT:
--infile_fastq <fastq>      : Sequence file in fastq format.
--infile_fasta <fasta>      : Sequence file in fasta format.
--ref_db <string>           : Reference fasta database for chimera hunting.
--barcodes <string>         : Barcodes file in fasta format. Optional: used
                              for puting names in final cluster tables
                              instead of barcodes sequences.
--num_threads <int>         : Number of threads. Only used in the 100%
                              dereplication step.
--first_round <float>       : First round of clustering id (default = 0.99)
--lowAbunCutoff <int>       : Only sequences >= lowAbunCutOff <int> will be
                              kept. So sequences < lowAbunCutOff are discarded.
--align_perc_cutoff <float> : For PacBio data (long reads), double check the chimera
                              output results - default parameters output weird things
                              when long reads.

OUTPUT:
--outdir <string>           : Output directory where files will be generated.

OPTIONS:
--verbose                   : print status messages to stdout

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

NOTE:
This clustering pipelne follows recommandations by Robert Edgar - drive5.com
But with dnaclust instead of usearch.

ENDHERE

my ($help, $infile_fastq, $infile_fasta, $ref_db, $start_at, $verbose, $outdir, $barcodes, $num_threads, $lowAbunCutOff, $first_round, $second_round, $align_perc_cutoff);
GetOptions(
  'infile_fasta=s'      => \$infile_fasta,
  'infile_fastq=s'      => \$infile_fastq,
  'barcodes=s'          => \$barcodes,
  'ref_db=s'            => \$ref_db,
  'outdir=s'            => \$outdir,
  'lowAbunCutOff=i'     => \$lowAbunCutOff,
  'align_perc_cutoff=f' => \$align_perc_cutoff,
  'num_threads=i'       => \$num_threads,
  'first_round=f'       => \$first_round,
  'second_round=f'      => \$second_round,
  'verbose'             => \$verbose,
  'start_at=i'          => \$start_at,
  'help'                => \$help
);

if ($help){ print $usage; exit;}

## Validate
die "--infile_fastq or --infile_fasta arg missing.\n" if(!defined($infile_fastq) and !defined($infile_fasta));
die "--ref_db arg is missing (fasta reference database).\n" unless($ref_db);
die "--outdir arg missing.\n" unless($outdir);
$num_threads = 1 unless($num_threads);
$first_round = 0.99 unless($first_round);
$second_round = 0.97 unless($second_round);
$lowAbunCutOff = 3 unless($lowAbunCutOff);
$align_perc_cutoff = 60 unless($align_perc_cutoff);
print STDERR "first round : $first_round\n";
print STDERR "second round : $second_round\n";
print STDERR "lowAbunCutOff round : $lowAbunCutOff\n";
print STDERR "second round : $align_perc_cutoff\n";

#print STDERR "Running clustering workflow with:\nFirst round perc. id = $first_round\nSecond round perc. id = $second_round\nLow abund. cutoff = $lowAbunCutOff\n";

## Check paths.
my $dereplicate = which('dereplicate.pl');
print "path of dereplicate:\t".$dereplicate."\n" if($verbose);
my $dnaclust = which('dnaclust'); chomp $dnaclust;
print "path of dnaclust:\t".$dnaclust."\n" if($verbose);
my $vsearch= which('vsearch'); chomp $vsearch;
print "path of vsearch:\t".$vsearch."\n" if($verbose);
my $chimera_parser= which('chimeraParser.pl'); chomp $chimera_parser;
print "path of chimera_parser:\t".$chimera_parser."\n" if($verbose);
die "Can't find dnaclust on path\n" if(!defined($dnaclust));
die "Can't find vsearch on path\n" if(!defined($vsearch));

$start_at = 0 unless($start_at);

##############
## BEG SUBS ##
##############
## Parse dclust table to generate a fasta file.
## Input dclust file
## Ref fasta file
## Fasta outfile
sub dclustToFasta{
  my ($dclust_file, $ref_fasta, $fasta_outfile) = @_;
  
  my %HoH = ();

  # Loop through dclust table.
  print STDERR "Looping through dclust file\n";
  open(IN, $dclust_file) or die "Can't open ".$dclust_file."\n";
  while(<IN>){
    chomp;
    
    my @row = split(/\t/, $_);

    my $totalSize = 0;
    my $i = 0;
    my $repId;
    foreach my $header (@row){
      
      my @header = split(/;/, $header);
      if($i == 0){
        $repId = shift @header;
        $HoH{$repId}{seqid} = $repId;
      }else{
        shift(@header);
      }
      my $clusterSize = pop @header; # We extract value, but not used.  

      # At this point only #ACCGT=1;#ACGGT=2;etc... is left
      foreach my $el (@header){

        $el       =~ m/#([ACGT]*)/;
        my $barcodeSeq   = $1;
        #$HoH{$repId}{$barcodeSeq} = $size;
  
        $el        =~ m/=(\d+)/;
        my $size     = $1;
        $HoH{$repId}{$barcodeSeq} += $size;
        
        $totalSize     = $totalSize + $size;
        
      }      
      $HoH{$repId}{size} = $totalSize;
      $i++;
    }
  }  
  close(IN);
  
  # Loop through fasta reference to grab actual sequence representative for each cluster.
  my $ref_fasta_db = Iterator::FastaDb->new($ref_fasta) or die("Unable to open Fasta file ".$ref_fasta."\n");
  while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/>//;
    my @row = split(/;/, $header);
    my $seqid = shift(@row);
    $HoH{$seqid}{seq} = $curr->seq if(exists $HoH{$seqid});
  }
  
  # Open fasta outfile for writing.
  open(OUT, ">".$fasta_outfile) or die "Can't open ".$fasta_outfile."\n";

  # Loop through hash and print.
  foreach my $cluster_id (sort {$HoH{$b}{size} <=> $HoH{$a}{size} } keys %HoH){
    print OUT ">".$cluster_id.";";
    
    for my $barcode_seq (keys %{ $HoH{$cluster_id} }){
      next if( ($barcode_seq eq "size") or ($barcode_seq eq "seqid") or ($barcode_seq eq "seq")  ); #if size or seqid or seq.
      print OUT "#".$barcode_seq."=".$HoH{$cluster_id}{$barcode_seq}.";";
    }
    
    #print OUT $HoH{$cluster_id}{seqid};
    print OUT "size=".$HoH{$cluster_id}{size}."\n";
    print OUT $HoH{$cluster_id}{seq}."\n";
  }
  close(OUT);  
}

## Write clusters table.
## Input: clusters fasta sequence rep.
## Output: cluster table outfile and cluster fasta file.
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
      if($_ =~ m/\#([ACGT]*)=(\d+)/){
        $HoH{$seqid}{$1} = $2;
        $barcode_hash{$1} = 1;
      }
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
      $header =~ s/\t//g;
      $names{$curr->seq} = $header;
    }    
  }
  
  print OUT "#CLUSTER";
  foreach my $ref_barcode (sort keys %barcode_hash) {
    if($barcodes){
      #if(defined $names{$ref_barcode} and $names{$ref_barcode} ne "" and exists ($names{$ref_barcode}) ){
        print OUT "\t".$names{$ref_barcode};
      #}
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
my $sr = $second_round;
$sr =~ s/\.//g;

## dereplicate reads at 100%
my $cmd = $dereplicate;
$cmd .= " --fasta ".$infile_fasta if($infile_fasta);
$cmd .= " --fastq ".$infile_fastq if($infile_fastq);
$cmd .= " --outfile ".$outdir."/derep1.fasta";
$cmd .= " --num_threads ".$num_threads;
print STDERR $cmd."\n" if($start_at <= 1);
system($cmd) if($start_at <= 1);
die "command failed: $!\n" if($? != 0);

## Cluster at 97%
$cmd = $dnaclust." --no-k-mer-filter";
$cmd .= " -t ".$num_threads;
$cmd .= " -s ".$first_round;
$cmd .= " -i ".$outdir."/derep1.fasta"; 
$cmd .= " > ".$outdir."/derep1_$fr.dc";
print STDERR $cmd."\n" if($start_at <= 2);
system($cmd) if($start_at <= 2);
die "command failed: $!\n" if($? != 0);

## Parse dnaclust outfile and write clusters
dclustToFasta($outdir."/derep1_$fr.dc", $outdir."/derep1.fasta", $outdir."/derep1_$fr.fasta") if($start_at <= 3);

## Sort by size
$cmd = $dereplicate;
$cmd .= " --fasta ".$outdir."/derep1_$fr.fasta";
$cmd .= " --outfile ".$outdir."/derep1_".$fr."_derep2.fasta";
$cmd .= " --minsize ".$lowAbunCutOff;
$cmd .= " --sort ";
$cmd .= " --num_threads ".$num_threads;
print STDERR $cmd."\n" if($start_at <= 3);
system($cmd) if($start_at <= 3);
die "command failed: $!\n" if($? != 0);

#die $outdir."/derep1_099_derep2.fasta is empty.\n" if(-s $outdir."/derep1_099_derep2.fasta is empty.");

#Chimera detection (UCHIME de novo mode):
#  usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
#     [-uchimeout results.uch] [-uchimealns results.alns]
#  Input is estimated amplicons with integer abundances specified using ";size=N".
#usearch -uchime_denovo

#$cmd = $usearch;
#$cmd .= " -uchime_denovo ".$outdir."/derep1_099_derep2.fasta";
#$cmd .= " -nonchimeras ".$outdir."/derep1_099_derep2_nonChimDeNovo.fasta";
##$cmd .= " -threads ".$num_threads;
#system($cmd) if($start_at <= 5);
#die "command failed: $!\n" if($? != 0);

#usearch -uchime_ref
#Chimera detection (UCHIME ref. db. mode):
#  usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
#    [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]
$cmd = $vsearch;
$cmd .= " --uchime_ref ".$outdir."/derep1_".$fr."_derep2.fasta";
$cmd .= " --db ".$ref_db;
$cmd .= " --nonchimeras ".$outdir."/derep1_".$fr."_derep2_nonChimRef.fasta";
$cmd .= " --threads ".$num_threads;
$cmd .= " --uchimeout $outdir/uchime_ref.tsv";
print STDERR $cmd."\n" if($start_at <= 4);
system($cmd) if($start_at <= 4);
die "command failed: $!\n" if($? != 0);

#$cmd = $chimera_parser;
$cmd = "chimeraParser.pl";
$cmd .= " --infile ".$outdir."/uchime_ref.tsv";
$cmd .= " --infile_fasta ".$outdir."/derep1_".$fr."_derep2_nonChimRef.fasta";
$cmd .= " --cutoff $align_perc_cutoff > ".$outdir."/derep1_".$fr."_derep2_nonChimRef_parsed.fasta 2> ".$outdir."/chimera.fasta";
print STDERR $cmd."\n" if($start_at <= 5);
system($cmd) if($start_at <= 5);
die "command failed: $!\n" if($? != 0);

## Sort by size
$cmd = $dereplicate;
$cmd .= " --fasta ".$outdir."/derep1_".$fr."_derep2_nonChimRef_parsed.fasta";
$cmd .= " --outfile ".$outdir."/derep1_".$fr."_derep2_nonChimRef_sorted.fasta";
$cmd .= " --sort";
$cmd .= " --num_threads ".$num_threads;
print STDERR $cmd."\n" if($start_at <= 6);
system($cmd) if($start_at <= 6);
die "command failed: $!\n" if($? != 0);

## Cluster at 97%
#$cmd = $dnaclust." --no-k-mer-filter";
#$cmd .= " -t ".$num_threads;
#$cmd .= " -s ".$second_round;
#$cmd .= " -i ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted.fasta"; 
#$cmd .= " > ".$outdir."/derep1_099_derep2_nonChimDenovoRef_sorted_097.dc";
#system($cmd) if($start_at <= 8);
#die "command failed: $!\n" if($? != 0);

## Parse Uclust file and write clusters
#dclustToFasta(
#  $outdir."/derep1_099_derep2_nonChimDenovoRef_sorted_097.dc",   # DC file 
#  $outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted.fasta",    # Ref file
#  $outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097.fasta" # new file
#) if($start_at <= 9);

## Sort by size
#$cmd = $dereplicate;
#$cmd .= " --fasta ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097.fasta";
#$cmd .= " --outfile ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097_sorted.fasta";
#$cmd .= " --sort";
#$cmd .= " --num_threads ".$num_threads;
#system($cmd) if($start_at <= 10);
#die "command failed: $!\n" if($? != 0);

## Generate Obs table   derep1_".$fr."_derep2_nonChimRef_sorted.fasta
writeClusters($outdir."/derep1_".$fr."_derep2_nonChimRef_sorted.fasta", $outdir."/obs.tsv", $outdir."/obs.fasta");

exit;
