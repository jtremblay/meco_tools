#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string>     : Sequence file
--annotation <string> : annotation file directly download from CCDS ftp.
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $annotations, $outfileFasta);
my $verbose = 0;

GetOptions(
    'infile=s' 	    => \$infile,
    'outfileFasta=s'=> \$outfileFasta,
    'annotations=s' => \$annotations,
    'verbose' 	    => \$verbose,
    'help' 		      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT, ">".$outfileFasta) or die "Can't open $outfileFasta\n";

# Parse annotation.
my %hash;
open(IN, "<".$annotations) or die "Can't open $annotations\n";
while(<IN>){
  chomp;
  my @row = split(/\t/, $_);
  my $chromo = $row[0];   #chromosome
  my $nc_access = $row[1];#nc_accession
  my $gene = $row[2];
  my $gene_id = $row[3];
  my $ccds_id = $row[4];
  my $ccds_status = $row[5];
  my $ccds_strand = $row[6];
  my $ccds_from = $row[7];
  my $ccds_to = $row[8];
  my $ccds_locations = $row[9];
  my $match_type = $row[10];
  
  $ccds_id .= ".X" if($chromo eq "X");
  $ccds_id .= ".Y" if($chromo eq "Y");
  $hash{$ccds_id}{'chr'}     = $chromo;
  $hash{$ccds_id}{'gene'}    = $gene;
  $hash{$ccds_id}{'geneID'}  = $gene_id;
  $hash{$ccds_id}{'strand'}  = $ccds_strand; 
}
close(IN);

print STDERR "Done building annotation hash...\n";

# Parse fasta file.
my %last;
my %start;
my %end;

my %seqs;

my $Ns = "";
for(my $i=0; $i<1000;$i++){
  $Ns .= "N";
}
#print STDERR length($Ns)."\n";
#exit;

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header;
  $header =~ s/>//g;
  my @row = split(/\|/, $header); 

  #print STDERR "header: $header\n"; 

  my $ccds_id = $row[0];
  my $name    = $row[1];
  my $chromo  = $row[2];

  my $chromoId = $chromo;
  $chromoId =~ s/chr//g;
  if(exists $seqs{$chromoId}){
    $seqs{$chromoId}{'seq'} .= $curr->seq.$Ns; 

  }else{
    $seqs{$chromoId}{'seq'} = $curr->seq.$Ns; 
  }

  #print STDERR "ccds_id $ccds_id\n";
  #print STDERR "name $name\n";
  #print STDERR "chromo $chromo\n";

  my $currChromo = $chromo;
  if(exists $start{$chromo}){
    $start{$chromo} = $last{$chromo};
    $last{$chromo}  = $start{$chromo} + length($curr->seq) + length($Ns);
    $end{$chromo}   = length($curr->seq) + $end{$chromo} + length($Ns); 
  }else{
    $start{$chromo} = 1; 
    $last{$chromo}  = length($curr->seq) + 1 + length($Ns);
    $end{$chromo}   = length($curr->seq) + length($Ns);
  }

  if($chromo =~ m/chrX/){
    $ccds_id .= ".X";
  }
  if($chromo =~ m/chrY/){
    $ccds_id .= ".Y"; 
  }
  if(exists $hash{$ccds_id}){
    print STDOUT "chr".$hash{$ccds_id}{'chr'}."\tprotein_coding\texon\t".$start{$chromo}."\t".($end{$chromo} - length($Ns))."\t.\t".$hash{$ccds_id}{'strand'}."\t.\t gene_id \"".$ccds_id."\"; transcript_id \"".$ccds_id."T\"; exon_number \"1\"; gene_name \"".$hash{$ccds_id}{'gene'}."\"; gene_biotype \"protein_coding\"; transcript_name \"".$hash{$ccds_id}{'gene'}."T\";\n";
  }    
}
close(IN);

foreach my $key (sort {$a <=> $b} keys %seqs){ 
  print OUT ">chr".$key."\n".$seqs{$key}{'seq'}."\n";
}

