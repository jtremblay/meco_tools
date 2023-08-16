#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateEmirgeOTUTable.pl

PURPOSE:

INPUT:
--fastas <string>     : Sequence file (list of, separated by a ",").
--taxonomies <string> : taxonomy files separated by a ","
--names <string>      : Name for each file provided in --infiles.
				            also separated by a comma.
--flag <string>       : Either 'rdp' or 'blastn'
--cutoff <float>        : cutoff for rdp = 0.50. blastn cutoff was already specified
                        elsewhere...

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $fastas, $taxonomies, $names, $flag, $cutoff);
my $verbose = 0;

GetOptions(
   'fastas=s'     => \$fastas,
   'taxonomies=s' => \$taxonomies,
   'names=s'      => \$names,
   'cutoff=f'     => \$cutoff,
   'flag=s'       => \$flag,
   'verbose' 	   => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

die "--flag has to either 'rdp' or 'blastn'" unless($flag eq "rdp" || $flag eq "blastn");

my %hash;

## MAIN
	
$cutoff = 0.50 unless($cutoff);

my @files = split(/,/, $taxonomies);
my @names = split(/,/, $names);
my @fastas = split(/,/, $fastas);
my $number_of_samples = (@files);

foreach my $file (@files){
   my $name = shift(@names);
   my $fasta = shift(@fastas);

   my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
   while( my $curr = $ref_fasta_db->next_seq() ) {
       my ($id, $abundance) = $curr->header =~ m/>(\d+\|\S+) Prior=(\d+\.\d+)/g;
       $id = "$name-$id";
       $abundance = $abundance * 100;
       $hash{$name}{$id}{abundance} = sprintf("%.1f", $abundance); 
       $hash{$name}{$id}{id} = $id;
   }
	
   if($flag eq "rdp"){
	   open(IN, "<".$file) or die "Can't open $file\n";
	   while(<IN>){
	      chomp;
	      my @row = split(/\t/, $_);
	      my $id = "$name-$row[0]";
	      # 5,6,7 = k__bact, Kingdom, score
	      # 8,9,10 = phy
	      # 11,12,13 = class
	      # 14,15,16  = order
	      # 17,18,19 = family
	      # 20,21,22 = genus
	      # 23,24,25 = specie
	      
	      my $lineage = "";
	      if($row[7] >= $cutoff){
	         $lineage .= "$row[5];";
	      }
	      if($row[10] >= $cutoff){
	         $lineage .= "$row[8];";
	      }
	      if($row[13] >= $cutoff){
	         $lineage .= "$row[11];";
	      }
	      if($row[16] >= $cutoff){
	         $lineage .= "$row[14];";
	      }
	      if($row[19] >= $cutoff){
	         $lineage .= "$row[17];";
	      }
	      if($row[22] >= $cutoff){
	         $lineage .= "$row[20];";
	      }
	      if($row[25] >= $cutoff){
	         $lineage .= "$row[23];";
	      }
	      next if $lineage eq "";
	
	      $hash{$name}{$id}{lineage} = $lineage;
	   }
	}elsif($flag eq "blastn"){
   
       #$cutoff = "1e-02" unless($cutoff);
	   
      open(IN, "<".$file) or die "Can't open $file\n";
	   while(<IN>){
	      chomp;
         next if($_ =~ m/^gene_id/);
	      my @row = split(/\t/, $_);
	      my $id = "$name-$row[0]";
         my $lineage = "k__$row[2];p__$row[3];c__$row[4];o__$row[5];f__$row[6];g__$row[7];s__$row[8];";
	      $hash{$name}{$id}{lineage} = $lineage;
      }
      close(IN);
   }
}   

#print STDERR Dumper(\%hash);
	
# print header
my $header = "#OTU ID";
for my $k1 (keys %hash) {
   $header .= "\t".$k1;
}
$header .= "\ttaxonomy\n";
print STDOUT $header;


my $curr_pos = 0;
#sample
for my $k1 (keys %hash) {

   # rRNA entries
   for my $k2 (keys %{ $hash{$k1} }){
      print STDOUT $hash{$k1}{$k2}{id};
      
      for(my $i=0; $i<$number_of_samples; $i++){
         if($curr_pos == $i){
            print STDOUT "\t$hash{$k1}{$k2}{abundance}";
         }else{
            print STDOUT "\t0";
         }
      }
      print STDOUT "\t".$hash{$k1}{$k2}{lineage};
      print STDOUT "\n";
   }
   $curr_pos++;
}

exit;
