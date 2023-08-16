#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateOTUTableReadsCentric.pl

PURPOSE:

INPUT:
--infiles <string> : Sequence file (list of, separated by a ",").
--names <string>   : Name for each file provided in --infiles.
				         also separated by a comma.
--cutoff <float>   : RDP score cutoff to reconstitute lineages.
                     lineages will be reconstructed up to the 
                     --cutoff value. Default = 0.5.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $names, $cutoff);
my $verbose = 0;

GetOptions(
   'infiles=s' => \$infiles,
   'names=s'   => \$names,
   'cutoff=f'  => \$cutoff,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

$cutoff = 0.5 unless($cutoff);

my %hash;

## MAIN
my @files = split(/,/, $infiles);
my @names = split(/,/, $names);
my $number_of_samples = (@files);


foreach my $file (@files){
   my $name = shift(@names);

   open(IN, "<".$file) or die "Can't open $file\n";
   while(<IN>){
      chomp;
      my @row = split(/\t/, $_);
      my $id = $row[0];
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

      if(exists $hash{$name}{$lineage}){
         $hash{$name}{$lineage}{abundance}++;
      }else{
         $hash{$name}{$lineage}{abundance} = 0;
         $hash{$name}{$lineage}{id} = $id;
      }
   }
}


#print STDERR Dumper(\%hash);

# print header
my $header = "#OTU ID";
for my $k1 (keys %hash) {
   $header .= "\t".$k1;
    #for my $k2 (keys %{ $hash{$k1} }){
    #}
}
$header .= "\ttaxonomy\n";
print STDOUT $header;


my $curr_pos = 0;
for my $k1 (keys %hash) {

   for my $k2 (keys %{ $hash{$k1} }){
      print STDOUT $hash{$k1}{$k2}{id};
      
      for(my $i=0; $i<$number_of_samples; $i++){
         if($curr_pos == $i){
            print STDOUT "\t$hash{$k1}{$k2}{abundance}";
         }else{
            print STDOUT "\t0";
         }
      }
      print STDOUT "\t".$k2;
      print STDOUT "\n";
   }
   $curr_pos++;
}
