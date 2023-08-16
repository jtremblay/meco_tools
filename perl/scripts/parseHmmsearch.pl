#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
parseHmmsearch.pl

PURPOSE:
Takes a hmmsearch output and parse it.

INPUT:
--infile <string> : .tblout file
--n <int>         : Number of hits to keep (default:1).
--e <int>         : Evalue cutoff. Default = 1e-10
--q <int>         : minimum query length. Default = 100
--a <int>         : minimum alignment length. Default = 100

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $n, $e, $q, $a);
my $verbose = 0;

GetOptions(
  'infile=s'  => \$infile,
  'n=i'       => \$n,
  'e=f'       => \$e,
  'q=i'       => \$q,
  'a=i'       => \$a,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $header = "";
my %hash;
#my $i=1;
#my $j=0;

$a = 100 unless($a);
$q = 100 unless($q);
$e = 1e-10 unless($e);
$n = 1 unless($n);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    if($_ =~ m/#/g){
        next;
    }

    my @row = split(/\s+/, $_);

    my $gene_id = $row[2];
    my $pfam = $row[1];
    my $evalue = $row[4];

    if($evalue < $e){

        if(exists $hash{$gene_id}){
            if($evalue < $hash{$gene_id}{evalue}){
                $hash{$gene_id}{id} = $_;
                $hash{$gene_id}{evalue} = $evalue;
            } 

        }else{
            $hash{$gene_id}{id} = $_;
            $hash{$gene_id}{evalue} = $evalue;
        }  
    }
}

print STDOUT $header;

for my $k1 (sort keys %hash) { 
  print STDOUT $hash{$k1}{id}."\n";
  #for my $k2 (sort keys %{ $hash{$k1} }){
  #  print STDOUT $hash{$k1}{$k2}."\n";

  #}
}
exit;
