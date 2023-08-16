#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $gene_id = $row[0];
   
   if($_ =~ m/\[gene=(\S+)\]/){
    my $gene_tag = $1;
    $hash{$gene_id}{gene_tag} = $gene_tag;
   }else{
    $hash{$gene_id}{gene_tag} = "undefined";
   
   }
   
   if($_ =~ m/\[locus_tag=(\S+)\]/){
    my $locus_tag = $1; 
    $hash{$gene_id}{locus_tag} = $locus_tag;
   }else{
    $hash{$gene_id}{locus_tag} = "undefined";
   }

}

print STDERR Dumper(\%hash);

close(IN);
foreach my $key (keys %hash){
    print STDOUT $key."\t".$hash{$key}{locus_tag}."\t".$hash{$key}{gene_tag}."\n";
}
exit;
