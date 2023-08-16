#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile_fasta <string> : Sequence file
--infile_tax <string>   : Tsv taxonomy file

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_tax);
my $verbose = 0;

GetOptions(
   'infile_fasta=s' => \$infile_fasta,
   'infile_tax=s'   => \$infile_tax,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile_tax) or die "Can't open $infile_tax\n";
while(<IN>){
   chomp;
   next if($. == 1);
   my @row = split(/\t/, $_);
   my $header = shift(@row);
   my $k = $row[0];
   my $p = $row[1];
   my $c = $row[2];
   my $o = $row[3];
   my $f = $row[4];
   my $g = $row[5];
   my $s = $row[6];
   my $lineage = "k__".$k.";";
   $lineage .= "p__".$p.";";
   $lineage .= "c__".$c.";";
   $lineage .= "o__".$o.";";
   $lineage .= "f__".$f.";";
   $lineage .= "g__".$g.";";
   $lineage .= "s__".$s."";

   $hash{$header} = $lineage;

}
close(IN);
print STDERR Dumper(\%hash);


my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my ($header) = $curr->header =~ m/^>(\S+)/;
    #my @header = split(/__/, $header);
    #$header = $header[0];
    #print STDERR $header."\n";
    if(exists $hash{$header}){
        print STDOUT ">".$header." ".$hash{$header}."\n".$curr->seq."\n";
    }
}
exit;

