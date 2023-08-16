#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
addTaxonomyToMiRNA.pl

PURPOSE:

INPUT:
--infile <string>          : Parsed mirna file
--infile_link <string>     : File linking each contig to their corresponding fna file.
--infile_taxonomy <string> : taxonomy file.
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $infile_link, $infile_taxonomy);
my $verbose = 0;

GetOptions(
   'infile=s'          => \$infile,
   'infile_link=s'     => \$infile_link,
   'infile_taxonomy=s' => \$infile_taxonomy,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my %hash_link;

open(IN, "<".$infile_link) or die "Can't open $infile_link\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $filename = $row[0];
   my $contig_id = $row[1];
   
   $hash_link{$contig_id} = $filename;

}
close(IN);
print STDERR "link done\n";

open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $filename = $row[0];
   my $p = $row[5];
   my $c = $row[6];
   my $o = $row[7];
   my $f = $row[8];
   my $g = $row[9];
   my $s = $row[1];

   $hash{$filename}{phylum} = $p;
   $hash{$filename}{class} = $c;
   $hash{$filename}{order} = $o;
   $hash{$filename}{family} = $f;
   $hash{$filename}{genus} = $g;
   $hash{$filename}{species} = $s;

}
close(IN);
print STDERR "Metadata done\n";

#print STDERR Dumper(\%hash);
#exit;

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    if($. == 1){
        print STDOUT $_."\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";
        next;
    }
    my @row = split(/\t/, $_);
    my $contig_id = $row[1];
    #my @contig_id = split(/_/, $contig_id);
    #pop(@contig_id);
    #$contig_id = join("_", @contig_id);
    my $filename = $hash_link{$contig_id};
    print STDOUT $_."\t";
    print STDOUT $hash{$filename}{phylum}."\t";
    print STDOUT $hash{$filename}{class}."\t";
    print STDOUT $hash{$filename}{order}."\t";
    print STDOUT $hash{$filename}{family}."\t";
    print STDOUT $hash{$filename}{genus}."\t";
    print STDOUT $hash{$filename}{species}."\n"; 
}
close(IN);

