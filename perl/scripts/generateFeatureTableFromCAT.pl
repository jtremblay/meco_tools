#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateOTUTableFromCAT.pl

PURPOSE:

INPUT:
--infile_taxonomy <string>  : Output of CAT
--infile_abundance <string> : Taxonomy abundance of each contig
				
OUTPUT:
STDOUT                      : OTU table format file.
--outfile_taxonomy <string> : More concised summary taxonomy file.

NOTES:
If -f paramter of CAT is too low, ambiguous classification can be given.
a simple score average will be perform here.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxonomy, $infile_abundance, $outfile_taxonomy);
my $verbose = 0;

GetOptions(
   'infile_taxonomy=s'  => \$infile_taxonomy,
   'infile_abundance=s' => \$infile_abundance,
   'outfile_taxonomy=s' => \$outfile_taxonomy,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT, ">".$outfile_taxonomy) or die "Can't open $outfile_taxonomy\n";
print OUT "contig_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";

my %hash;
open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    my @row = split(/\t/, $_);
    
    my $superkingdom = "k__NULL";
    my $kingdom = "k__NULL";
    my $phylum = "p__NULL";
    my $class = "c__NULL";
    my $order = "o__NULL";
    my $family = "f__NULL";
    my $genus = "g__NULL";
    my $species = "s__NULL";
    

    # Here superkingdom will = Bacteria, Archeaea. However Kingdom = Fungi when that fungal organism.

    if($row[1] eq "classified" || $row[1] eq "taxid assigned"){
        my ($match) = $_ =~ m/(root.*)$/;
        my @row2 = split(/\t/, $match);

        #print STDERR "Parsing $match\n";
        foreach my $el (@row2){
            if($el =~ m/\(superkingdom\)/){
                $superkingdom = "k__".(split /\(/, $el)[0];
                chop($superkingdom);
            }elsif($el =~ m/\(kingdom\)/){
                $kingdom = "k__".(split /\(/, $el)[0];
                chop($kingdom); 
            }elsif($el =~ m/\(phylum\)/){
                $phylum = "p__".(split /\(/, $el)[0];
                chop($phylum); 
            }elsif($el =~ m/\(class\)/){
                $class = "c__".(split /\(/, $el)[0];
                chop($class); 
            }elsif($el =~ m/\(order\)/){
                $order = "o__".(split /\(/, $el)[0];
                chop($order); 
            }elsif($el =~ m/\(family\)/){
                $family = "f__".(split /\(/, $el)[0];
                chop($family); 
            }elsif($el =~ m/\(genus\)/){
                $genus = "g__".(split /\(/, $el)[0];
                chop($genus); 
            }elsif($el =~ m/\(species\)/){
                $species = "s__".(split /\(/, $el)[0];
                chop($species); 
            }
        }
    }

    if($superkingdom eq "k__Bacteria" || $superkingdom eq "k__Archaea"){
        $kingdom = $superkingdom;
    }elsif($superkingdom eq "k__Viruses"){
        $kingdom = $superkingdom;
    }

    $hash{$row[0]} = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
    $superkingdom =~ s/k__//;
    $phylum =~ s/p__//;
    $class =~ s/c__//;
    $order =~ s/o__//;
    $family =~ s/f__//;
    $genus =~ s/g__//;
    $species =~ s/s__//;
    print OUT $row[0]."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
}
close(OUT);
close(IN);

#print STDERR Dumper(\%hash);

# Then print OTU table.
open(IN, "<".$infile_abundance) or die "Can't open $infile_abundance\n";
while(<IN>){
    chomp;
    if($. == 1){
        my @row = split(/\t/, $_);
        shift(@row);
        print STDOUT "#FEATURE_ID\t".join("\t", @row)."\ttaxonomy\n";
        next;
    }
    my @row = split(/\t/, $_);

    if(exists $hash{$row[0]}){
       print STDOUT $_."\t".$hash{$row[0]}."\n"; 
    }
}
close(IN);
exit;
