#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use List::Util qw(sum);

my $usage=<<'ENDHERE';
NAME:
generateTaxonomyFromCAT.pl

PURPOSE:

INPUT:
--infile_taxonomy <string>  : Output of CAT
				
OUTPUT:
STDOUT <string>             : More concised summary taxonomy file.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxonomy);
my $verbose = 0;

GetOptions(
   'infile_taxonomy=s'  => \$infile_taxonomy,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## MAIN
print STDOUT "contig_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";

my %hash;
open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";

my $seen = "true";
my $old_id = "initial";
my $old_average = -1;
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    my @row = split(/\t/, $_);
    
    my $kingdom = "k__NULL";
    my $phylum = "p__NULL";
    my $class = "c__NULL";
    my $order = "o__NULL";
    my $family = "f__NULL";
    my $genus = "g__NULL";
    my $species = "s__NULL";
    
    if($row[1] !~ m/no taxid assigned/){

        my ($match) = $_ =~ m/(root.*)$/;
        my @row2 = split(/\t/, $match);

        foreach my $el (@row2){
            if($el =~ m/superkingdom/){
                $kingdom = "k__".(split /\(/, $el)[0];
                chop($kingdom); 
            }elsif($el =~ m/phylum/){
                $phylum = "p__".(split /\(/, $el)[0];
                chop($phylum); 
            }elsif($el =~ m/class/){
                $class = "c__".(split /\(/, $el)[0];
                chop($class); 
            }elsif($el =~ m/order/){
                $order = "o__".(split /\(/, $el)[0];
                chop($order); 
            }elsif($el =~ m/family/){
                $family = "f__".(split /\(/, $el)[0];
                chop($family); 
            }elsif($el =~ m/genus/){
                $genus = "g__".(split /\(/, $el)[0];
                chop($genus); 
            }elsif($el =~ m/species/){
                $species = "s__".(split /\(/, $el)[0];
                chop($species); 
            }
        }
        $kingdom =~ s/k__//;
        $phylum =~ s/p__//;
        $class =~ s/c__//;
        $order =~ s/o__//;
        $family =~ s/f__//;
        $genus =~ s/g__//;
        $species =~ s/s__//;
        
        # Finally, chose if result is to be printed.
        my $curr_id = $row[0];
        my @values = split(/;/, $row[4]);
        my $curr_average = mean(@values);

        if($curr_id ne $old_id){ # means we have a new entry. So we populate the hash by default
            $hash{$curr_id} = $row[0]."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species;
            $old_average = $curr_average
        }else{ # means we have an existing entry, we populate the hash if avg is higher than the last value.
            if($curr_average >= $old_average){
                $hash{$curr_id} = $row[0]."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species;
            }
            $old_average = $curr_average;
        }
        $old_id = $curr_id;
    }
    #STDOUT $row[0]."\t".$kingdom."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus."\t".$species."\n";
}
close(IN);

#print STDERR Dumper(\%hash);
for my $key (keys %hash){
    print STDOUT $hash{$key}."\n";
}

exit;
    
sub mean {
    return sum(@_)/@_;
}
