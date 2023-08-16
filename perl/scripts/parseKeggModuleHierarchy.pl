#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
#use Iterator::FastaDb;
#use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseKeggModuleHierarchy.pl

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
   'infile=s'   => \$infile,
   'verbose'    => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
my ($one, $two, $three, $module, $module_desc);
while(<IN>){
    chomp;
    if($_ =~ m/^ (\S.*)$/){
        $one = $1;
        print STDERR "one:$one\n";
        #$hash{$one} = $1;
    }
    if($_ =~ m/^   (\S.*)$/){
        $two = $1;
        print STDERR "two:$two\n";
        #$hash{$one}{$two} = $1;
    }
    if($_ =~ m/^     (\S.*)$/){
        #$hash{$one}{$two}{$three} = $1;
        $three = $1;
        print STDERR "three:$three\n";
    }
    if($_ =~ m/^       (M\d{5})  (.*$)/){
        $module = $1;
        $module_desc = $2;
        print STDERR "module:$module\n";
        #$hash{$one}{$two}{$three}{$module} = $module;
    }
    #K00855  PRK, prkB; phosphoribulokinase [EC:2.7.1.19]
    if($_ =~ m/^         (K\d{5})  (.*) \[(EC:.*)\]$/){
        my $K = $1;
        my $desc = $2;
        my $EC = $3;
        
        my @el = split(/;/, $desc);
        my $genes = $el[0];
        $desc = $el[1];
        
        $hash{$one}{$two}{$three}{$module}{$K}{mdesc} = $module_desc;
        $hash{$one}{$two}{$three}{$module}{$K}{desc} = $desc;
        $hash{$one}{$two}{$three}{$module}{$K}{EC} = $EC;
        $hash{$one}{$two}{$three}{$module}{$K}{symbol} = $genes;
        print STDERR "K:$K\n";

    }elsif($_ =~ m/^         (K\d{5})  (.*)$/){
        my $K = $1;
        my $desc = $2;
        my $EC = $3;
        
        my @el = split(/;/, $desc);
        my $genes = $el[0];
        $desc = $el[1];
        
        $hash{$one}{$two}{$three}{$module}{$K}{mdesc} = $module_desc;
        $hash{$one}{$two}{$three}{$module}{$K}{desc} = $desc;
        $hash{$one}{$two}{$three}{$module}{$K}{EC} = $EC;
        $hash{$one}{$two}{$three}{$module}{$K}{symbol} = $genes;
        print STDERR "K:$K\n";
    }

}
print STDERR Dumper(\%hash);
close(IN);

print STDOUT "#level1\tlevel2\tlevel3\tMODULE\tMODULE_NAME\tKEGG_ENTRY\tK_DESC\tSymbol\tENZYME_CODE\n";
for my $k1 (sort keys %hash) {
    for my $k2 (sort keys %{ $hash{$k1} }){
        for my $k3 (sort keys %{ $hash{$k1}{$k2} }) {
            for my $module (sort keys %{ $hash{$k1}{$k2}{$k3} }){
                for my $K (sort keys %{ $hash{$k1}{$k2}{$k3}{$module} }){
                    print STDOUT "$k1\t$k2\t$k3\t$module\t";
                    print STDOUT $hash{$k1}{$k2}{$k3}{$module}{$K}{mdesc}."\t";
                    print STDOUT "$K\t";
                    print STDOUT $hash{$k1}{$k2}{$k3}{$module}{$K}{desc}."\t";
                    print STDOUT $hash{$k1}{$k2}{$k3}{$module}{$K}{symbol}."\t";
                    print STDOUT $hash{$k1}{$k2}{$k3}{$module}{$K}{EC}."\n";
                }
            }
        }
    }
}


open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    if($_ =~ m/(M\d{5})/){
        print STDERR "Module: $1\n";
    }
}
exit;
