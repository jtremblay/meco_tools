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
--infile <string> : kegg ref pathway file htext format
				
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
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
my ($A_name, $A_id, $B_name, $B_id, $C_name, $C_id, $D_name, $D_id);
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    #my @row = split(/\t/, $_);

    if($_ =~ m/^A$/){next;}
    if($_ =~ m/^B$/){next;}
    if($_ =~ m/^C$/){next;}
    if($_ =~ m/^D$/){next;}

    if($_ =~ m/^A(\d+) (.*)$/){
        $A_name = $2;
        $A_id = $1;
    }
    if($_ =~ m/^B  (\d+) (.*)$/){
        $B_name = $2;
        $B_id = $1;
    }
    if($_ =~ m/^C    (\d+) (.*)$/){
        $C_name = $2;
        $C_id = $1;
    }
    if($_ =~ m/^D      (K\d\d\d\d\d)  (.*)$/){
        $D_name = $2;
        $D_id = $1;
       
        $hash{$A_id}{$B_id}{"ko".$C_id}{$D_id} = $D_name;
        $hash{$A_id}{A_name} = $A_name;
        $hash{$A_id}{$B_id}{B_name} = $B_name;
        $hash{$A_id}{$B_id}{"ko".$C_id}{C_name} = $C_name;
        #$hash{$A_id}{$B_id}{$C_id}{$D_id}{D_name} = $D_name;


#        $hash{$A_id}{level1} = $A_id;
#        $hash{$A_id}{level1_desc} = $A_name;
#        $hash{$A_id}{level2} = $B_id;
#        $hash{$A_id}{level2_desc} = $B_name;
#        $hash{$A_id}{ko} = "ko".$C_id;
#        $hash{$A_id}{ko_desc} = $C_name;
#        $hash{$A_id}{KO} = $D_id;
#        $hash{$A_id}{KO_desc} = $D_name;
#
#        $hash{$D_id}{level1} = $A_id;
#        $hash{$D_id}{level1_desc} = $A_name;
#        $hash{$D_id}{level2} = $B_id;
#        $hash{$D_id}{level2_desc} = $B_name;
#        $hash{$D_id}{ko} = "ko".$C_id;
#        $hash{$D_id}{ko_desc} = $C_name;
#        $hash{$D_id}{KO} = $D_id;
#        $hash{$D_id}{KO_desc} = $D_name;
    }
}
close(IN);
print STDERR Dumper(\%hash);
print STDOUT "level1\tlevel1_desc\tlevel2\tlevel2_desc\tko\tko_desc\tKO\tKO_desc\n";

foreach my $A_id (sort {$a <=> $b} keys %hash){
    #if($C_id eq "B_name"){next;}
    foreach my $B_id (keys %{ $hash{$A_id} }) {
        if($B_id eq "A_name"){next;}
        foreach my $C_id (keys %{ $hash{$A_id}{$B_id} }) {
            if($C_id eq "B_name"){next;}
            foreach my $D_id (keys %{ $hash{$A_id}{$B_id}{$C_id} }) {
                if($D_id eq "C_name"){next;}
                print STDOUT $A_id."\t".$hash{$A_id}{A_name}."\t";
                print STDOUT $B_id."\t".$hash{$A_id}{$B_id}{B_name}."\t";
                print STDOUT $C_id."\t".$hash{$A_id}{$B_id}{$C_id}{C_name}."\t";
                print STDOUT $D_id."\t".$hash{$A_id}{$B_id}{$C_id}{$D_id}."\n";
            }
        }
    }
}

#foreach my $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
#foreach my $D_id (sort { $hash{$a}{level1} <=> $hash{$b}{level1} } keys %hash){
#    print STDOUT $hash{$D_id}{level1}."\t".$hash{$D_id}{level1_desc}."\t";
#    print STDOUT $hash{$D_id}{level2}."\t".$hash{$D_id}{level2_desc}."\t";
#    print STDOUT $hash{$D_id}{ko}."\t".$hash{$D_id}{ko_desc}."\t";
#    print STDOUT $hash{$D_id}{KO}."\t".$hash{$D_id}{KO_desc}."\n";
#}
exit;
