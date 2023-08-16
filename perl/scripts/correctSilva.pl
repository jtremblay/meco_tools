#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
correctSilva.pl

PURPOSE:

INPUT:
--infile_fasta <string>   : Sequence file.
--infile_remove <string>  : txt file of offending values.
--infile_keep <string>    : txt file of offending values.
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_remove, $infile_keep, $infile_correct);
my $verbose = 0;

GetOptions(
   'infile_fasta=s' => \$infile_fasta,
   'infile_remove=s' => \$infile_remove,
   'infile_keep=s' => \$infile_keep,
   #'infile_correct=s' => \$infile_correct,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash_remove;
open(IN, "<".$infile_remove) or die "Can't open $infile_remove\n";
while(<IN>){
    chomp;
    $hash_remove{$_} = 0;
}
close(IN);
#print STDERR Dumper(\%hash);
#exit;

my %hash_keep;
open(IN, "<".$infile_keep) or die "Can't open $infile_keep\n";
while(<IN>){
    chomp;
    $hash_keep{$_} = 0;
}
close(IN);
#print STDERR Dumper(\%hash);
#exit;

my %hash_correct;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    if($header =~ m/Bacteria|Archaea/){
        my($access) = $header =~ m/^>(\S+) /;
        my($lineage) = $header =~ m/^>\S+ (.*)/;
        my @row = split(/;/, $lineage);
    
        my ($parsed_lineage, $k, $p, $c, $o, $f, $g, $s);
        
        $k = lc($row[0]);
        $p = lc($row[1]);
        $c = lc($row[2]);
        $o = lc($row[3]);
        $f = lc($row[4]);
        $g = lc($row[5]);
        $s = lc($row[6]);

        $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
        #$s =~ s/s__//;
        $s =~ s/;$//;

        if(exists $hash_keep{$s}){
            #print STDERR $s."\n";
            $hash_correct{$s}{$parsed_lineage}++;
        }
    }
}
#print STDERR Dumper(\%hash_keep);
print STDERR Dumper(\%hash_correct);
#exit;

#loop through $hash{$s}
my %hash_correct2;
my %hash_correct2rev;
foreach my $key (keys %hash_correct) {
    my $max = 0;
    my $ref_key;
    my $ref_key2;
    foreach my $key2 (keys %{ $hash_correct{$key} }) {
        my $value = $hash_correct{$key}{$key2};
        if($value > $max){
            $max = $value;
            $ref_key = $key;
            $ref_key2 = $key2;
        }
    }
    $hash_correct2{$ref_key} = $ref_key2;
    $hash_correct2rev{$ref_key2} = $ref_key;
}

$ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    if($header =~ m/Bacteria|Archaea/){
        my($access) = $header =~ m/^>(\S+) /;
        my($lineage) = $header =~ m/^>\S+ (.*)/;
        my @row = split(/;/, $lineage);
    
        my ($parsed_lineage, $parsed_lineage2, $k, $p, $c, $o, $f, $g, $s,  $k2, $p2, $c2, $o2, $f2, $g2, $s2);
        
        $k = $row[0];
        $p = $row[1];
        $c = $row[2];
        $o = $row[3];
        $f = $row[4];
        $g = $row[5];
        $s = $row[6];
        
        $k2 = lc($row[0]);
        $p2 = lc($row[1]);
        $c2 = lc($row[2]);
        $o2 = lc($row[3]);
        $f2 = lc($row[4]);
        $g2 = lc($row[5]);
        $s2 = lc($row[6]);

        $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
        $parsed_lineage2 = "$k2;$p2;$c2;$o2;$f2;$g2;$s2";
        #$s =~ s/s__//;
        $s2 =~ s/;$//;

        if(!exists $hash_remove{$s2} && 
           !exists $hash_remove{$g2} && 
           !exists $hash_remove{$f2} && 
           !exists $hash_remove{$o2} && 
           !exists $hash_remove{$c2} && 
           !exists $hash_remove{$p2} && 
            !exists $hash_keep{$s2}){
            print STDOUT $header."\n".$curr->seq."\n";
        }
        
        if(exists $hash_keep{$s2}){
            if(exists $hash_correct2rev{$parsed_lineage2}){
                print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
                print STDERR ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
            }
        }
    }
}
print STDERR Dumper(\%hash_correct2);
