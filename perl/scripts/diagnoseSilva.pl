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
--infile_fasta <string>     : Sequence file.
--infile_conflicts <string> : txt file of offending values.
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_conflicts);
my $verbose = 0;

GetOptions(
   'infile_fasta=s' => \$infile_fasta,
   'infile_conflicts=s' => \$infile_conflicts,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile_conflicts) or die "Can't open $infile_conflicts\n";
while(<IN>){
    chomp;
    my $offending;
    if($_ =~ m/^s__(.*) species/){
        $offending = $1;
        $offending =~ s/^\s+//;
        $offending =~ s/\s+$//;
        $hash{$offending} = 0;
    }
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
        $s =~ s/s__//;
        $s =~ s/;$//;

        if(exists $hash{$s}){
            #print STDERR $s."\n";
            $hash_correct{$s}{$parsed_lineage}++;
        }
    }
}
print STDERR Dumper(\%hash_correct);
