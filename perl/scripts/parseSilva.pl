#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseSilva.pl

PURPOSE:

INPUT:
--infile_taxa <string>   : tsv file
--infile_access <string> : tsv file
--infile_fasta <string>  : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxa, $infile_access, $infile_fasta);
my $verbose = 0;

GetOptions(
   'infile_taxa=s'   => \$infile_taxa,
   'infile_access=s' => \$infile_access,
   'infile_fasta=s'  => \$infile_fasta,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

my %hash_access;
my %hash_id;
open(IN, "<".$infile_access) or die "Can't open $infile_access\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $access = $row[0];
   my $id = $row[1];
   $hash_access{$access}{id} = $id
}
close(IN);

open(IN, "<".$infile_taxa) or die "Can't open $infile_taxa\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $taxa = $row[0];
   my $id = $row[1];
   my $domain = $row[2];
   $hash_id{$id}{taxa} = $taxa;
   $hash_id{$id}{domain} = $domain;
}
close(IN);

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    if($header =~ m/Bacteria|Archaea/){
        my($access) = $header =~ m/^>(\S+) /;
        my($orig_lineage) = $header =~ m/^>\S+ (.*)/;
        my @orig_row = split(/;/, $orig_lineage);
    
        if(exists $hash_access{$access}){
            my $id = $hash_access{$access}{id};
            if(exists $hash_id{$id}){
                my $domain = $hash_id{$id}{domain}; 
                my $taxa = $hash_id{$id}{taxa}; 
            
                my @row = split(/;/, $taxa);
                my ($parsed_lineage, $k, $p, $c, $o, $f, $g, $g2, $s);
               
                # Correct for inconstencies. 
                if($domain eq "family" && $row[4] == "") {
                    print STDERR "Found inconsistent family domain : ".$taxa."\n";
                    $domain = "order";
                }
            
                if($domain eq "genus"){
                    $k = "k__".$row[0];
                    $p = "p__".$row[1];
                    $c = "c__".$row[2];
                    $o = "o__".$row[3];
                    $f = "f__".$row[4];
                    $g = "g__".$row[5];
                    $g2 = "$row[5]GE";
                    $s = "s__".$row[5]."GE";
                    if($g =~ m/uncultured/i){
                        $g = "g__uncultured-$row[4]";
                    }
                    #$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
                if($domain eq "family"){
                    $k = "k__".$row[0];
                    $p = "p__".$row[1];
                    $c = "c__".$row[2];
                    $o = "o__".$row[3];
                    $f = "f__".$row[4];
                    $g = "g__".$row[4]."FA";
                    $s = "s__".$row[4]."FA";
                    #$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
                if($domain eq "order"){
                    $k = "k__".$row[0];
                    $p = "p__".$row[1];
                    $c = "c__".$row[2];
                    $o = "o__".$row[3];
                    $f = "f__".$row[3]."OR";
                    $g = "g__".$row[3]."OR";
                    $s = "s__".$row[3]."OR";
                    #$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
                if($domain eq "class"){
                    $k = "k__".$row[0];
                    $p = "p__".$row[1];
                    $c = "c__".$row[2];
                    $o = "o__".$row[2]."CL";
                    $f = "f__".$row[2]."CL";
                    $g = "g__".$row[2]."CL";
                    $s = "s__".$row[2]."CL";
                    #$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
                if($domain eq "phylum"){
                    $k = "k__".$row[0];
                    $p = "p__".$row[1];
                    $c = "c__".$row[1]."PH";
                    $o = "o__".$row[1]."PH";
                    $f = "f__".$row[1]."PH";
                    $g = "g__".$row[1]."PH";
                    $s = "s__".$row[1]."PH";
                    #$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
                if($domain eq "domain"){
                    $k = "k__".$row[0];
                    $p = "c__".$row[0]."KI";
                    $c = "c__".$row[0]."KI";
                    $o = "o__".$row[0]."KI";
                    $f = "f__".$row[0]."KI";
                    $g = "g__".$row[0]."KI";
                    $s = "s__".$row[0]."KI";
                    #$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
        
                # Try to match a species as well.
                if(@orig_row == 7){
                    $s = $orig_row[6];
                    if($s =~ m/uncultured/i){
                        $s = "s__uncultured-$g2";
                    }else{
                        $s = "s__".$orig_row[6];
                    }
                    #:$parsed_lineage = $parsed_lineage.";".$s;
                    $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }else{
                    $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                }
                $hash_access{$access}{lineage} = $parsed_lineage;
                print STDOUT ">$access $parsed_lineage\n".$curr->seq."\n";
            }
        }
    }
}

#print STDERR Dumper(\%hash_access);

