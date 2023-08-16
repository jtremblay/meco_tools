#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
addTaxonomyToFasta.pl

PURPOSE:

INPUT:
--infile_fasta <string>    : Sequence file in fasta format.
--infile_taxonomy <string> : taxonomy file (tsv).
				
OUTPUT:
STDOUT <string>            : fasta file where header contains taxonomy.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxonomy, $infile_fasta);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'    => \$infile_fasta,
   'infile_taxonomy=s' => \$infile_taxonomy,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){
    chomp;
    next if($. == 1);
    my @row = split(/\t/, $_);
    my $k = $row[1];
    my $p = $row[2];
    my $c = $row[3];
    my $o = $row[4];
    my $f = $row[5];
    my $g = $row[6];
    my $s = $row[7];

    # Flag offending lineages;
    if($g eq "Thiomonas"){
        $f = "Comamonadaceae";
    }
    
    # Then correct specific offeding lineages
    if($f eq "Vannellidae"){
        $p = "Amoebozoa"; $c = "Discosea"; $o = "Vannellida";
    }
    if($g eq "Bachelotia"){
        $p = "Ochrophyta"; $c = "Phaeophyceae"; $o =  "Scytothamnales"; $f = "Bachelotiaceae";
    }
    if($f eq "Himatismenida"){
        $p = "Amoebozoa"; $c = "Discosea";
    }
    if($c eq "Xanthophyceae" || $c eq "Eustigmatophyceae" || $c eq "Raphidophyceae" || $c eq "Phaeophyceae" || $c eq "Raphidophyceae" || $c eq "Dictyochophyceae" || $c eq "Chrysophyceae" || $c eq "Synurophyceae"){
        $p = "Ochrophyta"
    }
    if($c eq "Oomycetes"){
        $p = "Oomycota";
    }
    if($c eq "Dinophyceae"){
        $p = "Myzozoa"
    }
    


    my $lineage;
    if($o eq "NULL" && $f eq "NULL" && $g eq "NULL"){
       $o = $c."-undef";
       $f = $c."-undef";     
       $g = $c."-undef";     
    }elsif($o eq "NULL" && $f eq "NULL"){
       $o = $c."-undef";
       $f = $c."-undef";
         
    }elsif($o eq "NULL"){
       $o = $c."-undef";
    }elsif($f eq "NULL" && $g eq "NULL"){
       $f = $o."-undef";
       $g = $o."-undef";
    }elsif($g eq "NULL"){
       $g = $f."-undef";
    }elsif($f eq "NULL"){
       $f = $o."-undef";
    }

    my @s = split(/\s/, $s);
    if(@s > 1){
        $s = $s[0]." ".$s[1];
    }

    if($p eq "NULL" && $c eq "Oomycetes"){
        $k = "Fungi";
        $p = "Oomycota";
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }elsif($p eq "NULL" && $c eq "Hyphochytriomycetes" && $f eq "Hyphochytriaceae"){
        $p = "Hyphochytriomycota";
        $o = "Hyphochytriales";
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }elsif($p eq "NULL" && $f eq "Hyphochytriomycetes" && $f eq "Rhizidiomycetaceae"){
        $p = "Hyphochytriomycota";
        $o = "Hyphochytriales";
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }elsif($p eq "NULL" && $c eq "Hyphochytriomycetes"){
        $p = "Hyphochytriomycota";
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }elsif($p eq "NULL" && $c eq "Chrysophyceae"){
        $k = "Eukaryota";
        $p = "Ochrophyta";
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }elsif($p eq "NULL" && $c eq "Dinophyceae"){
        $k = "Eukaryota";
        $p = "Myzozoa";
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }elsif($p eq "NULL" && $c ne "NULL"){
        $lineage = "k__".$k.";"."p__undef-PH-".$k.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s; 
    }elsif($p eq "NULL"){
        #$lineage = "k__Bacteria"; 
        $lineage = "k__".$k;
        #$lineage = "k__".$k.";"."p__undef-PH-".$k.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s; 
    }elsif($c eq "NULL"){
        $lineage = "k__".$k.";"."p__".$p;
        #$lineage = "k__$k;p__".$p.";"."c__undef-CL-".$p.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;;
    }elsif($o eq "NULL" && $f eq "NULL" && $g eq "NULL"){
        $lineage = "k__$k;p__".$p.";c__".$c.";o__".$o;
    }else{
        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    }
    
    $hash{$row[0]} = $lineage;
}

#print STDERR Dumper(\%hash);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my ($header) = $curr->header =~ m/^>(\S+)/;
    if(exists $hash{$header}){
        print STDOUT ">".$header." ".$hash{$header}."\n".$curr->seq."\n";        
    }#else{
    #   print STDERR "Did not find $header seq\n";
    #}
}
exit;


