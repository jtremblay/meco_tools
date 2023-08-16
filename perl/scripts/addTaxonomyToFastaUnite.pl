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
    my @row2 = split(/;/, $row[1]);
    my $k = $row2[0];
    $k =~ s/k__//;
    my $p = $row2[1];
    $p =~ s/p__//;
    my $c = $row2[2];
    $c =~ s/c__//;
    my $o = $row2[3];
    $o =~ s/o__//;
    my $f = $row2[4];
    $f =~ s/f__//;
    my $g = $row2[5];
    $g =~ s/g__//;
    my $s = $row2[6];
    $s =~ s/s__//;

    if($g eq "Dufourea"){$g = $g."(".$k.")"; }
    if($g eq "Rozella"){$g = $g."(".$k.")"; }
    if($g eq "Petrophila"){$g = $g."(".$k.")"; }
    if($g eq "Mertensia"){$g = $g."(".$k.")"; }
    if($g eq "Helicotylenchus"){ $o = "Tylenchida"; }
    if($g eq "Diplodiscus"){$g = $g."(".$k.")"; }
    if($g eq "Chondrilla"){$g = $g."(".$k.")"; }
    if($g eq "Rhytisma"){$g = $g."(".$k.")"; }
    if($g eq "Pilophorus"){$g = $g."(".$k.")"; }
    if($g eq "Leucoptera"){$g = $g."(".$k.")"; }
    if($g eq "Hypnum"){ $f = "Hypnaceae"; }
    if($g eq "Hymenolepis"){$g = $g."(".$k.")"; }
    if($g eq "Setaria"){$g = $g."(".$k.")"; }
    if($g eq "Siphonaria"){$g = $g."(".$k.")"; }
    if($g eq "Leptobryum"){ $o = "Bryales";}
    if($g eq "Tetracladium"){$g = $g."(".$k.")"; }
    if($g eq "Acanthella"){$g = $g."(".$k.")"; }
    if($g eq "Opercularia"){$g = $g."(".$k.")"; }
    if($g eq "Rhabdias"){$o = "Rhabditida"; $f = "Rhabdiasidae";}
    if($g eq "Lacrymaria"){$g = $g."(".$k.")"; }
    if($g eq "Sciuro-hypnum"){$o = "Hypnales"; }
    if($g eq "Stokesia"){$g = $g."(".$k.")"; }
    if($g eq "Alaria"){$g = $g."(".$k.")"; }
    if($g eq "Gomphus"){$g = $g."(".$k.")"; }
    if($g eq "Bonamia"){$g = $g."(".$k.")"; }
    if($g eq "Pyramidula"){$g = $g."(".$k.")"; }
    if($g eq "Discula"){$g = $g."(".$k.")"; }
    if($g eq "Rhyacodrilus"){$o = "Haplotaxida";}
    if($g eq "Sarcophyton"){$g = $g."(".$k.")"; }
    if($g eq "Brachythecium"){$o = "Hypnales";}
    if($g eq "Herzogiella"){$o = "Hypnales";}
    if($g eq "Eisenia"){$g = $g."(".$k.")"; }
    if($g eq "Anomodon"){$f = "Thuidiaceae";}
    if($g eq "Stenostomum"){$g = $g."(".$k.")"; }
    if($g eq "Stagnicola"){$g = $g."(".$k.")"; }
    if($g eq "Crambe"){$g = $g."(".$k.")"; }
    if($g eq "Uronema"){$g = $g."(".$k.")"; }
    if($f eq "Hoplolaimidae"){$o = "Rhabditida"; }
    if($f eq "Brachytheciaceae"){$o = "Hypnales"; }
    if($f eq "Tubificidae"){$o = "Tubificida"; }
    if($g eq "unidentified"){$g = $g."(".$f.")";}
    if($f eq "Enchytraeidae"){$o = "Haplotaxida"; }
    if($g eq "Mesenchytraeus"){ $o = "Haplotaxida"; }

    my  $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
    
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

#
    # Flag offending lineages;
#    if($g eq "Thiomonas"){
#        $f = "Comamonadaceae";
#    }
#    
#    # Then correct specific offeding lineages
#    if($f eq "Vannellidae"){
#        $p = "Amoebozoa"; $c = "Discosea"; $o = "Vannellida";
#    }
#    if($g eq "Bachelotia"){
#        $p = "Ochrophyta"; $c = "Phaeophyceae"; $o =  "Scytothamnales"; $f = "Bachelotiaceae";
#    }
#    if($f eq "Himatismenida"){
#        $p = "Amoebozoa"; $c = "Discosea";
#    }
#    if($c eq "Xanthophyceae" || $c eq "Eustigmatophyceae" || $c eq "Raphidophyceae" || $c eq "Phaeophyceae" || $c eq "Raphidophyceae" || $c eq "Dictyochophyceae" || $c eq "Chrysophyceae" || $c eq "Synurophyceae"){
#        $p = "Ochrophyta"
#    }
#    if($c eq "Oomycetes"){
#        $p = "Oomycota";
#    }
#    if($c eq "Dinophyceae"){
#        $p = "Myzozoa"
#    }
#    
#    my $lineage;
#    if($o eq "NULL" && $f eq "NULL" && $g eq "NULL"){
#       $o = $c."-undef";
#       $f = $c."-undef";     
#       $g = $c."-undef";     
#    }elsif($o eq "NULL" && $f eq "NULL"){
#       $o = $c."-undef";
#       $f = $c."-undef";
#         
#    }elsif($o eq "NULL"){
#       $o = $c."-undef";
#    }elsif($f eq "NULL" && $g eq "NULL"){
#       $f = $o."-undef";
#       $g = $o."-undef";
#    }elsif($g eq "NULL"){
#       $g = $f."-undef";
#    }elsif($f eq "NULL"){
#       $f = $o."-undef";
#    }
#
#    my @s = split(/\s/, $s);
#    if(@s > 1){
#        $s = $s[0]." ".$s[1];
#    }
#
#    if($p eq "NULL" && $c eq "Oomycetes"){
#        $k = "Fungi";
#        $p = "Oomycota";
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }elsif($p eq "NULL" && $c eq "Hyphochytriomycetes" && $f eq "Hyphochytriaceae"){
#        $p = "Hyphochytriomycota";
#        $o = "Hyphochytriales";
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }elsif($p eq "NULL" && $f eq "Hyphochytriomycetes" && $f eq "Rhizidiomycetaceae"){
#        $p = "Hyphochytriomycota";
#        $o = "Hyphochytriales";
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }elsif($p eq "NULL" && $c eq "Hyphochytriomycetes"){
#        $p = "Hyphochytriomycota";
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }elsif($p eq "NULL" && $c eq "Chrysophyceae"){
#        $k = "Eukaryota";
#        $p = "Ochrophyta";
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }elsif($p eq "NULL" && $c eq "Dinophyceae"){
#        $k = "Eukaryota";
#        $p = "Myzozoa";
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }elsif($p eq "NULL" && $c ne "NULL"){
#        $lineage = "k__".$k.";"."p__undef-PH-".$k.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s; 
#    }elsif($p eq "NULL"){
#        #$lineage = "k__Bacteria"; 
#        $lineage = "k__".$k;
#        #$lineage = "k__".$k.";"."p__undef-PH-".$k.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s; 
#    }elsif($c eq "NULL"){
#        $lineage = "k__".$k.";"."p__".$p;
#        #$lineage = "k__$k;p__".$p.";"."c__undef-CL-".$p.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;;
#    }elsif($o eq "NULL" && $f eq "NULL" && $g eq "NULL"){
#        $lineage = "k__$k;p__".$p.";c__".$c.";o__".$o;
#    }else{
#        $lineage = "k__".$k.";"."p__".$p.";"."c__".$c.";"."o__".$o.";"."f__".$f.";"."g__".$g.";"."s__".$s;
#    }
