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
--infile_fasta <string>                             : Sequence file.
--infile_replace_species_with_genus <string>        : txt file. contains quirky species name
                                                      that needs to be replace with generic
                                                      species names based on their genus.
--infile_change_species_name_if_genus_diff <string> : txt file contianing species name
                                                      that contains conflicts. If species name
                                                      prefix diff form genus, assign new species
                                                      name based on the genus name.
 
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_remove, $infile_keep, $infile_replace_species_with_genus, $infile_change_species_name_if_genus_diff);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'                              => \$infile_fasta,
   'infile_remove=s'                             => \$infile_remove,
   'infile_keep=s'                               => \$infile_keep,
   'infile_replace_species_with_genus=s'         => \$infile_replace_species_with_genus,
   'infile_change_species_name_if_genus_diff=s'  => \$infile_change_species_name_if_genus_diff,
   'verbose'                                     => \$verbose,
   'help'                                        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash_quirky;
open(IN, "<".$infile_replace_species_with_genus) or die "Can't open $infile_replace_species_with_genus\n";
while(<IN>){
    chomp;
    $hash_quirky{lc($_)} = lc($_);
}
close(IN);

my %hash_replace;
open(IN, "<".$infile_change_species_name_if_genus_diff) or die "Can't open $infile_change_species_name_if_genus_diff\n";
while(<IN>){
    chomp;
    $hash_replace{lc($_)} = lc($_);
}
close(IN);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
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
        
        $k =~ s/<//;
        $k =~ s/>//;
        $p =~ s/<//;
        $p =~ s/>//;
        $c =~ s/<//;
        $c =~ s/>//;
        $o =~ s/<//;
        $o =~ s/>//;
        $f =~ s/<//;
        $f =~ s/>//;
        $g =~ s/<//;
        $g =~ s/>//;
        $s =~ s/<//;
        $s =~ s/>//;
        
        $k2 = lc($row[0]);
        $p2 = lc($row[1]);
        $c2 = lc($row[2]);
        $o2 = lc($row[3]);
        $f2 = lc($row[4]);
        $g2 = lc($row[5]);
        $s2 = lc($row[6]);

        $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
        $parsed_lineage2 = "$k2;$p2;$c2;$o2;$f2;$g2;$s2";
        $k2 =~ s/k__//;
        $k2 =~ s/;$//;
        $p2 =~ s/p__//;
        $p2 =~ s/;$//;
        $c2 =~ s/c__//;
        $c2 =~ s/;$//;
        $o2 =~ s/o__//;
        $o2 =~ s/;$//;
        $f2 =~ s/f__//;
        $f2 =~ s/;$//;
        $g2 =~ s/g__//;
        $g2 =~ s/;g//;
        $s2 =~ s/s__//;
        $s2 =~ s/;$//;
        
        $k2 =~ s/<//;
        $k2 =~ s/>//;
        $p2 =~ s/<//;
        $p2 =~ s/>//;
        $c2 =~ s/<//;
        $c2 =~ s/>//;
        $o2 =~ s/<//;
        $o2 =~ s/>//;
        $f2 =~ s/<//;
        $f2 =~ s/>//;
        $g2 =~ s/<//;
        $g2 =~ s/>//;
        $s2 =~ s/<//;
        $s2 =~ s/>//;

        # Parse entry name that is too long.
        if($p =~ m/Marine Hydrothermal Vent Group/){
            #print STDERR "Found ".$p."\n";
            if($p =~ m/MHVG-1/){
                #print STDERR "Found group 1\n";
                $k = "k__Archaea";
                $p = "p__MHVG-1";
                $c = "c__MHVG-1";
                $o = "o__MHVG-1";
                $f = "f__MHVG-1";
                $g = "g__MHVG-1";
                $s = "s__MHVG-1";
                $parsed_lineage = "k__Archaea;p__MHVG-1;c__MHVGPH-1;o__MHVGPH-1;f__MHVGPH-1;g__MHVGPH-1;s__MHVGPH-1"
            
            }elsif($p =~ m/MHVG-2/){
                $k = "k__Archaea";
                $p = "p__MHVG-2";
                $c = "c__MHVG-2";
                $o = "o__MHVG-2";
                $f = "f__MHVG-2";
                $g = "g__MHVG-2";
                $s = "s__MHVG-2";
                $parsed_lineage = "k__Archaea;p__MHVG-2;c__MHVGPH-2;o__MHVGPH-2;f__MHVGPH-2;g__MHVGPH-2;s__MHVGPH-2"
            
            }else{
                $k = "k__Archaea";
                $p = "p__MHVG";
                $c = "c__MHVG";
                $o = "o__MHVG";
                $f = "f__MHVG";
                $g = "g__MHVG";
                $s = "s__MHVG";
                $parsed_lineage = "k__Archaea;p__MHVG;c__MHVGPH;o__MHVGPH;f__MHVGPH;g__MHVGPH;s__MHVGPH"
            }
        }
        if($lineage =~ m/possible/i){
            next;
        }

        ####
        # First parse species
        # Remove strain name. full length 16S should be able to identify to up to the species level, but certainly not strain level.
        #print STDERR "s__:".$s."\n" if($s =~ m/Bacillus sp\./);
        if($s =~ m/^(s__\S+ \S+)/){
            $s = $1;
			$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
        }elsif($s =~ m/^(s__\S+)/){
            $s = $1;
		    $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
        }else{
            print STDERR "Could not parse species s__\n";
            exit(1);
        }
        # Bacillus sp. Many species are labelled bacillus sp. even if they do not belong to the Bacillus genus. Correct that...
        if($s eq "s__Bacillus sp."){
            if($g ne "g__Bacillus"){
                $s = "s__".$g2." sp.";
                #print STDERR "s__:".$s."\n";
				$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            }
        }
        if($s eq "s__Shewanella sp."){
            if($g ne "g__Shewanella"){
                $s = "s__".$g2." sp.";
                #print STDERR "s__:".$s."\n";
				$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            }   
        } 
        if($s eq "s__clostridium sp."){
            if($g ne "g__Clostridium" || $g ne "g__Clostridium sensu stricto"){
                $s = "s__".$g2." sp.";
                #print STDERR "s__:".$s."\n";
				$parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            }   
        } 


        if($o eq "o__Unknown Order"){
            $o = $p."-".$c."-".$o; 
            $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        
        }elsif($f eq "f__Family XI" || $f eq "f__Family XII" || $f eq "f__Family XI" || $f eq "f__Family XII" || $f eq "f__Family XIII" ||  $f eq "f__Family XIV"){
            
            if($o eq "o__Bacillales"){
                $parsed_lineage = "$k;$p;$c;$o;$f(Bacillales);$g;$s";
            }elsif($o eq "o_Clostridiales"){
                $parsed_lineage = "$k;$p;$c;$o;$f(Clostridiales);$g;$s";
            }

        }elsif($g eq "g__uncultured-Unknown Family" && $p eq "p__Bacteroidetes"){
            $g = "g__uncultured-Unknown Family(Bacteroidetes)";
            $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        
        }elsif($g eq "g__uncultured-Unknown Family" && $p eq "p__Proteobacteria"){
            $g = "g__uncultured-Unknown Family(Proteobacteria)";
            $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            if($o eq "o__Gammaproteobacteria Incertae Sedis" && $f eq "f__Unknown Family"){
                $o = "";
                $parsed_lineage = "$k;$p;$c;$o;$f($p2-$c2-$o2);$g;$s";
            }
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        }elsif($c eq "c__Cyanobacteria" && $o =~ m/o__subsection/i){
            $o = "o__Subsection";
            $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";

            if($s eq "s__Synechococcus sp." && $g eq "g__Synechococcus"){
                #print STDERR "Found subsection!\n";
                $s = "s__Synechococcus sp.(Synechococcus)";
                $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
             }elsif($s eq "s__Synechococcus sp." && $g eq "g__Phormidium"){
                $s = "s__Synechococcus sp.(Phormidium)";
                $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
             }elsif($s eq "s__Calothrix sp. PCC 7103" && $g eq "g__Calothrix"){
                $s = "s__Calothrix sp. PCC 7103(Calothrix)";
                $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
                #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
             }
        }elsif($g eq "g__Bacillus" && $s eq "s__Synechocystis sp. PCC 6803"){
            $s = "s__Synechocystis sp. PCC 6803-bacillus";
            $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        }elsif($f eq "f__Unknown Family"){
            $parsed_lineage = "$k;$p;$c;$o;$f($p2-$c2-$o2);$g;$s";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        }elsif($f eq "f__Holophagaceae" && $g eq "g__marine group"){
            $parsed_lineage = "$k;$p;$c;$o;$f;$g($f2);$s";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        
        }elsif($g =~ m/"g__Clostridium sensu stricto"/){
			print STDERR "FOUND CLOSTRIDIUM SENSU STRICTO!\n";
            $g = "g__Clostridium";
            $parsed_lineage = "$k;$p;$c;$o;$f;$g($f2);$s";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";  
        }
		elsif($s eq "s__Archangium disciforme"){ 
            $parsed_lineage = "$k;$p;$c;$o;$f;$g;$s($g2)";
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        }else{
            #print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
        }

        if($g eq "g__unidentified"){
            $parsed_lineage = "$k;$p;$c;$o;$f;$g($f2)";
        }

        # After all modifications, loop through both hash and make final corrections. Do an exception for s__Escherichia coli as it is classified under the g__Shigella-Escheria group.
        my @lineage = split(/;/, $parsed_lineage);

        # E.coli
        if($lineage[6] eq "s__Escherichia coli" && $lineage[5] eq "g__Escherichia-Shigella"){
            # Basically leave as is.
            $parsed_lineage = "$lineage[0];$lineage[1];$lineage[2];$lineage[3];$lineage[4];$lineage[5];$lineage[6]";
        
        }else{
            
            # Then the rest.
            if(exists $hash_quirky{lc($lineage[6])}){
                #print STDERR "found hash_quirky\n";
                my @lineage = split(/;/, $parsed_lineage);
                my $curr_species = $lineage[6];
                my $curr_genus = $lineage[5];
                $curr_genus =~ s/g__//;
                if(lc($curr_species) eq lc($hash_quirky{lc($lineage[6])})){
                    
                    $s = "s__".ucfirst($curr_genus);
                    $parsed_lineage = "$lineage[0];$lineage[1];$lineage[2];$lineage[3];$lineage[4];$lineage[5];$s";
                }
            }
            
            # After all modifications, loop through both hash and make final corrections.
            if(exists $hash_replace{lc($lineage[6])}){
                #print STDERR "found hash_replace: ".$hash_replace{lc($lineage[6])}."\n";
                my $curr_species = $lineage[6];
                my $curr_genus = $lineage[5];
                my $curr_genus_parsed = $curr_genus;
                $curr_genus_parsed =~ s/g__//;
                my $curr_species_parsed = $curr_species;
                $curr_species_parsed =~ s/s__//;
                $curr_species_parsed =~ s/ sp\.//;
                #print STDERR "curr_species: ".$curr_species_parsed."\n";
                #print STDERR "curr_genus: ".$curr_genus_parsed."\n";
    
                if(lc($curr_genus_parsed) ne lc($curr_species_parsed)){
                    $s = "s__".ucfirst($curr_genus)." sp.";
                    $parsed_lineage = "$lineage[0];$lineage[1];$lineage[2];$lineage[3];$lineage[4];$lineage[5];$s";
                }
            }
    
            @lineage = split(/;/, $parsed_lineage);
            for my $lineage (@lineage){
                my $curr_species = $lineage[6];
                $curr_species =~ s/s__//;
                $curr_species =~ s/G__//;
                $curr_species = "s__".ucfirst($curr_species);
                if($curr_species eq "s__Escherichia-shigella sp."){
                    $curr_species = "s__Escherichia-Shigella sp.";
                }
                $parsed_lineage = "$lineage[0];$lineage[1];$lineage[2];$lineage[3];$lineage[4];$lineage[5];$curr_species";
            }
    
            # One last check for one offending sensu stricto and clostridium in general:
            # Another check for species labelled something like: s__Bacterium T12kex. Assign to genus.
            @lineage = split(/;/, $parsed_lineage);
            for my $lineage (@lineage){
                my $curr_family = $lineage[4];
                my $curr_genus = $lineage[5];
                my $curr_species = $lineage[6];
                if($curr_family eq "f__Clostridiaceae 4" && $curr_genus eq "g__Clostridium sensu stricto"){
                   $curr_genus = "g__Clostridium sensu stricto(Clostridiaceae 4)"; 
                   $curr_species = "s__Clostridium sensu stricto(Clostridiaceae 4)"; 
                }elsif($curr_family eq "f__Clostridiaceae 2" && $curr_genus eq "g__Clostridium sensu stricto"){
                   $curr_genus = "g__Clostridium sensu stricto(Clostridiaceae 2)"; 
                   $curr_species = "s__Clostridium sensu stricto(Clostridiaceae 2)"; 
                }

                # s__Bacterium
                if($curr_species =~ m/s__Bacterium.*/i){
                   $curr_species = $curr_genus;
                   $curr_genus =~ s/g__/s__/;
                   $curr_genus = $curr_genus." sp.";
                }
                
                # uncultured in lc
                if($curr_species =~ m/s__Uncultured-(.*)/){
                    my $suffix = ucfirst($1);
                    $curr_species = "s__uncultured-$suffix";
                }
                $parsed_lineage = "$lineage[0];$lineage[1];$lineage[2];$lineage[3];$curr_family;$curr_genus;$curr_species";
            }
        }
        print STDOUT ">".$access." ".$parsed_lineage."\n".$curr->seq."\n";
    }
}
#print STDERR Dumper(\%hash_correct2);
