#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use List::Util qw(max);
use List::UtilsBy qw(max_by);

my $usage=<<'ENDHERE';
NAME:
generateConsensusTaxonomy.pl

PURPOSE:
Implementation of the methodology described in CAT paper - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1817-x
The implementation is supplemented (optional) with integration of blastn results of contigs vs ncbi complete genome only with CAT workflow (not BAT).
Basically, once CAT has completed, go through all contigs taxonomy and replace lineage with the one obtained by blastn vs ncbi complete genome
if evalue > 1e-80 and qlen / slen > 0.5. Can be parametrized.

INPUT:
--infile_taxonomy_blastp <string>     : diamond blastp results against nr
--infile_taxonomy_blastn <string>     : blastn results against ncbi complete genomes (optional). This file has to be generated using blastn --outfmt "6 qlen slen".
                                        File is expected to not contain header. 
--type <string>                       : CAT (default) or BAT.
--r <int>                             : r parameter from the CAT paper. Default=10 for CAT and 5 for BAT. the lower = more permissive. higher = more conservative.
--f <float>                           : f parameter from the CAT paper. Default=0.5
--ratio <float>                       : if --infile_taxonomy_blastn is provided. Will replace taxonomic lineage obtained by CAT if contig hit align_len/qlen > 0.5 and < --evalue
--evalue <float>                      : if --infile_taxonomy_blastn is provided. Will replace taxonomic lineage obtained by CAT if contig hit align_len/qlen > 0.5 and < --evalue


OUTPUT:
STDOUT <string>                       : consensus taxonomy with score for each contig
--out_stats <string>                  : comprehensive statistics spreadsheet. 
--out_with_blastn_correction <string> : outfile with consensus taxonomy corrected with blastn results. Only specify in combination with --infile_taxonomy_blastn.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_taxonomy_blastp, $infile_taxonomy_blastn, $out_stats, $type, $r, $f, $ratio, $evalue, $out_with_blastn_correction);
my $verbose = 0;

GetOptions(
   'infile_taxonomy_blastp=s'     => \$infile_taxonomy_blastp,
   'infile_taxonomy_blastn=s'     => \$infile_taxonomy_blastn,
   'out_with_blastn_correction=s' => \$out_with_blastn_correction,
   'out_stats=s'                  => \$out_stats,
   'type=s'                       => \$type,
   'r=i'                          => \$r,
   'f=f'                          => \$f,
   'ratio=f'                      => \$ratio,
   'evalue=f'                     => \$evalue,
   'verbose'                      => \$verbose,
   'help'                         => \$help
);
if ($help) { print $usage; exit; }

## MAIN

open(STATS, ">".$out_stats) or die "Can't open $out_stats\n";
$type = "CAT" unless($type);
if($type ne "CAT" and $type ne "BAT"){
    die "--type has to be 'CAT' or 'BAT'\n";
}
if(!defined($r) && $type eq "BAT"){
    $r = 5;
}elsif(!defined($r) && $type eq "CAT"){
    $r = 10;
}
$f = 0.5 unless($f);
$ratio = 0.5 unless($ratio);
$evalue = 1e-80 unless($evalue);

if($out_with_blastn_correction){
    open(OUT, ">".$out_with_blastn_correction) or die "Can't open $out_with_blastn_correction for writing\n";
}

print STDERR "Running CAT/BAT with --type $type r=$r; f=$f\n";

my %hash;
my %hash_lca;
my $last_gene_id = "";
my $last_contig_id = "";
my $last_bitscore = 0;
my $top_bitscore;
my $curr_distance;
#my (@kingdom, @phylum, @class, @order, @family, @genus, @species);
# Loop through taxonomy file. Consider hits ranging within 50% of top hit bit-score. Hits are in order of bitscores. So
# assume first hit has the highest bitscore.
open(IN, "<".$infile_taxonomy_blastp) or die "Can't open $infile_taxonomy_blastp\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    next if($. == 1);
    my $curr_contig_id = $row[0];
    my $curr_gene_id = $row[1];
    my $kingdom = $row[3];
    my $phylum  = $row[4];
    my $class   = $row[5];
    my $order   = $row[6];
    my $family  = $row[7];
    my $genus   = $row[8];
    my $species = $row[9];
    my $bitscore = $row[12];
    if($type eq "BAT"){
        $curr_contig_id = "";
    }

    if($curr_gene_id eq $last_gene_id){
        #continue current gene...    
        $curr_distance = ($last_bitscore/$top_bitscore) * 100;
        if($curr_distance > $r){ #register for vote.
            $hash{kingdom}{$kingdom}++;
            $hash{phylum}{$phylum}++;
            $hash{class}{$class}++;
            $hash{order}{$order}++;
            $hash{family}{$family}++;
            $hash{genus}{$genus}++;
            $hash{species}{$species}++;
         }else{
            print STDERR "Reject: $curr_distance lower than 10%\n";
         }
         
         if(eof){# if eof, conclude here as well to consider the last gene.
             #print STDERR "gene_id:$last_gene_id\n";
             #print STDERR Dumper(\%hash);
            foreach my $key (keys %hash) {
                my $most_abun = 0;
                my $unanimous_taxon = "undef";
                
                foreach my $key2 (keys %{ $hash{$key} }) {
                    # here could implement a majority vote. 
                    #if($hash{$key}{$key2} > $most_abun){
                    #    $most_abun = $hash{$key}{$key2};
                    #    $most_abun_taxon = $key2;
                    #}
                    if($hash{$key}{$key2} == 10){
                        $unanimous_taxon = $key2;
                    }
                }
                if($key eq "kingdom"){ $hash_lca{$last_contig_id}{$last_gene_id}{kingdom}{taxon} = $unanimous_taxon; $hash_lca{$last_contig_id}{$last_gene_id}{kingdom}{bitscore} = $top_bitscore; }
                if($key eq "phylum"){  $hash_lca{$last_contig_id}{$last_gene_id}{phylum}{taxon} = $unanimous_taxon;  $hash_lca{$last_contig_id}{$last_gene_id}{phylum}{bitscore} = $top_bitscore; }
                if($key eq "class"){   $hash_lca{$last_contig_id}{$last_gene_id}{class}{taxon} = $unanimous_taxon;   $hash_lca{$last_contig_id}{$last_gene_id}{class}{bitscore} = $top_bitscore; }
                if($key eq "order"){   $hash_lca{$last_contig_id}{$last_gene_id}{order}{taxon} = $unanimous_taxon;   $hash_lca{$last_contig_id}{$last_gene_id}{order}{bitscore} = $top_bitscore;}
                if($key eq "family"){  $hash_lca{$last_contig_id}{$last_gene_id}{family}{taxon} = $unanimous_taxon;  $hash_lca{$last_contig_id}{$last_gene_id}{family}{bitscore} = $top_bitscore; }
                if($key eq "genus"){   $hash_lca{$last_contig_id}{$last_gene_id}{genus}{taxon} = $unanimous_taxon;   $hash_lca{$last_contig_id}{$last_gene_id}{genus}{bitscore} = $top_bitscore;  }
                if($key eq "species"){ $hash_lca{$last_contig_id}{$last_gene_id}{species}{taxon} = $unanimous_taxon; $hash_lca{$last_contig_id}{$last_gene_id}{species}{bitscore} = $top_bitscore; }
            }
         }

    }else{
        #Complete/conclude current gene...
        #Find LCA (Lowest Common Ancestor). If no common ancestor or all equal hits
        if($. != 2){
            #print STDERR "gene_id:$last_gene_id\n";
            #print STDERR Dumper(\%hash);
            
            foreach my $key (keys %hash) {
                my $most_abun = 0;
                my $unanimous_taxon = "undef";
                my $final_top_bitscore = 0;
                
                foreach my $key2 (keys %{ $hash{$key} }) {
                    # here could implement a more lenient majority vote...? 
                    #if($hash{$key}{$key2} > $most_abun){
                    #    $most_abun = $hash{$key}{$key2};
                    #    $most_abun_taxon = $key2;
                    #}
                    if($hash{$key}{$key2} == 10){
                        $unanimous_taxon = $key2;
                        $final_top_bitscore = $top_bitscore;
                    }
                }
                if($key eq "kingdom"){ $hash_lca{$last_contig_id}{$last_gene_id}{kingdom}{taxon} = $unanimous_taxon; $hash_lca{$last_contig_id}{$last_gene_id}{kingdom}{bitscore} = $final_top_bitscore; }
                if($key eq "phylum"){  $hash_lca{$last_contig_id}{$last_gene_id}{phylum}{taxon} = $unanimous_taxon;  $hash_lca{$last_contig_id}{$last_gene_id}{phylum}{bitscore} = $final_top_bitscore; }
                if($key eq "class"){   $hash_lca{$last_contig_id}{$last_gene_id}{class}{taxon} = $unanimous_taxon;   $hash_lca{$last_contig_id}{$last_gene_id}{class}{bitscore} = $final_top_bitscore; }
                if($key eq "order"){   $hash_lca{$last_contig_id}{$last_gene_id}{order}{taxon} = $unanimous_taxon;   $hash_lca{$last_contig_id}{$last_gene_id}{order}{bitscore} = $final_top_bitscore;}
                if($key eq "family"){  $hash_lca{$last_contig_id}{$last_gene_id}{family}{taxon} = $unanimous_taxon;  $hash_lca{$last_contig_id}{$last_gene_id}{family}{bitscore} = $final_top_bitscore; }
                if($key eq "genus"){   $hash_lca{$last_contig_id}{$last_gene_id}{genus}{taxon} = $unanimous_taxon;   $hash_lca{$last_contig_id}{$last_gene_id}{genus}{bitscore} = $final_top_bitscore;  }
                if($key eq "species"){ $hash_lca{$last_contig_id}{$last_gene_id}{species}{taxon} = $unanimous_taxon; $hash_lca{$last_contig_id}{$last_gene_id}{species}{bitscore} = $final_top_bitscore; }
            }
        } 

        #And start new gene... assume first hit = top bitscore.
        #set contig id
        if($curr_contig_id ne $last_contig_id){
            $last_contig_id = $curr_contig_id;
        }
        $last_gene_id = $curr_gene_id;
        $last_bitscore = $bitscore;
        $top_bitscore = $last_bitscore;
        $curr_distance = ($last_bitscore/$top_bitscore) * 100;
        undef(%hash);
        %hash = ();
        if($curr_distance > $r){ #register for vote.
            $hash{kingdom}{$kingdom}++;
            $hash{phylum}{$phylum}++;
            $hash{class}{$class}++;
            $hash{order}{$order}++;
            $hash{family}{$family}++;
            $hash{genus}{$genus}++;
            $hash{species}{$species}++;
        
        }else{
            print STDERR "Reject: $curr_distance lower than 10%\n" if($verbose);
        }
    }
}
close(IN);
#print STDERR Dumper(\%hash_lca);

# Then proceed to compute sum, fraction of Bsum and mbs for each taxon.
# and store in final hash
my %final_hash;
foreach my $contig_id (keys %hash_lca) {
    # compute scores for each taxon of each orf of each contig.
    my (@kingdom, @phylum, @class, @order, @family, @genus, @species);

    my %tmp_hash;
    foreach my $gene_id (keys %{ $hash_lca{$contig_id} }) {

        # get taxon string
        my $kingdom = $hash_lca{$contig_id}{$gene_id}{kingdom}{taxon}; 
        my $phylum  = $hash_lca{$contig_id}{$gene_id}{phylum}{taxon}; 
        my $class   = $hash_lca{$contig_id}{$gene_id}{class}{taxon}; 
        my $order   = $hash_lca{$contig_id}{$gene_id}{order}{taxon}; 
        my $family  = $hash_lca{$contig_id}{$gene_id}{family}{taxon}; 
        my $genus   = $hash_lca{$contig_id}{$gene_id}{genus}{taxon}; 
        my $species = $hash_lca{$contig_id}{$gene_id}{species}{taxon};

        #get bitscore
        my $kb = $hash_lca{$contig_id}{$gene_id}{kingdom}{bitscore};
        my $pb = $hash_lca{$contig_id}{$gene_id}{phylum}{bitscore};
        my $cb = $hash_lca{$contig_id}{$gene_id}{class}{bitscore};
        my $ob = $hash_lca{$contig_id}{$gene_id}{order}{bitscore};
        my $fb = $hash_lca{$contig_id}{$gene_id}{family}{bitscore};
        my $gb = $hash_lca{$contig_id}{$gene_id}{genus}{bitscore};
        my $sb = $hash_lca{$contig_id}{$gene_id}{species}{bitscore};
        
        $tmp_hash{kingdom}{$gene_id}{$kingdom} += $kb;
        $tmp_hash{phylum}{$gene_id}{$phylum}   += $pb;
        $tmp_hash{class}{$gene_id}{$class}     += $cb;
        $tmp_hash{order}{$gene_id}{$order}     += $ob;
        $tmp_hash{family}{$gene_id}{$family}   += $fb;
        $tmp_hash{genus}{$gene_id}{$genus}     += $gb;
        $tmp_hash{species}{$gene_id}{$species} += $sb;
        
    }
    my ($bsum_kingdom, $bsum_phylum, $bsum_class, $bsum_order, $bsum_family, $bsum_genus, $bsum_species);
    my ($bsumf_kingdom, $bsumf_phylum, $bsumf_class, $bsumf_order, $bsumf_family, $bsumf_genus, $bsumf_species);
    my ($mbs_kingdom, $mbs_phylum, $mbs_class, $mbs_order, $mbs_family, $mbs_genus, $mbs_species);
    #print STDERR "{----\n";
    #print STDERR Dumper(%tmp_hash); 
    #print STDERR "----}\n";

    my ($kingdom, $phylum, $class, $order, $family, $genus, $species);
   
    #########################
    ## KINGDOM
    ######################### 
    my %tmp_hash2; 
    foreach my $gene_id (keys %{ $tmp_hash{kingdom} }) {
        foreach my $taxon (keys %{ $tmp_hash{kingdom}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{kingdom}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    my $max = 0; my $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $bsum_kingdom = $max;
    $kingdom = $taxon;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon kingdom: $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);

    #########################
    ## PHYLUM
    ######################### 
    undef(%tmp_hash2);
    foreach my $gene_id (keys %{ $tmp_hash{phylum} }) {
        foreach my $taxon (keys %{ $tmp_hash{phylum}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{phylum}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    $max = 0; $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $bsum_phylum = $max;
    $phylum = $taxon;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon phylum : $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);

    #########################
    ## CLASS
    ######################### 
    undef(%tmp_hash2);
    foreach my $gene_id (keys %{ $tmp_hash{class} }) {
        foreach my $taxon (keys %{ $tmp_hash{class}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{class}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    $max = 0; $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $bsum_class = $max;
    $class = $taxon;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon class : $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);

    #########################
    ## ORDER
    ######################### 
    undef(%tmp_hash2);
    foreach my $gene_id (keys %{ $tmp_hash{order} }) {
        foreach my $taxon (keys %{ $tmp_hash{order}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{order}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    $max = 0; $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $bsum_order = $max;
    $order = $taxon;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon order : $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);

    #########################
    ## FAMILY
    ######################### 
    undef(%tmp_hash2);
    foreach my $gene_id (keys %{ $tmp_hash{family} }) {
        foreach my $taxon (keys %{ $tmp_hash{family}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{family}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    $max = 0; $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $bsum_family = $max;
    $family = $taxon;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon family : $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);

    #########################
    ## GENUS
    ######################### 
    undef(%tmp_hash2);
    foreach my $gene_id (keys %{ $tmp_hash{genus} }) {
        foreach my $taxon (keys %{ $tmp_hash{genus}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{genus}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    $max = 0; $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $bsum_genus = $max;
    $genus = $taxon;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon genus : $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);
    
    #########################
    ## SPECIES
    ######################### 
    undef(%tmp_hash2);
    foreach my $gene_id (keys %{ $tmp_hash{species} }) {
        foreach my $taxon (keys %{ $tmp_hash{species}{$gene_id} }) {
            #print STDERR "Taxon: $taxon\n";
            if($taxon ne "undef"){
                $tmp_hash2{$taxon} += $tmp_hash{species}{$gene_id}{$taxon}; #get bitscore;
            }
        }
    }
    $max = 0; $taxon = "";
    for my $key (keys %tmp_hash2){
        if($tmp_hash2{$key} > $max){
            $max = $tmp_hash2{$key};
            $taxon = $key;
        }
    }
    $species = $taxon;
    $bsum_species = $max;
    print STDERR "{====\n" if($verbose);
    print STDERR Dumper(%tmp_hash2) if($verbose);
    print STDERR "Taxon species : $taxon ; with bitscore of $max was voted\n" if($verbose);
    print STDERR "====}\n" if($verbose);

    # then process all orfs.
    if($bsum_kingdom != 0){
        
        $bsumf_kingdom = sprintf("%.2f", $bsum_kingdom/$bsum_kingdom);
        $bsumf_phylum  = sprintf("%.2f", $bsum_phylum/$bsum_kingdom);
        $bsumf_class   = sprintf("%.2f", $bsum_class/$bsum_kingdom);
        $bsumf_order   = sprintf("%.2f", $bsum_order/$bsum_kingdom);
        $bsumf_family  = sprintf("%.2f", $bsum_family/$bsum_kingdom);
        $bsumf_genus   = sprintf("%.2f", $bsum_genus/$bsum_kingdom);
        $bsumf_species = sprintf("%.2f", $bsum_species/$bsum_kingdom);
        
        my $mbs = $f * $bsum_kingdom;
        print STDOUT $contig_id."\t";
        print STATS $contig_id."\t";
       
        if($bsum_kingdom > $mbs){
            print STDOUT "\t".$kingdom;
            print STATS "\t".$bsumf_kingdom;
            $final_hash{$contig_id}{kingdom} = $kingdom;
            $final_hash{$contig_id}{kingdom_bsumf} = $bsumf_kingdom;
            
            if($bsum_phylum > $mbs){
                print STDOUT "\t".$phylum;
                print STATS  "\t".$bsumf_phylum;
                $final_hash{$contig_id}{phylum} = $phylum;
                $final_hash{$contig_id}{phylum_bsumf} = $bsumf_phylum;
                
                if($bsum_class > $mbs){
                    print STDOUT "\t".$class;
                    print STATS "\t".$bsumf_class;
                    $final_hash{$contig_id}{class} = $class;
                    $final_hash{$contig_id}{class_bsumf} = $bsumf_class;
                
                    if($bsum_order > $mbs){
                        print STDOUT "\t".$order;
                        print STATS "\t".$bsumf_order;
                        $final_hash{$contig_id}{order} = $order;
                        $final_hash{$contig_id}{order_bsumf} = $bsumf_order;
                
                        if($bsum_family > "\t".$mbs){
                            print STDOUT "\t".$family;
                            print STATS "\t".$bsumf_family;
                            $final_hash{$contig_id}{family} = $family;
                            $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
                    
                            if($bsum_genus > $mbs){
                                print STDOUT "\t".$genus;
                                print STATS "\t".$bsumf_genus;
                                $final_hash{$contig_id}{genus} = $genus;
                                $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
                        
                                if($bsum_species > $mbs){
                                    print STDOUT "\t".$species."\n";
                                    print STATS "\t".$bsumf_species."\n";
                                    $final_hash{$contig_id}{species} = $species;
                                    $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
                                }else{
                                    print STDOUT "\tundef\n";
                                    print STATS "\t$bsumf_species\n";
                                    $final_hash{$contig_id}{species} = "undef";
                                    $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
                                }
                            }else{
                                print STDOUT "\tundef\tundef\n";
                                print STATS "\t$bsumf_genus\t$bsumf_species\n";
                                $final_hash{$contig_id}{genus} = "undef";
                                $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
                                $final_hash{$contig_id}{species} = "undef";
                                $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
                            }
                        }else{
                            print STDOUT "\tundef\tundef\tundef\n";
                            print STATS "\t$bsumf_family\t$bsumf_genus\t$bsumf_species\n";
                            $final_hash{$contig_id}{family} = "undef";
                            $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
                            $final_hash{$contig_id}{genus} = "undef";
                            $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
                            $final_hash{$contig_id}{species} = "undef";
                            $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
                        }
                    }else{
                        print STDOUT "\tundef\tundef\tundef\tundef\n";
                        print STATS "\t$bsumf_order\t$bsumf_family\t$bsumf_genus\t$bsumf_species\n";
                        $final_hash{$contig_id}{order} = "undef";
                        $final_hash{$contig_id}{order_bsumf} = $bsumf_order;
                        $final_hash{$contig_id}{family} = "undef";
                        $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
                        $final_hash{$contig_id}{genus} = "undef";
                        $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
                        $final_hash{$contig_id}{species} = "undef";
                        $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
                    }
                }else{
                    print STDOUT "\tundef\tundef\tundef\tundef\tundef\n";
                    print STATS "\t$bsumf_class\t$bsumf_order\t$bsumf_family\t$bsumf_genus\t$bsumf_species\n";
                    $final_hash{$contig_id}{class} = "undef";
                    $final_hash{$contig_id}{class_bsumf} = $bsumf_class;
                    $final_hash{$contig_id}{order} = "undef";
                    $final_hash{$contig_id}{order_bsumf} = $bsumf_order;
                    $final_hash{$contig_id}{family} = "undef";
                    $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
                    $final_hash{$contig_id}{genus} = "undef";
                    $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
                    $final_hash{$contig_id}{species} = "undef";
                    $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
                }
            }else{
                print STDOUT "\tundef\tundef\tundef\tundef\tundef\tundef\n";
                print STATS "\t$bsumf_phylum\t$bsumf_class\t$bsumf_order\t$bsumf_family\t$bsumf_genus\t$bsumf_species\n";
                $final_hash{$contig_id}{phylum} = "undef";
                $final_hash{$contig_id}{phylum_bsumf} = $bsumf_phylum;
                $final_hash{$contig_id}{class} = "undef";
                $final_hash{$contig_id}{class_bsumf} = $bsumf_class;
                $final_hash{$contig_id}{order} = "undef";
                $final_hash{$contig_id}{order_bsumf} = $bsumf_order;
                $final_hash{$contig_id}{family} = "undef";
                $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
                $final_hash{$contig_id}{genus} = "undef";
                $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
                $final_hash{$contig_id}{species} = "undef";
                $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
            }
        }else{
            print STDOUT "\tundef\tundef\tundef\tundef\tundef\tundef\tundef\n"; 
            print STATS "\t$bsumf_kingdom\t$bsumf_phylum\t$bsumf_class\t$bsumf_order\t$bsumf_family\t$bsumf_genus\t$bsumf_species\n";
            $final_hash{$contig_id}{kingdom} = "undef";
            $final_hash{$contig_id}{kingdom_bsumf} = $bsumf_kingdom;
            $final_hash{$contig_id}{phylum} = "undef";
            $final_hash{$contig_id}{phylum_bsumf} = $bsumf_phylum;
            $final_hash{$contig_id}{class} = "undef";
            $final_hash{$contig_id}{class_bsumf} = $bsumf_class;
            $final_hash{$contig_id}{order} = "undef";
            $final_hash{$contig_id}{order_bsumf} = $bsumf_order;
            $final_hash{$contig_id}{family} = "undef";
            $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
            $final_hash{$contig_id}{genus} = "undef";
            $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
            $final_hash{$contig_id}{species} = "undef";
            $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
        }

    }else{
        print STDOUT "$contig_id\tundef\tundef\tundef\tundef\tundef\tundef\tundef\n"; 
        print STATS "$contig_id\t0\t0\t0\t0\t0\t0\t0\n";
        $final_hash{$contig_id}{kingdom} = "undef";
        $final_hash{$contig_id}{kingdom_bsumf} = $bsumf_kingdom;
        $final_hash{$contig_id}{phylum} = "undef";
        $final_hash{$contig_id}{phylum_bsumf} = $bsumf_phylum;
        $final_hash{$contig_id}{class} = "undef";
        $final_hash{$contig_id}{class_bsumf} = $bsumf_class;
        $final_hash{$contig_id}{order} = "undef";
        $final_hash{$contig_id}{order_bsumf} = $bsumf_order;
        $final_hash{$contig_id}{family} = "undef";
        $final_hash{$contig_id}{family_bsumf} = $bsumf_family;
        $final_hash{$contig_id}{genus} = "undef";
        $final_hash{$contig_id}{genus_bsumf} = $bsumf_genus;
        $final_hash{$contig_id}{species} = "undef";
        $final_hash{$contig_id}{species_bsumf} = $bsumf_species;
    }
}

# If a contig blastn of bins contigs seqs vs ncbi genomes is provided:
my %hash_blastn1;
my %hash_blastn2;
if($infile_taxonomy_blastn){
    open(IN, "<".$infile_taxonomy_blastn) or die "Can't open $infile_taxonomy_blastn\n";
    print STDERR "Complement CAT taxonomy with blastn results.\n";

    my $i = 0;
    while(<IN>){
        chomp;
        next if($. == 1);
        my @row = split(/\t/, $_);
        my $contig_id = $row[0];
        my $ratio_alignlength_on_qlen = $row[10] / $row[12]; 
        
        if($row[9] < $evalue && $ratio_alignlength_on_qlen > $ratio){
            $hash_blastn1{$contig_id}{$i}{al_length}  = $row[10];  
            $hash_blastn1{$contig_id}{$i}{evalue}     = $row[9];  
            $hash_blastn1{$contig_id}{$i}{qlen}       = $row[12];  
            $hash_blastn1{$contig_id}{$i}{kingdom}    = $row[2];
            $hash_blastn1{$contig_id}{$i}{phylum}     = $row[3];
            $hash_blastn1{$contig_id}{$i}{class}      = $row[4];
            $hash_blastn1{$contig_id}{$i}{order}      = $row[5];
            $hash_blastn1{$contig_id}{$i}{family}     = $row[6];
            $hash_blastn1{$contig_id}{$i}{genus}      = $row[7];
            $hash_blastn1{$contig_id}{$i}{species}    = $row[8];
            $hash_blastn1{$contig_id}{$i}{ratio}      = $ratio_alignlength_on_qlen;
            $i++;
        } 
    }
    close(IN);
    #print STDERR Dumper(\%hash_blastn1);
    
    #Then loop through these hits and create a consensus taxonomy
    foreach my $contig_id (keys %hash_blastn1){
        my %tmp_hash;

        foreach my $i (keys %{ $hash_blastn1{$contig_id} }) {
            my $kingdom = $hash_blastn1{$contig_id}{$i}{kingdom};
            my $phylum  = $hash_blastn1{$contig_id}{$i}{phylum};
            my $class   = $hash_blastn1{$contig_id}{$i}{class};
            my $order   = $hash_blastn1{$contig_id}{$i}{order};
            my $family  = $hash_blastn1{$contig_id}{$i}{family};
            my $genus   = $hash_blastn1{$contig_id}{$i}{genus};
            my $species = $hash_blastn1{$contig_id}{$i}{species};

            $tmp_hash{kingdom}{$kingdom}++;
            $tmp_hash{phylum}{$phylum}++;
            $tmp_hash{class}{$class}++;
            $tmp_hash{order}{$order}++;
            $tmp_hash{family}{$family}++;
            $tmp_hash{genus}{$genus}++;
            $tmp_hash{species}{$species}++;
        }

        #Then take generate consensus lineage starting from kingdom->species by majority rule voting process. Maybe consider that if no winners, stop constructing the lineage.
        my $blastn_kingdom = "undef";
        my $blastn_phylum  = "undef";
        my $blastn_class   = "undef";
        my $blastn_order   = "undef";
        my $blastn_family  = "undef";
        my $blastn_genus   = "undef";
        my $blastn_species = "undef";
        
        #my $blastn_kingdom_count = 0; my $blastn_phylum_count = 0; my $blastn_class_count = 0; my $blastn_order_count = 0; my $blastn_family_count = 0; my $blastn_genus_count = 0; my $blastn_species_count = 0;
        
        #print STDERR Dumper(\%tmp_hash);
        
        $blastn_kingdom = max_by { $tmp_hash{kingdom}{$_} } keys %{ $tmp_hash{kingdom} };
        #print STDERR $blastn_kingdom."\n";
        
        $blastn_phylum = max_by { $tmp_hash{phylum}{$_} } keys %{ $tmp_hash{phylum} };
        #print STDERR $blastn_phylum."\n";
        
        $blastn_class = max_by { $tmp_hash{class}{$_} } keys %{ $tmp_hash{class} };
        #print STDERR $blastn_class."\n";
        
        $blastn_order = max_by { $tmp_hash{order}{$_} } keys %{ $tmp_hash{order} };
        #print STDERR $blastn_order."\n";
        
        $blastn_family = max_by { $tmp_hash{family}{$_} } keys %{ $tmp_hash{family} };
        #print STDERR $blastn_family."\n";
        
        $blastn_genus = max_by { $tmp_hash{genus}{$_} } keys %{ $tmp_hash{genus} };
        #print STDERR $blastn_genus."\n";
        
        $blastn_species = max_by { $tmp_hash{species}{$_} } keys %{ $tmp_hash{species} };
        #print STDERR $blastn_species."\n";
       
        #$hash_blastn2{$contig_id}{al_length} = $row[10];  
        #$hash_blastn2{$contig_id}{evalue}    = $row[9];  
        #$hash_blastn2{$contig_id}{qlen}      = $row[12];  
        #$hash_blastn2{$contig_id}{ratio}     = $ratio_alignlength_on_qlen;
        $hash_blastn2{$contig_id}{kingdom}   = $blastn_kingdom;
        $hash_blastn2{$contig_id}{phylum}    = $blastn_phylum;
        $hash_blastn2{$contig_id}{class}     = $blastn_class;
        $hash_blastn2{$contig_id}{order}     = $blastn_order;
        $hash_blastn2{$contig_id}{family}    = $blastn_family;
        $hash_blastn2{$contig_id}{genus}     = $blastn_genus;
        $hash_blastn2{$contig_id}{species}   = $blastn_species;
       
        my $blastn_depth = 0;
        if($hash_blastn2{$contig_id}{kingdom} ne "undef"){ $blastn_depth = 1; }
        if($hash_blastn2{$contig_id}{phylum} ne "undef"){  $blastn_depth = 2; }
        if($hash_blastn2{$contig_id}{class} ne "undef"){   $blastn_depth = 3; }
        if($hash_blastn2{$contig_id}{order} ne "undef"){   $blastn_depth = 4; }
        if($hash_blastn2{$contig_id}{family} ne "undef"){  $blastn_depth = 5; }
        if($hash_blastn2{$contig_id}{genus} ne "undef"){   $blastn_depth = 6; }
        if($hash_blastn2{$contig_id}{species} ne "undef"){ $blastn_depth = 7; }
        $hash_blastn2{$contig_id}{depth} = $blastn_depth;

     }

#    while(<IN>){
#        chomp;
#        next if($. == 1);
#        my @row = split(/\t/, $_);
#        my $contig_id = $row[0];
#        my $ratio_alignlength_on_qlen = $row[10] / $row[12]; 
#        
#        if($row[9] < $evalue && $ratio_alignlength_on_qlen > $ratio){
#            $hash_blastn2{$contig_id}{al_length} = $row[10];  
#            $hash_blastn2{$contig_id}{evalue}    = $row[9];  
#            $hash_blastn2{$contig_id}{qlen}      = $row[12];  
#            $hash_blastn2{$contig_id}{ratio}     = $ratio_alignlength_on_qlen;
#            $hash_blastn2{$contig_id}{kingdom}   = $row[2];
#            $hash_blastn2{$contig_id}{phylum}    = $row[3];
#            $hash_blastn2{$contig_id}{class}     = $row[4];
#            $hash_blastn2{$contig_id}{order}     = $row[5];
#            $hash_blastn2{$contig_id}{family}    = $row[6];
#            $hash_blastn2{$contig_id}{genus}     = $row[7];
#            $hash_blastn2{$contig_id}{species}   = $row[8];
#              
#            my $blastn_depth = 0;
#            if($hash_blastn2{$contig_id}{kingdom} ne "undef"){ $blastn_depth = 1; }
#            if($hash_blastn2{$contig_id}{phylum} ne "undef"){  $blastn_depth = 2; }
#            if($hash_blastn2{$contig_id}{class} ne "undef"){   $blastn_depth = 3; }
#            if($hash_blastn2{$contig_id}{order} ne "undef"){   $blastn_depth = 4; }
#            if($hash_blastn2{$contig_id}{family} ne "undef"){  $blastn_depth = 5; }
#            if($hash_blastn2{$contig_id}{genus} ne "undef"){   $blastn_depth = 6; }
#            if($hash_blastn2{$contig_id}{species} ne "undef"){ $blastn_depth = 7; }
#            $hash_blastn2{$contig_id}{depth} = $blastn_depth;
#        }
#    }
#    close(IN);

    my $potentially_better_classified_contigs = 0;
    for my $contig_id (keys %final_hash){
        my $cat_depth = 0;
        if(exists $hash_blastn2{$contig_id}){
           if($final_hash{$contig_id}{kingdom} ne "undef"){ $cat_depth = 1; }
           if($final_hash{$contig_id}{phylum} ne "undef"){  $cat_depth = 2; }
           if($final_hash{$contig_id}{class} ne "undef"){   $cat_depth = 3; }
           if($final_hash{$contig_id}{order} ne "undef"){   $cat_depth = 4; }
           if($final_hash{$contig_id}{family} ne "undef"){  $cat_depth = 5; }
           if($final_hash{$contig_id}{genus} ne "undef"){   $cat_depth = 6; }
           if($final_hash{$contig_id}{species} ne "undef"){ $cat_depth = 7; }
           
           if($cat_depth == $hash_blastn2{$contig_id}{depth}){ #if equal go with CAT results.
               print OUT $contig_id."\t".$final_hash{$contig_id}{kingdom}."\t".$final_hash{$contig_id}{phylum}."\t".$final_hash{$contig_id}{class}."\t".$final_hash{$contig_id}{order}."\t".$final_hash{$contig_id}{family}."\t".$final_hash{$contig_id}{genus}."\t".$final_hash{$contig_id}{species}."\n";
           }elsif($cat_depth > $hash_blastn2{$contig_id}{depth}){
               print OUT $contig_id."\t".$final_hash{$contig_id}{kingdom}."\t".$final_hash{$contig_id}{phylum}."\t".$final_hash{$contig_id}{class}."\t".$final_hash{$contig_id}{order}."\t".$final_hash{$contig_id}{family}."\t".$final_hash{$contig_id}{genus}."\t".$final_hash{$contig_id}{species}."\n";

           }elsif($cat_depth < $hash_blastn2{$contig_id}{depth}){
               print OUT $contig_id."\t".$hash_blastn2{$contig_id}{kingdom}."\t".$hash_blastn2{$contig_id}{phylum}."\t".$hash_blastn2{$contig_id}{class}."\t".$hash_blastn2{$contig_id}{order}."\t".$hash_blastn2{$contig_id}{family}."\t".$hash_blastn2{$contig_id}{genus}."\t".$hash_blastn2{$contig_id}{species}."\n";
               $potentially_better_classified_contigs++;
               #print STDERR "{---\n";
               #print STDERR $contig_id."\t".$hash_blastn2{$contig_id}{kingdom}."\t".$hash_blastn2{$contig_id}{phylum}."\t".$hash_blastn2{$contig_id}{class}."\t".$hash_blastn2{$contig_id}{order}."\t".$hash_blastn2{$contig_id}{family}."\t".$hash_blastn2{$contig_id}{genus}."\t".$hash_blastn2{$contig_id}{species}."\n";
               #print STDERR "ratio: $hash_blastn2{$contig_id}{ratio}\n";
               #print STDERR "evalue: $hash_blastn2{$contig_id}{evalue}\n";
               #print STDERR "al_length: $hash_blastn2{$contig_id}{al_length}\n";
               #print STDERR "qlen: $hash_blastn2{$contig_id}{qlen}\n";
               #print STDERR "---}\n";
           }else{
               print STDERR "contig_id: $contig_id         cat_depth:$cat_depth     blastn_depth:$hash_blastn2{$contig_id}{depth}\n" if($verbose);
           }
        }else{
            print OUT $contig_id."\t".$final_hash{$contig_id}{kingdom}."\t".$final_hash{$contig_id}{phylum}."\t".$final_hash{$contig_id}{class}."\t".$final_hash{$contig_id}{order}."\t".$final_hash{$contig_id}{family}."\t".$final_hash{$contig_id}{genus}."\t".$final_hash{$contig_id}{species}."\n";
        }
    }
    print STDERR "Potentially improved contigs classification: $potentially_better_classified_contigs\n";
}
close(OUT);
