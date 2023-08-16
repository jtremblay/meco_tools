#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
addTaxToObsClosed.pl

PURPOSE:
This script takes a obs abundance table and it's corresponding 
taxonomy (generated with blastn, besthit.).
It appends the Greengenes lineages to each corresponding obs clusters 
and format the output as a Feature table compatible with QIIME.

INPUT:
--seqobs <tab_file>         :  obs output tablablature file (cluster abundance/indexes)
--blast_output <string>     :  output from blast (best hit only).
--cutoff <int>              :  Selection cutoff, value between 0.0 and 1. Default is 0.8.
--evalue <float>            :  evalue cutoff. evalue > than < --evalue will be discarded
--al_length <int>           :  alignment length cutoff. alignment length < --al_length will be discarded.
--perc_id <float>           :  percentage id cutoff. hit with perc_id < --perc_id will be discarded.

OUTPUT:
--outfile <string>          :  Feature table compatible with the Qiime tools.
--outfile_failed <string>   :  Feature table compatible with the Qiime tools. 
                               Contains values that didn't meet the cutoff requirement.
--link <string>             :  Tsv file that links the original open ref Feature id (1st column) with 
                               the closed ref Feature id (2nd column).

NOTES:

BUGS/LIMITATIONS:
  
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $seqobs, $blast_output, $outfile, $outfile_failed, $co_evalue, $co_al_length, $co_perc_id, $link);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'seqobs=s'           => \$seqobs,
  'blast_output=s'     => \$blast_output,
  'evalue=f'           => \$co_evalue,
  'al_length=i'        => \$co_al_length,
  'perc_id=f'          => \$co_perc_id,
  'outfile=s'          => \$outfile,  
  'outfile_failed=s'   => \$outfile_failed,  
  'link=s'             => \$link,  
  'verbose'            => \$verbose,
  'help'               => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--seqobs seqobs table (qiime format) file required.\n") unless $seqobs;
die("--blast_output blast_output table required.\n") unless $blast_output;
die("--outfile outfile required.\n") unless $outfile;
#$cutoff = 0.80 unless $cutoff;
#$tax_level = 0 unless $tax_level;
#die("--cutoff cutoff value must be between 0 and 1.0 inclusively.\n") unless($cutoff >= 0.00 or $cutoff <= 1.00);
#$cutoff = 0.00 if($cutoff == 0 or $cutoff == 0.0);
#$cutoff = int($cutoff);

## MAIN
open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";
open(OUT_FAILED, ">".$outfile_failed) or die "Can't open file ".$outfile_failed."\n";
open(LINK, ">".$link) or die "Can't open file ".$link."\n";
open(SEQOBS, $seqobs) or die "Can't open file ".$seqobs."\n";
open(BLAST, $blast_output) or die "Can't open file ".$blast_output."\n";

print OUT "#Full Feature Counts - closed ref\n";
print OUT_FAILED "#Full Feature Counts - closed ref\n";
#print OUT_SPECIFIC "#Full Feature Counts - closed ref\n";

my %hash = ();
my %hash_prob = ();
my %hash_prob_best = ();

#Put relevant rdp entries in a hash table
while(<BLAST>){
    chomp($_);
    my @row = split(/\t/, $_);
    #print $_."\n" if($verbose);
    my $query          = $row[0];
    my $subject        = $row[1];
    my $perc_id        = $row[2];
    my $al_length      = $row[3];
    my $evalue         = $row[10];
    my $id_and_lineage = $row[12];
    #foreach my $el (@row){
    #    print STDERR $el."\n";
    #}
    
    #4484200  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__BacteroidesSP;
    my $last_char = substr($id_and_lineage, -1);
    if($last_char eq ";"){ chop($id_and_lineage); }
    my @id_and_lineage = split(/;/, $id_and_lineage);
    my $id_and_k = shift(@id_and_lineage);
    my $k;
    my $id;
    #print STDERR "idandk1:"."\t".$id_and_k."\n";
    if($id_and_k =~ m/(^\d+)[ ]{1,}(\S+)/){
        $id = $1;
        $k = $2;
    }
    #print STDERR "idandk2:"."\t".$id_and_k."\n";
    #print STDERR "id"."\t".$id."\n";
    #print STDERR "k"."\t".$k."\n";

    my $p = $id_and_lineage[0];
    my $c = $id_and_lineage[1];
    my $o = $id_and_lineage[2];
    my $f = $id_and_lineage[3];
    my $g = $id_and_lineage[4];
    my $s = $id_and_lineage[5];

    $id =~ s/@//;
    $id =~ s/>//;
    if(substr($id, -1) eq "#"){
        $id = substr($id, 0, -1);
    }
  
    #$id = $1 if($id =~ m/(\d+)/);
 
    $hash{$query}{query} = $query;
    $hash{$query}{subject} = $subject;
    $hash{$query}{lineage} = $k.";".$p.";".$c.";".$o.";".$f.";".$g.";".$s;
    $hash{$query}{id} = $id;  
    $hash{$query}{evalue} = $evalue;  
    $hash{$query}{perc_id} = $perc_id;  
    $hash{$query}{al_length} = $al_length;  
  

}
close(BLAST);
#print STDERR (Dumper(\%hash));

#Loop through seqobs table and add rdp taxonomy accordingly
my $i=0; #for the header.
print LINK "FEATURE_ID_open\tFEATURE_ID_closed\n";
while(<SEQOBS>){
    chomp($_);
    if($i == 0){ #Print Qiime-style OTU table header.
        $_ =~ s/CLUSTER/FEATURE_ID/;
        print OUT $_."\ttaxonomy\n"; 
        print OUT_FAILED $_."\ttaxonomy\n"; 
        #print OUT_SPECIFIC $_."\ttaxonomy\n"; 
    }else{
        my @row = split(/\t/, $_);
        my $id = shift(@row);
        $_ = join("\t", @row); #reconstruct $_ string without original open ref otu id. 
        $id =~ s/@//;
        $id =~ s/>//;
        if(substr($id, -1) eq "#"){
            $id = substr($id, 0, -1);
        }
        #print STDERR $id."\n";
        if(exists ($hash{$id}) ){
            my $passed_cutoffs = 0;
            if($hash{$id}{evalue} <= $co_evalue && $hash{$id}{perc_id} >= $co_perc_id && $hash{$id}{al_length} >= $co_al_length){
                $passed_cutoffs = 1;
            }
            #if($hash{$id}{evalue} <= $co_evalue){
            #    $passed_cutoffs = 1;
            #}
            #if($hash{$id}{perc_id} >= $co_perc_id){
            #    $passed_cutoffs = 1;
            #}
            #if($hash{$id}{al_length} >= $co_al_length){
            #    $passed_cutoffs = 1;
            #}

            if($passed_cutoffs == 1){
                print OUT $hash{$id}{subject}."\t".$_."\t".$hash{$id}{lineage}."\n";
                print LINK $id."\t".$hash{$id}{subject}."\n";
                #if($hash{$id}{lineage} =~ m/k__Bacteria|k__Archaea/){
                    #my $curr_row = $_;
                    #my @curr_row = split(/\t/, $_);
                    #shift(@curr_row);
                    #$curr_row = join("\t", @curr_row);
                    #print OUT_SPECIFIC $hash{$id}{subject}."\t".$curr_row."\t".$hash{$id}{lineage}."\n" 
                #}
            }else{
                my $curr_row = $_;
                my @curr_row = split(/\t/, $_);
                shift(@curr_row);
                $curr_row = join("\t", @curr_row);
                print OUT_FAILED $hash{$id}{subject}."\t".$curr_row."\t".$hash{$id}{lineage}."\n" 
            }
            
        }else{
            print STDERR "Undefined taxonomy for : ".$id."==> no corresponding best blast hit on selected db?\n";
        }
    }
    $i=1;
}
close(SEQOBS);
close(OUT);
close(OUT_FAILED);
close(LINK);
exit;

