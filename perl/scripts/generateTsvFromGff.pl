#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateTsvFromGff.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT <string>

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
   'infile=s'  => \$infile,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    my @row = split(/\t/, $_);
    my $last_row = pop(@row);
    my @last_row = split(/;/, $last_row);


    my $gene_id;
    my $dbxref;      
    my $name;       
    my $chr;         
    my $gbkey;        
    my $genome;       
    my $mol_type;    
    my $product;      
    my $transcript_id;
    my $gene;
   
    foreach my $el (@last_row){
        print STDERR $el."\n";
        $gene_id = $1       if($el =~ /gene_id=(.*)/);
        $dbxref = $1        if($el =~ /Dbxref=(.*)/);
        $name = $1          if($el =~ /Name=(.*)/);
        $chr = $1           if($el =~ /chromosome=(.*)/);
        $gbkey = $1         if($el =~ /gbkey=(.*)/);
        $genome = $1        if($el =~ /genome=(.*)/);
        $mol_type = $1      if($el =~ /mol_type=(.*)/);
        $product = $1       if($el =~ /product=(.*)/);
        $transcript_id = $1 if($el =~ /transcript_id=(.*)/);  
        $gene = $1          if($el =~ /gene=(.*)/);  
    }
    #print STDERR $gene_id."\n";
    $hash{$gene_id}{dbxref}          = $dbxref;      
    $hash{$gene_id}{name}            = $name;       
    $hash{$gene_id}{chr}             = $chr;         
    $hash{$gene_id}{gbkey}           = $gbkey;        
    $hash{$gene_id}{genome}          = $genome;       
    $hash{$gene_id}{mol_type}        = $mol_type;    
    $hash{$gene_id}{product}         = $product;     
    $hash{$gene_id}{transcript_id}   = $transcript_id;
    $hash{$gene_id}{gene}            = $gene;

    #print STDERR $gene_id."\n";
}
print STDERR Dumper(\%hash);

for my $k (sort keys %hash) {
    #print STDERR "$k\n";
   print STDOUT "$k";
   if(defined $hash{$k}{dbxref})        {print STDOUT "\t".$hash{$k}{dbxref};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{name})          {print STDOUT "\t".$hash{$k}{name};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{chr})           {print STDOUT "\t".$hash{$k}{chr};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{gbkey})         {print STDOUT "\t".$hash{$k}{gbkey};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{genome})        {print STDOUT "\t".$hash{$k}{genome};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{mol_type})      {print STDOUT "\t".$hash{$k}{mol_type};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{product})       {print STDOUT "\t".$hash{$k}{product};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{transcript_id}) {print STDOUT "\t".$hash{$k}{transcript_id};} else {print STDOUT "\tNA";}
   if(defined $hash{$k}{gene})          {print STDOUT "\t".$hash{$k}{gene};} else {print STDOUT "\tNA";}
   print STDOUT "\n";
}

#c("GeneID", "dbxref", "name", "chr", "gbkey", "genome", "mol_type", "product", "transcript_id")
