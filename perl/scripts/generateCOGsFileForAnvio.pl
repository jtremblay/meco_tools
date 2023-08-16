#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
generateCOGsFileForAnvio.pl

PURPOSE:

INPUT:
--infile <string>   : COGs RPSBLAST annotations file (shotgunMG pipeline).

OUTPUT:
<STDOUT>            : function file compatible with anvio file import.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

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
print STDOUT join("\t", ("gene_callers_id", "source", "accession", "function", "e_value")); print STDOUT "\n";
my %hash;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    my @row = split(/\t/, $_);
    
    my $gene_id = $row[0];
    $gene_id =~ s/gene_id_//;
    my $evalue = $row[10];
    
    my $annotation = $row[12];

    my $cog;
    my $cog_category;
    my $desc1;
    my $name;

    if($annotation =~ m/\[/){
        my @annotation2 = split(/\[/, $annotation);
        $cog_category = pop(@annotation2);
        my $desc2 = shift(@annotation2);
        #print STDERR $desc2."\n";
        my @desc2 = split(/, /, $desc2);
        $cog = shift(@desc2);
        $name = shift(@desc2);
        $desc1 = join(", ", @desc2);
    }else{
        #print STDERR $desc2."\n";
        my @desc2 = split(/, /, $annotation);
        $cog = shift(@desc2);
        $name = shift(@desc2);
        $desc1 = join(", ", @desc2);
        $cog_category = "unknown"
    }
    
    
    print STDOUT $gene_id."\tCOG_category\t".$cog."\t".$cog_category."\t".$evalue."\n";
    print STDOUT $gene_id."\tCOG\t".$cog."\t".$desc1."\t".$evalue."\n";
    print STDOUT $gene_id."\tCOG_name\t".$cog."\t".$name."\t".$evalue."\n";


}
close(IN);

