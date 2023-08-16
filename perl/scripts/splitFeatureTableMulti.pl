#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
splitFeatureTable.pl

PURPOSE:
Split Feature table based on the input --select <string> value. All lineage matching to that value will be kept in the --matched <string> outfile.

INPUT:
--infile <string>           : Feature table in Qiime format (qiime.org)
--select <string>           : Either 'bacteriaArchaea' or 'fungi' or 'bacteria' or 'eukaryotes' or 'all'
--keepChloro <string>       : 'yes' or 'no'. Default ='yes'      
--keepMito <string>         : 'yes' or 'no'. Default = 'yes'

OUTPUT:
--selection.tsv <string>    : Feature table having only <selection>.
--unmatched <string>        : Feature table having non-Bacteria and Non-Archeal organisms.

NOTES:


BUGS/LIMITATIONS:

 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $unmatched, $select, $keepChloro, $keepMito, $outdir);
my $verbose = 0;

GetOptions(
  'infile=s'     => \$infile,
  'outdir=s'     => \$outdir,
  'unmatched=s'  => \$unmatched,
  'select=s'     => \$select,
  'keepChloro=s' => \$keepChloro,
  'keepMito=s'   => \$keepMito,
  'verbose'      => \$verbose,
  'help'         => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile <string> required\n" unless($infile);
die "--unmatched <string> required\n" unless($unmatched);
#die "--select must be either 'bactArch' or 'fungi' or 'bacteria' or 'archaea' or 'all'" if($select ne "bacteriaArchaea" && $select ne "fungi" && $select ne "bacteria" && $select ne "archaea" && $select ne "all" && $select ne "eukaryotes");
$keepChloro = "yes" unless($keepChloro);
$keepChloro = lc($keepChloro);
$keepMito = "yes" unless($keepMito);
$keepMito = lc($keepMito);

## MAIN

# TODO Put that in a separate script.
# This subroutine takes a feature_table in input and gives a otu table minus every lineage that does not start with k__Bacteria or k__Archaea in output.
# @input : feature_table.tab
# @output : feature_table.tab (only containing bacteria/archaea organisms).
open(IN, "<".$infile) or die "Can't open input file ".$infile."\n";
#open(OUT, ">".$matched) or die "Can't open input file ".$matched."\n";
open(OUT_F, ">".$unmatched) or die "Can't open input file ".$unmatched."\n";

my @regexes = split(/,/, $select);
#$regex = "k__bacteria|k__archaea" if($select eq "bacteriaArchaea");
#$regex = "k__fungi|o__fungi" if($select eq "fungi");
#$regex = "k__bacteria" if($select eq "bacteria");
#$regex = "k__archaea" if($select eq "archaea");
#$regex = "^" if($select eq "everything");
#$regex = "k__eukaryota" if($select eq "eukaryotes");

print STDERR "Keepchloro? ".$keepChloro."\n";
print STDERR "KeepMito? ".$keepMito."\n";

my $name = $infile;
$name =~ s{^.*/}{};     # remove the leading path  
$name =~ s{\.[^.]+$}{}; # remove the extension
my %handles;
for my $regex (@regexes){
    open my $handle, ">", $outdir."/".$name."_".$regex.".tsv" or die $!;
    $handles{$outdir."/".$name."_".$regex.".tsv"} = $handle;
}

while(<IN>){
   chomp;
   # header
   if($_ =~ m/#/){
       for my $regex (@regexes){
           #if($_ =~ m/($regex)/i){
           print { $handles{$outdir."/".$name."_".$regex.".tsv"} } $_."\n";
           #}
       }
       print OUT_F $_."\n";
       next;
    }

   # selection
   my $redFlag = "true";
   for my $regex (@regexes){
       if($_ =~ m/($regex)/i){
           print { $handles{$outdir."/".$name."_".$regex.".tsv"} } $_."\n";
           $redFlag = "false";
       }
    }
    if($redFlag eq "true"){
       print OUT_F $_."\n";
    }    
}
close(IN);
#close(OUT);
close(OUT_F);
close $_ foreach values %handles;

exit;
