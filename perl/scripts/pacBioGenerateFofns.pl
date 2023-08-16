#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
pacBioGenerateFofns.pl

PURPOSE:

INPUT:
--sample_name <string> : Sequence file
--readset <string>     : readset.tsv file.

OUTPUT:
outdir <string>        : outdir where will be written the fofn for
                         that specific <sample_name>

NOTES:
readset.tsv file should look like this with each field separated by a tab:

Sample    Smartcell     EstimatedGenomesize   BAS                BAX
LB501T    1             6500000               <leave empty>      ./raw_reads/m170517_223644_42158_c101204432550000001823285710241730_s1_p0.3.bax.h5
LB501T    2             6500000               <leave empty>      ./raw_reads/m170517_223644_42158_c101204432550000001823285710241730_s1_p0.2.bax.h5
LB501T    3             6500000               <leave empty>      ./raw_reads/m170517_223644_42158_c101204432550000001823285710241730_s1_p0.1.bax.h5


BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $sample_name, $readset, $outdir);
my $verbose = 0;

GetOptions(
   'sample_name=s'   => \$sample_name,
   'readset=s'       => \$readset,
   'outdir=s'        => \$outdir,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN

#system("mkdir -p ".$outdir) or die "Can't mkdir -p $outdir\n";
mkdir $outdir unless(-d $outdir);

my %hash;
my $i = 0;
open(IN, "<".$readset) or die "Can't open $readset\n";
while(<IN>){
    chomp;
    next if($_ =~ m/^Sample\t/);
    my @row = split(/\t/, $_);

    my $curr_sample_name = $row[0];
    if($curr_sample_name eq $sample_name){
        my $file_path = $row[4];

        #get basename of movie file name.
        my $base = $file_path;
        $base =~ s{^.*/}{};     # remove the leading path  
        $base =~ s{\.[^.]+$}{}; # remove the extension
        $base =~ s/_s1_p0\..*//;
        
        print STDERR $base."\n";

        $hash{$curr_sample_name}{$base}{$i} = $file_path;
        $i++;
    } 
}

#print STDERR Dumper(\%hash);
my $j = 1;
foreach my $key (keys %hash){
    foreach my $key2 (keys %{ $hash{$key} }) {
        open(OUT, ">".$outdir."/".$key."_".$j."_fofn.txt") or die "Can't open ".$outdir."/".$key."_fofn.txt\n";

        foreach my $key3 (keys %{ $hash{$key}{$key2} }) {
            print OUT $hash{$key}{$key2}{$key3}."\n";
        }
        $j++;
    }
    close(OUT);
}


#my %hash;
#my $i = 0;
#open(IN, "<".$readset) or die "Can't open $readset\n";
#while(<IN>){
#    next if($_ =~ m/^Sample\t/);
#    my @row = split(/\t/, $_);
#    my $curr_sample_name = $row[0];
#    my $file_path = $row[5];
#    $hash{$curr_sample_name}{$i} = $file_path;
#
#    $i++; 
#}
#
#foreach my $key (keys %hash){
#    open(OUT, ">".$outdir."/".$key."_fofn.txt") or die "Can't open ".$outdir."/".$key."_fofn.txt\n";
#    foreach my $key2 (keys %{ $hash{$key}{$key2} }) {
#        print OUT $hash{$key}{$key2}."\n";
#    }
#}

