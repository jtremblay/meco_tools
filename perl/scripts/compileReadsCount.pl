#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
compileReadsCount.pl

PURPOSE:

INPUT:
--infiles <string>  : list of rRNATagger count report files separated by commas.
--prefixes <string> : list of identifier prefixes corresponding to identify the infiles <string>

OUTPUT:
STDOUT <string>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $prefixes);
my $verbose = 0;

GetOptions(
   'infiles=s'  => \$infiles,
   'prefixes=s' => \$prefixes,
   'verbose'    => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my @infiles = split(/,/, $infiles);
my @prefixes = split(/,/, $prefixes);

my %hash;
foreach my $infile (@infiles){
    my $prefix = shift(@prefixes);

    open(IN, "<".$infile) or die "Can't open $infile\n";
    while(<IN>){
        chomp;
        my @row = split(/\t/, $_);

        if($row[0] eq "total_reads"){
            $hash{$prefix}{total_reads} = $row[1];
        
        }elsif($row[0] eq "contaminants_reads"){
            $hash{$prefix}{contaminants_reads} = $row[1];
        
        }elsif($row[0] eq "phix_reads"){
            $hash{$prefix}{phix_reads} = $row[1];
        
        }elsif($row[0] eq "non_contam_non_phix_reads"){
            $hash{$prefix}{non_contam_non_phix_reads} = $row[1];
        
        }elsif($row[0] eq "non_contam_non_phix_reads_1"){
            $hash{$prefix}{non_contam_non_phix_reads_1} = $row[1];
        
        }elsif($row[0] eq "reads_1_QC_passed"){
            $hash{$prefix}{reads_1_QC_passed} = $row[1];
        
        #}elsif($_ =~ m/Cluster counts/){
        
        }elsif($row[0] eq "non_contam_non_phix_reads_2"){
            $hash{$prefix}{non_contam_non_phix_reads_2} = $row[1];
        
        }elsif($row[0] eq "assembled_reads"){
            $hash{$prefix}{assembled_reads} = $row[1];
        
        }elsif($row[0] eq "assembled_reads_QC_passed"){
            $hash{$prefix}{assembled_reads_QC_passed} = $row[1];
        
        ######################################### 
        }elsif($row[0] eq "Clustered:"){
            #385,900 sequences give a total of 5,180 clusters
            if($row[1] =~ m/(\S+) sequences.* of (\S+) clusters/){
                $hash{$prefix}{clustered_sequences} = $1;
                $hash{$prefix}{number_of_clusters} = $2;
            }
        
        }elsif($row[0] eq "Classifed clusters:"){
            if($row[1] =~ m/(\S+) sequences.*/){
                $hash{$prefix}{clustered_sequences_classified} = $1;
            }
        }
    }
    close(IN);
}
print STDERR Dumper(\%hash);

print STDOUT "project_name\ttotal_reads\tcontaminants_reads\tphix_reads\tnon_contam_non_phix_reads\tnon_contam_non_phix_reads_1\tnon_contam_non_phix_reads_2\treads_1_QC_passed\tassembled_reads\tassembled_reads_QC_passed\tclustered_sequences\tnumber_of_clusters\tclustered_sequences_classified\n";
foreach my $key (keys %hash){
    print STDOUT $key."\t".$hash{$key}{total_reads}."\t";
    print STDOUT $hash{$key}{contaminants_reads}."\t";
    print STDOUT $hash{$key}{phix_reads}."\t";
    print STDOUT $hash{$key}{non_contam_non_phix_reads}."\t";
    print STDOUT $hash{$key}{non_contam_non_phix_reads_1}."\t";
    print STDOUT $hash{$key}{non_contam_non_phix_reads_2}."\t";
    print STDOUT $hash{$key}{reads_1_QC_passed}."\t";
    print STDOUT $hash{$key}{assembled_reads}."\t";
    print STDOUT $hash{$key}{assembled_reads_QC_passed}."\t";
    print STDOUT $hash{$key}{clustered_sequences}."\t";
    print STDOUT $hash{$key}{number_of_clusters}."\t";
    print STDOUT $hash{$key}{clustered_sequences_classified}."\n";

    #foreach($key2 (keys %{$hash{$key} })) {
    #     
    #}
}

