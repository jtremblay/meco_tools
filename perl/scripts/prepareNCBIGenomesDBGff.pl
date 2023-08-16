#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
prepareNCBIGenomesDB.pl

PURPOSE:
To setup the database containing all 'reference' and 'representative' genomes found here:
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
where column 4 is labeled with reference or representative (and more).

INPUT:
--infile <string> : Assembly summary file.
--outdir <string> : Directory were genomes files will be stored.
                    Each downloaded file will be merged to only 
                    get one 'big' file containing all relevant genomes
--concatenate     : Set flag if all genomes are to be concatenated into one single file.

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outdir, $concatenate);
my $verbose = 0;

GetOptions(
   'infile=s'      => \$infile,
   'outdir=s'      => \$outdir,  
   'concatenate'   => \$concatenate,
   'verbose'       => \$verbose,
   'help'          => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my $i = 1;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;

    next if($_ =~ m/^#/);

    my @row = split(/\t/, $_);
    if($row[4] eq "representative genome" || $row[4] eq "reference genome"){
        
        my $basename = $row[19];
        $basename =~ s{.*/}{};         # removes path  
        #$basename =~ s{\.[^.]+$}{};    # removes extension
        my $filename = $basename."_genomic.fna.gz";
       
        print STDERR $basename."\n";

        my $genome_url = $row[19]."/".$basename."_genomic.fna.gz";

        # Get the file
        system("wget -P $outdir $genome_url");
        die "command failed: $!\n" if($? != 0);


        if($concatenate){
            # put the file into the large concatenated DB.
            if($i == 1){
                system("mv $outdir/$filename $outdir/ncbi_genomes.fna.gz");
                die "command failed: $!\n" if($? != 0);
            }else{
                system("cat $outdir/$filename >> $outdir/ncbi_genomes.fna.gz"); 
                die "command failed: $!\n" if($? != 0);
                system("rm $outdir/$filename");
                die "command failed: $!\n" if($? != 0);
            }
            
            $i++;
        }
    }
}
