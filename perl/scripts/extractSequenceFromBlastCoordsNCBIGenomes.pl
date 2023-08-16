#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
extractSequenceFromBlastCoordsNCBIGenomes.pl

PURPOSE:

INPUT:
--infile_fasta <string>           : Sequence file (eg NCBI genome fasta file)
--infile_blastout <string>        : blast out results of blast against the--infile_fasta file.
--outfile_fasta_complete <string> : Fasta file with accession - complete record. 

OUTPUT:
STDOUT <string>                   : Fasta file with accession + coords in header.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_blastout, $outfile_fasta_complete);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'           => \$infile_fasta,
   'infile_blastout=s'        => \$infile_blastout,
   'outfile_fasta_complete=s' => \$outfile_fasta_complete,
   'verbose'                  => \$verbose,
   'help'                     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $OUT;
if($outfile_fasta_complete){
    open($OUT, ">".$outfile_fasta_complete) or die "Can't open $outfile_fasta_complete for writing\n";
}

my %hash;
open(IN, "<".$infile_blastout) or die "Can't open $infile_blastout\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   $hash{$row[1]}{start} = $row[8];
   $hash{$row[1]}{end} = $row[9];
}
close(IN);

print STDERR Dumper(\%hash);

my $i = 0;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my ($header) = $curr->header =~ m/^>(\S+)/;
   

    if(exists $hash{$header}){

        if($outfile_fasta_complete){
            print $OUT ">".$header."\n".$curr->seq."\n"
        }

        print STDERR "Found ".$header."\n";    
        my $start = $hash{$header}{start};
        my $end = $hash{$header}{end};

        # Orient in proper orientation.
        my $reverse;
        my $start2;
        my $end2;
        if($start > $end){
            $reverse = "true";
            $start2 = $end;
            $end2 = $start;
        }else{
            $reverse = "false";
            $start2 = $start;
            $end2 = $end;        
        }
        my $matched_seq = substr( $curr->seq, $start2, (($end2)-($start2)) );
        if($reverse eq "false"){
            print STDOUT ">".$header."\n".$matched_seq."\n";
        }else{
            print STDOUT ">".$header."\n".reverse($matched_seq)."\n";
        }
        
        #$i++;
    }#else{
    #   print STDERR $header." does not exist...\n";
    #}
}
exit;

