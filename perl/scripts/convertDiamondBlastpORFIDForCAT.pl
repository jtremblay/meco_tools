#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
convertDiamondBlastpORFIDForCAT.pl

PURPOSE:

INPUT:
--infile_blastp <string> : blastp diamond output
--infile_gff <string>    : gff file given by prodigal
--infile_faa <string>    : faa fasta file

OUTPUT:
<STDOUT>                 : renamed blastp table
--outfile_fasta <string> : renamed faa sequences


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_blastp, $infile_gff, $infile_faa, $outfile_faa);
my $verbose = 0;

GetOptions(
   'infile_gff=s'    => \$infile_gff,
   'infile_blastp=s' => \$infile_blastp,
   'infile_faa=s'    => \$infile_faa,
   'outfile_faa=s'   => \$outfile_faa,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## MAIN
if($infile_faa){
    open(OUT, ">".$outfile_faa) or die "Can't open $outfile_faa for writing...\n";
}
my %hash;
my $last_contig_id = "";
my $last_gene_id = "";
my $i = 1;
open(IN, "<".$infile_gff) or die "Can't open $infile_gff\n";
while(<IN>){
    chomp;
    next if($_ =~ /^#/);
    my @row = split(/\t/, $_);
    my $curr_contig_id = $row[0];
    my @field = split(/;/, $row[8]);
    my $curr_gene_id = $field[0];

    if(eof){
        $i++;
        #print STDERR $curr_contig_id."_".$i."\t".$curr_gene_id."\n"; 
        $hash{$curr_gene_id} = $curr_contig_id."_".$i;
        
    }elsif($curr_contig_id eq $last_contig_id){
        #continue
        $i++;
        #print STDERR $curr_contig_id."_".$i."\t".$curr_gene_id."\n"; 
        $hash{$curr_gene_id} = $curr_contig_id."_".$i;

    }elsif($curr_contig_id ne $last_contig_id){
        # mean new contig.
        # 1st: complete last contig
        if($i == 1){
            #print STDERR $curr_contig_id."_".$i."\t".$curr_gene_id."\n";
            $hash{$curr_gene_id} = $curr_contig_id."_".$i;
        }else{
            $i = 1;
            #print STDERR $curr_contig_id."_".$i."\t".$curr_gene_id."\n";
            $hash{$curr_gene_id} = $curr_contig_id."_".$i;
        } 

        # 2nd start new one.
        $last_contig_id = $curr_contig_id;
        $last_gene_id = $curr_gene_id;
    }
}
close(IN);


open(IN, "<".$infile_blastp) or die "Can't open $infile_blastp\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $gene_id = shift(@row);
    if(exists $hash{$gene_id}){
        print STDOUT $hash{$gene_id}."\t".join("\t", @row)."\n";
    }else{
        print STDERR "$gene_id does not exist...\n";
    }
}
close(IN);

if($infile_faa){
    my $ref_fasta_db = Iterator::FastaDb->new($infile_faa) or die("Unable to open Fasta file, $infile_faa\n");
    while( my $curr = $ref_fasta_db->next_seq() ) {
        #my $header = $curr->header;
        my ($header) = $curr->header =~ m/^>(\S+)/;
        if(exists $hash{$header}){
           print OUT $header."\n".$curr->seq."\n";
        }else{
            print STDERR "$header does not exist...\n";
        }
    }
    close(OUT);
}

exit;
