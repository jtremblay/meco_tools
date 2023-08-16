#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
cleanMAGBasedOnTaxonomy.pl

PURPOSE:
Takes indirectory of metabin and generates a tsv file with each bins + percentages of each 
different order in each bins. If --split_bins_by_taxonomy <string> is set, only to up to
3 most abundant bins of that particular <string> will be kept.

INPUT:
--infile_fasta <string>    : MAG/bin multi fasta file
--infile_taxonomy <string> : taxonomy.tsv file (contigs tax file.) given by CAT
--split <string>           : Either 'k__', 'p__', 'c__', 'o__', 'f__', 'g__' or 's__' (default: 'o__')
--taxon <string>           : A taxon. Ex: 'Oceanospirillales'

OUTPUT:
--outfile_fasta <string>   : Cleaned MAG multi-fasta outfile.
--outfile_stats <string>   : Stats. how many contigs were left out vs kept. 
STDOUT <string>            : Cleaned taxonomy file. 

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_taxonomy, $split, $outfile_stats, $taxon, $outfile_fasta);
my $verbose = 0;

GetOptions(
   'infile_taxonomy=s'  => \$infile_taxonomy,
   'infile_fasta=s'     => \$infile_fasta,
   'taxon=s'            => \$taxon,
   'outfile_stats=s'    => \$outfile_stats,
   'outfile_fasta=s'    => \$outfile_fasta,
   'split=s'            => \$split,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile_fasta is missing\n" unless($infile_fasta);
die "--infile_taxonomy is missing\n" unless($infile_taxonomy);
die "--outfile_stats is missing\n" unless($outfile_stats);


## MAIN
open(OUT_F, ">".$outfile_fasta) or die "Can't open $outfile_fasta\n";
open(OUT_S, ">".$outfile_stats) or die "Can't open $outfile_stats\n";

# Here, if split = 'o__', split at the order level and --taxon Oceanospirillales. All lineage who matched Oceanospirillales will be kept.
# If taxon at 'o__' does not match Oceanospirillales and does not match undef, it is discarded.
# If taxon at 'o__' matches undef, contig is kept only if Kingdom also equals Bacteria.
my %hash;
my $rejected_because_undef_kingdom = 0;
my $rejected_because_not_matching_split = 0;
my $kept = 0;
my $total_seqs = 0;
open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $contig_id = $row[0];
    #my $species = join(".", $row[8]);
    my $kingdom = $row[1];
    my $phylum = $row[2];
    my $class = $row[3];
    my $order = $row[4];
    my $family = $row[5];
    my $genus = $row[6];
    my $species = $row[7];
   
    my $depth_to_compare;
    if($split eq "o__"){
        $depth_to_compare = $order;
    }elsif($split eq "c__"){
        $depth_to_compare = $order;
    }else{
        die "only o__ and c__ implemented\n";
    }

    if($depth_to_compare eq $taxon){
        $hash{$contig_id}{kingdom} = $kingdom;
        $hash{$contig_id}{phylum} = $phylum;
        $hash{$contig_id}{class} = $class;
        $hash{$contig_id}{order} = $order;
        $hash{$contig_id}{family} = $family;
        $hash{$contig_id}{genus} = $genus;
        $hash{$contig_id}{species} = $species;
        $kept++;
        print STDOUT $_."\n";
    }elsif($depth_to_compare eq "undef" && $kingdom ne "undef"){
        $hash{$contig_id}{kingdom} = $kingdom;
        $hash{$contig_id}{phylum} = $phylum;
        $hash{$contig_id}{class} = $class;
        $hash{$contig_id}{order} = $order;
        $hash{$contig_id}{family} = $family;
        $hash{$contig_id}{genus} = $genus;
        $hash{$contig_id}{species} = $species;
        $kept++;
        print STDOUT $_."\n";
    }elsif($kingdom eq "undef"){
        $rejected_because_undef_kingdom++;
        print STDERR "Rejecting: ".$_."\n";
    }else{
        $rejected_because_not_matching_split++; 
        print STDERR "Rejecting: ".$_."\n";
    }
    $total_seqs++;
}
close(IN);


# Parse fasta.
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my ($contig_id) = $curr->header =~ /^>(\S+)/;
   $contig_id =~ s/^>//;
   
   if(exists $hash{$contig_id}){
       #print OUT_F $contig_id."\t".$hash{$contig_id}{kingdom}."\t".$hash{$contig_id}{phylum}."\t".$hash{$contig_id}{class}."\t".$hash{$contig_id}{order}."\t".$hash{$contig_id}{family}."\t".$hash{$contig_id}{genus}."\t".$hash{$contig_id}{species}."\n";
       print OUT_F $curr->output;
   }
}

print OUT_S "Total seqs: $total_seqs\n";
print OUT_S "Kept sequences: $kept\n";
print OUT_S "Sequences rejected because undef kingdom: $rejected_because_undef_kingdom\n";
print OUT_S "Sequenced rejected because $split not matching $taxon : $rejected_because_not_matching_split\n";
close(OUT_S);
close(OUT_F);
exit;
