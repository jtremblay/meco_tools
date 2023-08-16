#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
#use Bio::SeqIO;
use Bio::DB::Query::GenBank;
use Bio::DB::GenBank;
#use lib "/house/homedirs/j/jtremblay/myperl/lib/";
#use FastaDb;

my $usage=<<'ENDHERE';
NAME:
get_genbank_taxonomy.pl

PURPOSE:
From a list of GI entries, fetch taxonomy lineage from NCBI

INPUT:
--infile <fasta_file> : fasta file containing only one primer seq

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.comjtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $infile, $outfile);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'infile=s'  => \$infile,
    'verbose'   => \$verbose,
    'help'      => \$help
);
if ($help) { print $usage; exit; }

#VALIDATE

#MAIN=====================================================================================================================================================================

my %seen;

my @queries;

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $geneId = $row[0];
   #my $query = $row[1];
   my $query = $row[0];

   if(!exists $seen{$query}){
      push(@queries, $query);  
   }
}
close(IN);

print STDERR join("\n", @queries);
my $i = 0;

#for my $query (@queries){
my $query = join(",", @queries);
print STDERR "Processing current query: $query\n";
#my $query = "CP000828.1";
my $query_obj = Bio::DB::Query::GenBank->new(-db => 'nucleotide', -query => $query );

my $gb_obj = Bio::DB::GenBank->new;
 
my $stream_obj = $gb_obj->get_Stream_by_query($query_obj);
 

while (my $seq_object = $stream_obj->next_seq) {  
	   
   my $species_object = $seq_object->species;
   my $species_string = $species_object->node_name;
    
   # Perlish
   #my $species_string = $seq_object->species->node_name;
   # either way, $species_string is "Homo sapiens"
    
   # get all taxa from the ORGANISM section in an array
    #print STDOUT "$query\t";
	my @classification = $seq_object->species->classification;
    @classification = reverse(@classification);
	foreach $_ (@classification){
		print STDOUT $_."\t";
	}
    print STDOUT "\n";
}   
#}

exit;
