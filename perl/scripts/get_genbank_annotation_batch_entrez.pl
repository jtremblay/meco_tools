#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Bio::DB::Query::GenBank;
use Bio::DB::GenBank;
use Iterator::FastaDb;
use Data::Dumper;
use LWP::Simple;

my $usage=<<'ENDHERE';
NAME:
get_genbank_annotation.pl

PURPOSE:
From a GI, get taxonomic info from ncbi

INPUT:
--infile <string> : list file having one reference per line.
--fasta <string>  : fasta file 

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $fasta);
my $verbose = 0;

## SCRIPTS
GetOptions(
   'infile=s'  => \$infile,
   'fasta=s'   => \$fasta,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

#VALIDATE
die "--infile missing\n" unless($infile);

#MAIN

#my $query = "255342900";
#my $query2 = "190714214";
#my @query = (255342900, 190714214);
#my $query_obj = Bio::DB::Query::GenBank->new(-db => 'nucleotide',  -query => join(" ", @query));

#my %seen;
#my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
#while( my $curr = $ref_fasta_db->next_seq() ) {
#   my $header = $curr->header;
#   print STDERR "HEADER:$header\n";
#   my ($id) = $header =~ m/^>(\S\+) /;
#   $seen{$id} = "";
#}

#print STDERR Dumper(\%seen);
#exit;

#my @queries;

open(IN, "<".$infile) or die "Can't open $infile\n";
my $i = 0;
my @queries;
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $geneId = $row[0];
   #my $query = $row[1];
   my $query = $row[0];
   push(@queries, $query);
}
close(IN);
my $id_list = join(",", @queries);

# Download protein records corresponding to a list of GI numbers.

my $db = 'protein';
#$id_list = '194680922,50978626,28558982,9507199,6678417';

#assemble the epost URL
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "epost.fcgi?db=$db&id=$id_list";

#post the epost URL
my $output = get($url);

#parse WebEnv and QueryKey
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

### include this code for EPost-ESummary
#assemble the esummary URL
$url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";

#post the esummary URL
my $docsums = get($url);
print "$docsums";

### include this code for EPost-EFetch
#assemble the efetch URL
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
$url .= "&rettype=fasta&retmode=text";

#post the efetch URL
my $data = get($url);
print STDOUT "$data";



exit;
