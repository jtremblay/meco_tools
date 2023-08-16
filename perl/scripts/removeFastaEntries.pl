#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
script_template.pl

PURPOSE:

INPUT:
--infile <string>    : Sequence file
--to_remove <string> : terms to remove. Will search for regular expression (not exact match).
				       one term per row.
OUTPUT:
<STDOUT>             : input fasta file minus matched fasta entries for which the headers
                       matched the remove list.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $to_remove);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'to_remove=s' => \$to_remove,
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my @terms;
open(IN, "<".$to_remove) or die "Can't open $to_remove\n";
while(<IN>){
   chomp;
   push(@terms, $_ );
}
close(IN);

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $match = "false";
    foreach my $term (@terms){
        if($curr->header =~ m/$term/){
            $match = "true";
        }
    }
    if($match eq "false"){
        print STDOUT $curr->header."\n".$curr->seq."\n";
    }
}

exit;
