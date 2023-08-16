#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
my $usage=<<'ENDHERE';
NAME:
parseGFFForMetatranscriptomics.pl

PURPOSE:

INPUT:
--infile_gff <string>  : gff
--infile_fna <string>  : genome fasta file
--outfile_fna <string> : gene fasta nucl file
--outfile_gff <string> : parsed gff file
--outfile_faa <string> : gene fasta aa file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_gff, $infile_fna, $outfile_link, $outfile_fna, $outfile_gff);
my $verbose = 0;

GetOptions(
   'infile_gff=s'   => \$infile_gff,
   'infile_fna=s'   => \$infile_fna,
   'outfile_gff=s'  => \$outfile_gff,
   'outfile_link=s' => \$outfile_link,
   'outfile_fna=s'  => \$outfile_fna,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT_GFF, ">".$outfile_gff) or die "Can't open $outfile_gff\n";
open(OUT_FNA, ">".$outfile_fna) or die "Can't open $outfile_fna\n";
open(OUT_LINK, ">".$outfile_link) or die "Can't open $outfile_link\n";

my %hash;
my %hash_geneid_old_to_new;
my $counter = 1;
open(GFF, "<".$infile_gff) or die "Can't open $infile_gff\n";
while(<GFF>){
   chomp;
   if($_ =~ m/^#/){
      print OUT_GFF $_."\n";
      next;
   }
   next if($_ =~ m/^$/);
   # Replace all " characters with ' chars. HtSeq has a problem with " chars...
   $_ =~ s/\"/\'/g;
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];
   my $type = $row[2];
   if($type ne "gene"){next;}
   my $originalId = $row[0];
   my @field = split(/;/, $row[8]);
   my $mgmId = $field[0];
   #$mgmId =~ s/ID=/gene_id=/;
   $mgmId =~ s/ID=//;
   #$mgmId =~ s/=/_/;
   #$mgmId =~ s/gene_id_\d+_\d+/gene_id_$counter/;
   $originalId =~ s/\s/\./g;
   #$hash{$mgmId} = $originalId."=".$mgmId;
   $hash{$row[0]}{$mgmId}{contig_id} = $row[0];
   # Determine which is start and which is end.
   my $start; my $end;
   if($row[3] < $row[4]){
      $start = $row[3];
      $end = $row[4];
   }elsif($row[3] > $row[4]){
      $start = $row[4];
      $end = $row[3];
   }else{
      die("Something is wrong with gene coordinates...\n");
   }
   $hash{$row[0]}{$mgmId}{start} = $start;
   $hash{$row[0]}{$mgmId}{end} = $end;
   $hash{$row[0]}{$mgmId}{strand} = $row[6];;
   $hash_geneid_old_to_new{$field[0]} = $mgmId;

   # Then modify only gene_id in the gff output
   my $new_line = $_;
   #$new_line =~ s/ID=/gene_id_/;
   $new_line =~ s/ID=//;
   #$new_line =~ s/gene_id_\d+_\d+/gene_id_$counter/;
   print OUT_GFF $new_line."\n";
   print OUT_LINK $contig_id."\t".$start."\t".$end."\t".$mgmId."\n";

   $counter++;
}
close(GFF);
close(OUT_GFF);
close(OUT_LINK);

print STDERR Dumper(\%hash);

# Genome file in input
# extract genes.
my $ref_fasta_db = Iterator::FastaDb->new($infile_fna) or die("Unable to open Fasta file, $infile_fna\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    #$header =~ s/>//;
    my ($contig_id) = $header =~ m/>(\S+) /;
    print STDERR "Contig id ". $contig_id."\n";
    if(exists $hash{$contig_id}){
        #foreach my $contig_id (keys %hash) {
            foreach my $gene_id (keys %{ $hash{$contig_id} }) {
                print STDERR "Contig id:  ".$contig_id."     Gene id: ". $gene_id."\n";
                my $start = $hash{$contig_id}{$gene_id}{start};
                my $end   = $hash{$contig_id}{$gene_id}{end};
                my $gene_seq = substr $curr->seq, ($start - 1), (($end - $start) + 1);
                
                if($hash{$contig_id}{$gene_id}{strand} eq "+"){
                    print OUT_FNA ">".$gene_id."\n".$gene_seq."\n";
                }elsif($hash{$contig_id}{$gene_id}{strand} eq "-"){
                    print OUT_FNA ">".$gene_id."\n".reverse_complement_IUPAC($gene_seq)."\n";
                }else{
                    die("problem with gene orientation on strand...\n");
                }
            }
         #}
    }
}
close(OUT_FNA);

sub reverse_complement_IUPAC {
   my $dna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $revcomp;
}

