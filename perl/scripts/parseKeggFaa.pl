#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseKeggFaa.pl

PURPOSE:

INPUT:
--infile_fasta <string>     : Sequence file. Kegg pep faa file.
--infile_ko <string>        : ko ref file.
				
OUTPUT:
--outfile_ko_terms <string> : outfile that links ko and terms.
STDOUT
parsed faa.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_ko, $outfile_ko_terms);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'    => \$infile_fasta,
   'infile_ko=s'       => \$infile_ko,
   'outfile_ko_terms=s'=> \$outfile_ko_terms,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash_ko;
my %hash_link;
my $counter = 0;
my $i = 0;
open(IN, "<".$infile_ko) or die "Can't open $infile_ko\n";
open(OUT, ">".$outfile_ko_terms) or die "Can't open $outfile_ko_terms\n";
my $currKO;
my $currCat;
while(<IN>){
   chomp;

   if($_ =~ m/^ENTRY/){
      $currCat = "ENTRY";
   }elsif($_ =~ m/^DEFINITION/){
      $currCat = "DEFINITION";    
   }elsif($_ =~ m/^PATHWAY/){
      $currCat = "PATHWAY";
   }elsif($_ =~ m/^MODULE/){
      $currCat = "MODULE";
   }elsif($_ =~ m/DBLINKS/){
      $currCat = "DBLINKS";
   }elsif($_ =~ m/GENES/){
      $currCat = "GENES";
   }elsif($_ =~ m/BRITE/){
      $currCat = "BRITE";
   }
   #print STDERR $currCat."\t".$_."\n";

   if($currCat eq "ENTRY"){
      if($_ =~ m/ENTRY\s+(K\d{5})/){
         $currKO = $1;
         #$hash_link{$currKO}{ENTRY} = $1;
      }
   }
#   
#   if($currCat eq "DEFINITION"){
#      if($_ =~ m/DEFINITION\s+(.*)\s\[.*\]$/){
#          $hash_ko{$currKO}{DEFINITION} = $1;
#      }elsif($_ =~ m/DEFINITION\s+(.*)/){
#          $hash_ko{$currKO}{DEFINITION} = $1;
#      }
#   }
#
#   if($currCat eq "PATHWAY"){
#      if($_ =~ m/^.*(ko\d{5})\s+(\S+).*$/){
#          if(exists $hash_ko{$currKO}{PATHWAY_ID}){
#            $hash_ko{$currKO}{PATHWAY_ID} .= "==".$1;
#          }else{
#            $hash_ko{$currKO}{PATHWAY_ID} = $1;
#          }
#          if(exists $hash_ko{$currKO}{PATHWAY_DESC}){
#            $hash_ko{$currKO}{PATHWAY_DESC} .= "==".$2;
#          }else{
#            $hash_ko{$currKO}{PATHWAY_DESC} = $2;
#          }
#      }
#   }
#
#   if($currCat eq "MODULE"){
#      if($_ =~ m/^.*(M\d{5})\s+(\S+).*$/){
#          if(exists $hash_ko{$currKO}{MODULE_ID}){
#            $hash_ko{$currKO}{MODULE_ID} .= "==".$1;
#          }else{
#            $hash_ko{$currKO}{MODULE_ID} = $1;
#          }
#          if(exists $hash_ko{$currKO}{MODULE_DESC}){
#            $hash_ko{$currKO}{MODULE_DESC} .= "==".$2;
#          }else{
#            $hash_ko{$currKO}{MODULE_DESC} = $2;
#          }
#      }
#   }
   
   if($currCat eq "GENES"){
       #print STDERR "Found GENES\n";
      if($_ =~ m/^\s+(\S{4}:)\s(.*)$/){
         #print STDERR "Found $1 and $2\n";
         my $prefix = lc($1);
         my @row = split(/\s/, $2);
         
         foreach my $el (@row){
            $el =~ s/\(.*\)//g;
            $hash_ko{$prefix."".$el} = $prefix."".$el;
            $hash_link{$prefix."".$el} = $currKO;
            #if(exists $hash_link{$currKO}{GENES}){
            #   $hash_link{$currKO}{GENES} .= "==".$2;
            #}else{
            #   $hash_link{$currKO}{GENES} = $2;
            #}
         }
      }elsif($_ =~ m/^\s+(\S{3}:)\s(.*)$/){
         #print STDERR "Found $1 and $2\n";
         my $prefix = lc($1);
         my @row = split(/\s/, $2);
         
         foreach my $el (@row){
            $el =~ s/\(.*\)//g;
            $hash_ko{$prefix."".$el} = $prefix."".$el;
            $hash_link{$prefix."".$el} = $currKO;
         }

      }elsif($_ =~ m/^\s+(\S{2}:)\s(.*)$/){
         #print STDERR "Found $1 and $2\n";
         my $prefix = lc($1);
         my @row = split(/\s/, $2);
         
         foreach my $el (@row){
            $el =~ s/\(.*\)//g;
            $hash_ko{$prefix."".$el} = $prefix."".$el;
            $hash_link{$prefix."".$el} = $currKO;
         }
      }
   }

   if($_ =~ m/^\\\\\\$/){
      $currCat = "";
      $counter++;
      print STDERR "counter: ".$counter."\n";
   }
      
   #last if($counter >= 3);
   $i++;
   #last if($i >= 2000);

}
close(IN);

print STDERR Dumper(\%hash_ko);
#print STDERR Dumper(\%hash_link);
#exit;
print OUT "kegg_gene_id\tko_term\tsource\n";
for my $key (keys %hash_link){
   print OUT $key."\tKO:".$hash_link{$key}."\t\n";
}
close(OUT);
#exit;

$counter = 0;
my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $header = $curr->header();
   my($gene_id) = $header =~ /^(>\S+)/;
   $gene_id =~ s/>//;
   #print $gene_id."\n";
   if(exists $hash_ko{$gene_id}){
      print STDOUT $header." symbol=$gene_id\n".$curr->seq."\n";
   }else{
     #print STDERR $header."\n";
     print STDERR $header." symbol=$gene_id\n".$curr->seq."\n";
   }
   #$counter++;
   #last if($counter >= 100000);
}

