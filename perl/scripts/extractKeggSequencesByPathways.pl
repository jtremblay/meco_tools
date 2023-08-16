#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
extractKeggSequencesByPathways.pl

PURPOSE:

INPUT:
--ko <string>       : ko main file
--outdir <string>   : outdir

OUTPUT:
STDOUT


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $outdir, $ko);
my $verbose = 0;

GetOptions(
   'outdir=s' => \$outdir,
   'ko=s'     => \$ko,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;


my %hash_ko;
my $counter = 0;
open(IN, "<".$ko) or die "Can't open $ko\n";
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
         $hash_ko{$currKO}{ENTRY} = $1;
      }
   }
   
   if($currCat eq "DEFINITION"){
      if($_ =~ m/DEFINITION\s+(.*)\s\[.*\]$/){
          $hash_ko{$currKO}{DEFINITION} = $1;
      }elsif($_ =~ m/DEFINITION\s+(.*)/){
          $hash_ko{$currKO}{DEFINITION} = $1;
      }
   }

   if($currCat eq "PATHWAY"){
      if($_ =~ m/^.*(ko\d{5})\s+(\S+).*$/){
          if(exists $hash_ko{$currKO}{PATHWAY_ID}){
            $hash_ko{$currKO}{PATHWAY_ID} .= "==".$1;
          }else{
            $hash_ko{$currKO}{PATHWAY_ID} = $1;
          }
          if(exists $hash_ko{$currKO}{PATHWAY_DESC}){
            $hash_ko{$currKO}{PATHWAY_DESC} .= "==".$2;
          }else{
            $hash_ko{$currKO}{PATHWAY_DESC} = $2;
          }
      }
   }

   if($currCat eq "MODULE"){
      if($_ =~ m/^.*(M\d{5})\s+(\S+).*$/){
          if(exists $hash_ko{$currKO}{MODULE_ID}){
            $hash_ko{$currKO}{MODULE_ID} .= "==".$1;
          }else{
            $hash_ko{$currKO}{MODULE_ID} = $1;
          }
          if(exists $hash_ko{$currKO}{MODULE_DESC}){
            $hash_ko{$currKO}{MODULE_DESC} .= "==".$2;
          }else{
            $hash_ko{$currKO}{MODULE_DESC} = $2;
          }
      }
   }
   
   if($currCat eq "GENES"){
      if($_ =~ m/^.*(ko\d{5})\s+(\S+).*$/){
          if(exists $hash_ko{$currKO}{PATHWAY_ID}){
            $hash_ko{$currKO}{PATHWAY_ID} .= "==".$1;
          }else{
            $hash_ko{$currKO}{PATHWAY_ID} = $1;
          }
          if(exists $hash_ko{$currKO}{PATHWAY_DESC}){
            $hash_ko{$currKO}{PATHWAY_DESC} .= "==".$2;
          }else{
            $hash_ko{$currKO}{PATHWAY_DESC} = $2;
          }
      }
   }


   if($_ =~ m/^\\\\\\$/){
      $currCat = "";
      $counter++;
   }
   
   #last if($counter == 20);

}
close(IN);

#print STDERR Dumper(\%hash_ko);

# Print final table
print STDOUT "#query\tkegg_gene_id\tKO_id\tENTRY\tDEFINITION\tPATHWAY_ID\tPATHWAY_DESC\tMODULE_ID\tMODULE_DESC\n";
for my $geneId (keys %hash){

   if(exists $hash_link{$geneId}){
      print STDOUT $hash{$geneId}."\t".$geneId;
      my $ko = $hash_link{$geneId};
      print STDOUT "\t".$ko;
      if(exists $hash_ko{$ko}{ENTRY}){        print STDOUT "\t".$hash_ko{$ko}{ENTRY};            }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{DEFINITION}){   print STDOUT "\t".$hash_ko{$ko}{DEFINITION};       }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{PATHWAY_ID}){   print STDOUT "\t".$hash_ko{$ko}{PATHWAY_ID};       }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{PATHWAY_DESC}){ print STDOUT "\t".$hash_ko{$ko}{PATHWAY_DESC};     }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{MODULE_ID}){    print STDOUT "\t".$hash_ko{$ko}{MODULE_ID};        }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{MODULE_DESC}){  print STDOUT "\t".$hash_ko{$ko}{MODULE_DESC}."\n"; }else{ print STDOUT "\tNULL\n"; }
   }else{
       #print STDERR "ko with no id: ".$geneId."\n";
   }
}

exit;

