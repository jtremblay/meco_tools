#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseKegg.pl

PURPOSE:

INPUT:
--ko                            : KEGG ko file.
--pathways <string>             : tsv file containing all three hiearchichal
                                  levels of KEGGs pathways.
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $ko, $pathways_file);
my $verbose = 0;

GetOptions(
   'ko=s' 	          => \$ko,
   'pathways_file=s'  => \$pathways_file,
   'verbose' 	       => \$verbose,
   'help'             => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %ref_hash;
my %hash;
my %hash_ko;

# loop through selected pathways (3 levels formated tsv)
open(IN, "<".$pathways_file) or die "Can't open $pathways_file\n";
while(<IN>){
   chomp;
   my @pathways = split(/\t/, $_);
   $ref_hash{$pathways[2]} = $pathways[2];
}
close(IN);


my $counter = 0;
open(IN, "<".$ko) or die "Can't open $ko\n";
my $currKO;
my $currCat;
while(<IN>){
   chomp;

   if($_ =~ m/^ENTRY/){
      $currCat = "ENTRY";
   }elsif($_ =~ m/^NAME/){
      $currCat = "NAME";
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
   
   if($currCat eq "NAME"){
      if($_ =~ m/NAME\s+(.*)$/){
         $hash_ko{$currKO}{NAME} = $1;
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
      if($_ =~ m/^.*(ko\d{5})\s+(.*)$/){
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
      if($_ =~ m/^.*(M\d{5})\s+(.*)$/){
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

   if($_ =~ m/^\\\\\\$/){
      $currCat = "";
      $counter++;
   }
   
   #last if($counter == 20);

}
close(IN);

#print STDOUT "#ENTRY\tNAME\tDEFINITION\tPATHWAY_DESC\tMODULE_DESC\n";
print STDOUT "#ENTRY\tNAME\tDEFINITION\n";
for my $key (keys %hash_ko){
   if(exists $hash_ko{$key}{PATHWAY_DESC}){
       my @pathways = split(/==/, $hash_ko{$key}{PATHWAY_DESC});
       foreach my $pathway (@pathways){
         if(exists $ref_hash{$pathway}){
            print STDOUT $hash_ko{$key}{ENTRY}."\t";
            print STDOUT $hash_ko{$key}{NAME}."\t";
            print STDOUT $hash_ko{$key}{DEFINITION}."\t";
            print STDOUT $hash_ko{$key}{PATHWAY_ID}."\t";
            print STDOUT $hash_ko{$key}{PATHWAY_DESC}."\t";
            print STDOUT $hash_ko{$key}{MODULE_DESC}."\n";
         }
      }
   }
}
exit;
