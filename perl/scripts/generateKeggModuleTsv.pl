#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateKeggModuleTsv.pl

PURPOSE:

INPUT:
--ko <string>   : ko main file. Module version of the ko file.

OUTPUT:
STDOUT          : ko main file in tsv format. Mainly for easier
                  manipulation with R.


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $ko);
my $verbose = 0;

GetOptions(
   'ko=s'       => \$ko,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;

my $counter = 0;
open(IN, "<".$ko) or die "Can't open $ko\n";
my $currMO;
my $currCat;
while(<IN>){
   chomp;

   if($_ =~ m/^ENTRY/){
      $currCat = "ENTRY";
   }elsif($_ =~ m/^NAME/){
      $currCat = "NAME";
   }elsif($_ =~ m/^DEFINITION/){
      $currCat = "DEFINITION";    
   }elsif($_ =~ m/^ORTHOLOGY/){
      $currCat = "ORTHOLOGY";
    }   
    #print STDERR $currCat."\t".$_."\n";

   if($currCat eq "ENTRY"){
      if($_ =~ m/ENTRY\s+(M\d{5})/){
         $currMO = $1;
         $hash{$currMO}{ENTRY} = $1;
      }
   }
   
   if($currCat eq "NAME"){
      if($_ =~ m/NAME\s+(.*)$/){
         $hash{$currMO}{NAME} = $1;
      }
   }
   
   #if($currCat eq "DEFINITION"){
   #   if($_ =~ m/DEFINITION\s+(.*)\s\[.*\]$/){
   #       $hash{$currMO}{DEFINITION} = $1;
   #   }elsif($_ =~ m/DEFINITION\s+(.*)/){
   #       $hash{$currMO}{DEFINITION} = $1;
   #   }
   #}

   # K00844,K12407,K00845  hexokinase/glucokinase [EC:2.7.1.1 2.7.1.2] [RN:R01786]
   if($currCat eq "ORTHOLOGY"){
      my ($Ks, $EC, $RN);
      if($_ =~ m/\s+(K.*) /){
        $Ks = $1;
        my @curr_Ks;
        if($Ks =~ m/K\d{5}\,/){
            my @Ks = split(/,/, $Ks);
            foreach my $K (@Ks){
                my($curr_Ks) = $K =~ m/(K\d+)/;
                $Ks = $curr_Ks;
                push(@curr_Ks, $curr_Ks);
            }
        }elsif($Ks =~ m/K\d{5}\+/){
            my @Ks = split(/\+/, $Ks);
            foreach my $K (@Ks){
                my($curr_Ks) = $K =~ m/(K\d+)/;
                $Ks = $curr_Ks;
                push(@curr_Ks, $curr_Ks);
            }
        }elsif($Ks =~ m/K\d{5}\ /){
            my @Ks = split(/\s/, $Ks);
            foreach my $K (@Ks){
                my($curr_Ks) = $K =~ m/(K\d+)/;
                $Ks = $curr_Ks;
                push(@curr_Ks, $curr_Ks);
            }
        }else{
            $Ks = $Ks
        }

        $Ks = join(",", @curr_Ks);
        $Ks =~ s/\,\,//g;
        $Ks =~ s/\,\,\,//g;
        $Ks =~ s/\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,\,\,\,\,\,//g;
        $Ks =~ s/\,\,\,\,\,\,\,\,\,\,\,//g;

        if($_ =~ m/.*\[(EC:\S+)\]\s\[(RN:R\d+)\]$/){
            $EC = $1;
            $RN = $2;
            
            $hash{$currMO}{$RN}{Ks} = $Ks;
            $hash{$currMO}{$RN}{EC} = $EC;
            $hash{$currMO}{$RN}{RN} = $RN;
        }
      }
   }

   #if($_ =~ m/^\\\\\\$/){
   #   $currCat = "";
   #   $counter++;
   #}
   
   #last if($counter == 1000);

}
close(IN);

print STDERR Dumper(\%hash);

print STDOUT "#MODULE\tMODULE_NAME\tREACTION\tENZYME_CODE\tKEGG_ENTRY\n";

foreach my $currMO (keys %hash) {
    my $module_name = $hash{$currMO}{NAME};
    foreach my $RN (keys %{ $hash{$currMO} }) {
        if($RN =~ m/RN:/){
            my $Ks = $hash{$currMO}{$RN}->{Ks};
            my @Ks = split(/,/, $Ks);
            foreach my $K (@Ks){
                print STDOUT "$currMO\t$module_name\t$RN\t$hash{$currMO}{$RN}{EC}\t$K\n";
            }
        }
    }
}
#
#   if(exists $hash_link{$kegg_gene_id}){
#       #print STDOUT $hash{$geneId}."\t".$geneId;
#      print STDOUT $geneId."\t".$kegg_gene_id;
#      my $ko = $hash_link{$kegg_gene_id};
#      print STDOUT "\t".$ko;
#      if(exists $hash_ko{$ko}{NAME}){         print STDOUT "\t".$hash_ko{$ko}{NAME};             }else{ print STDOUT "\tNULL";   }
#      if(exists $hash_ko{$ko}{ENTRY}){        print STDOUT "\t".$hash_ko{$ko}{ENTRY};            }else{ print STDOUT "\tNULL";   }
#      if(exists $hash_ko{$ko}{DEFINITION}){   print STDOUT "\t".$hash_ko{$ko}{DEFINITION};       }else{ print STDOUT "\tNULL";   }
#      if(exists $hash_ko{$ko}{PATHWAY_ID}){   print STDOUT "\t".$hash_ko{$ko}{PATHWAY_ID};       }else{ print STDOUT "\tNULL";   }
#      if(exists $hash_ko{$ko}{PATHWAY_DESC}){ print STDOUT "\t".$hash_ko{$ko}{PATHWAY_DESC};     }else{ print STDOUT "\tNULL";   }
#      if(exists $hash_ko{$ko}{MODULE_ID}){    print STDOUT "\t".$hash_ko{$ko}{MODULE_ID};        }else{ print STDOUT "\tNULL";   }
#      if(exists $hash_ko{$ko}{MODULE_DESC}){  print STDOUT "\t".$hash_ko{$ko}{MODULE_DESC}."\n"; }else{ print STDOUT "\tNULL\n"; }
#   }else{
#       print STDERR "ko with no id: ".$geneId." => ".$kegg_gene_id."\n";
#       $hash_second_pass{$geneId} = $kegg_gene_id;
#   }
#}
#print STDERR Dumper(\%hash_link_noprefix);
#exit;

exit;

