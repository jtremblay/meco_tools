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
--infile <string>         : output of kofmascan or hmmsearch against kofamscan or hmmsearch. 
--ref_database <string>   : KEGG reference database in tsv format. 
--hmmsearch               ; set flag if hmmsearch was used instead of kofamscan

OUTPUT:
STDOUT


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $ref_database, $hmmsearch);
my $verbose = 0;

GetOptions(
   'infile=s'       => \$infile,
   'ref_database=s' => \$ref_database,
   'hmmsearch'      => \$hmmsearch,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN

# Parse kofamscan table
my %hash;
my %hash_KO;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    $_ =~ s/^\s+//;
    my @row = split(/\s+/, $_);
    my $gene_id;
    my $KO;
    my $evalue;
    if($hmmsearch){
        $gene_id = $row[2];
        $KO = $row[0];
        $evalue = $row[4];
    }else{
        $gene_id = $row[0];
        $KO = $row[1];
        $evalue = $row[4];
    }
    if($evalue <= 1e-10){
        if(exists $hash{$gene_id}){
            if($evalue < $hash{$gene_id}{evalue}){
                $hash{$gene_id}{KO} = $KO;
                $hash{$gene_id}{evalue} = $evalue;
            }
        }else{
            $hash{$gene_id}{KO} = $KO;
            $hash{$gene_id}{evalue} = $evalue;
        }
    }
}
close(IN);

open(IN, "<".$ref_database) or die "Can't open $ref_database\n";
while(<IN>){
    #chomp;
    $_ =~ s/$//;
    my @row = split(/\t/, $_);
    my $KO = $row[0];
    my $KO_desc = $row[1];
    my $plevel1 = $row[2];
    my $plevel1_desc = $row[3];
    my $plevel2 = $row[4];
    my $plevel2_desc = $row[5];
    my $pathway = $row[6];
    my $pathway_desc = $row[7];
    my $mlevel1_desc = $row[8];
    my $mlevel2_desc = $row[9];
    my $mlevel3_desc = $row[10];
    my $module = $row[11];
    my $module_desc = $row[12];
    
    if($KO_desc =~ m/\|/){$KO_desc =~ s/\|/==/g;}
    if($plevel1 =~ m/\|/){$plevel1 =~ s/\|/==/g;}
    if($plevel1_desc =~ m/\|/){$plevel1_desc =~ s/\|/==/g;}
    if($plevel2 =~ m/\|/){$plevel2 =~ s/\|/==/g;}
    if($plevel2_desc =~ m/\|/){$plevel2_desc =~ s/\|/==/g;}
    if($pathway =~ m/\|/){$pathway =~ s/\|/==/g;}
    if($pathway_desc =~ m/\|/){$pathway_desc =~ s/\|/==/g;}
    if($mlevel1_desc =~ m/\|/){$mlevel1_desc =~ s/\|/==/g;}
    if($mlevel2_desc =~ m/\|/){$mlevel2_desc =~ s/\|/==/g;}
    if($mlevel3_desc =~ m/\|/){$mlevel3_desc =~ s/\|/==/g;}
    if($module =~ m/\|/){$module =~ s/\|/==/g;}
    if($module_desc =~ m/\|/){$module_desc =~ s/\|/==/g;}
    chomp($module_desc);

    my $symbol = "";
    my $KO_desc2;
    $KO_desc =~ s/^ //;
    if($KO_desc =~ m/^(.*); (.*)$/){
        $symbol = $1;
        $KO_desc2 = $2;
    }

    $hash_KO{$KO}{NAME} = $symbol;
    $hash_KO{$KO}{DEFINITION} = $KO_desc2;
    $hash_KO{$KO}{ENTRY} = $KO;
    $hash_KO{$KO}{PATHWAY_ID} = $pathway;
    $hash_KO{$KO}{PATHWAY_DESC} = $pathway_desc;
    $hash_KO{$KO}{MODULE_ID} = $module;
    $hash_KO{$KO}{MODULE_DESC} = $module_desc;
}
close(IN);

# Then loop through KEGG ref file

#print STDERR Dumper(\%hash_KO);

# Print final table

print STDOUT "#query\tkegg_gene_id\tKO_id\tNAME\tENTRY\tDEFINITION\tPATHWAY_ID\tPATHWAY_DESC\tMODULE_ID\tMODULE_DESC\tevalue\n";
for my $gene_id (keys %hash){
   my $kegg_gene_id = ""; # gene_id is not available with kofamscan

   print STDOUT $gene_id."\t".$kegg_gene_id;
   my $KO = $hash{$gene_id}{KO};
   print STDOUT "\t".$KO;
   if(exists $hash_KO{$KO}{NAME}){         print STDOUT "\t".$hash_KO{$KO}{NAME};             }else{ print STDOUT "\tNULL";   }
   if(exists $hash_KO{$KO}{ENTRY}){        print STDOUT "\t".$hash_KO{$KO}{ENTRY};            }else{ print STDOUT "\tNULL";   }
   if(exists $hash_KO{$KO}{DEFINITION}){   print STDOUT "\t".$hash_KO{$KO}{DEFINITION};       }else{ print STDOUT "\tNULL";   }
   if(exists $hash_KO{$KO}{PATHWAY_ID}){   print STDOUT "\t".$hash_KO{$KO}{PATHWAY_ID};       }else{ print STDOUT "\tNULL";   }
   if(exists $hash_KO{$KO}{PATHWAY_DESC}){ print STDOUT "\t".$hash_KO{$KO}{PATHWAY_DESC};     }else{ print STDOUT "\tNULL";   }
   if(exists $hash_KO{$KO}{MODULE_ID}){    print STDOUT "\t".$hash_KO{$KO}{MODULE_ID};        }else{ print STDOUT "\tNULL";   }
   if(exists $hash_KO{$KO}{MODULE_DESC}){  print STDOUT "\t".$hash_KO{$KO}{MODULE_DESC}."\t"; }else{ print STDOUT "\tNULL\t"; }
   print STDOUT $hash{$gene_id}{evalue}."\n";
}
#print STDERR Dumper(\%hash_link_noprefix);

exit;

