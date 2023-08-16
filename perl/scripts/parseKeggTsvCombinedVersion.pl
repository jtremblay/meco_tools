#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
parseKeggTsvCombinedVersion.pl

PURPOSE:

INPUT:
--infile <string> : kegg_ref_pathways_modules.tsv
				
OUTPUT:
<STDOUT>          : Combined version: One KO per line, but one or many pathways and/or modules per line.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' => \$infile,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;

    my @row = split(/\t/, $_);
    my $plevel1      = $row[0];
    my $plevel1_desc = $row[1];
    my $plevel2      = $row[2];
    my $plevel2_desc = $row[3];
    my $plevel3      = $row[4];
    my $plevel3_desc = $row[5];
    my $KO           = $row[6];
    my $KO_desc      = $row[7];
    my $mlevel1_desc = $row[8];
    my $mlevel2_desc = $row[9];
    my $mlevel3_desc = $row[10];
    my $module       = $row[11];
    my $module_desc  = $row[12];

    if(exists $hash{$KO}){
        $hash{$KO}{KO}            = $KO;
        $hash{$KO}{KO_desc}       = $KO_desc; 
        if($plevel3 eq "undefined"){ 
            $hash{$KO}{plevel1}       .= "|";
            $hash{$KO}{plevel1_desc}  .= "|";
            $hash{$KO}{plevel2}       .= "|";
            $hash{$KO}{plevel2_desc}  .= "|";
            $hash{$KO}{plevel3}       .= "|";
            $hash{$KO}{plevel3_desc}  .= "|";
        }else{
            if($hash{$KO}{plevel1} eq $plevel1 || $hash{$KO}{plevel1} eq ""){
                $hash{$KO}{plevel1}        = $plevel1;
            }else{
                $hash{$KO}{plevel1}       .= "|".$plevel1;
            }

            if($hash{$KO}{plevel1_desc} eq $plevel1_desc || $hash{$KO}{plevel1_desc} eq ""){
                $hash{$KO}{plevel1_desc}   = $plevel1_desc;
            }else{
                $hash{$KO}{plevel1_desc}  .= "|".$plevel1_desc;
            }

            if($hash{$KO}{plevel2} eq $plevel2 || $hash{$KO}{plevel2} eq ""){
                $hash{$KO}{plevel2}        = $plevel2;
            }else{
                $hash{$KO}{plevel2}       .= "|".$plevel2;
            }

            if($hash{$KO}{plevel2_desc} eq $plevel2_desc || $hash{$KO}{plevel2_desc} eq ""){
                $hash{$KO}{plevel2_desc}   = $plevel2_desc;
            }else{
                $hash{$KO}{plevel2_desc}  .= "|".$plevel2_desc;
            }

            if($hash{$KO}{plevel3} eq $plevel3 || $hash{$KO}{plevel3} eq ""){
                $hash{$KO}{plevel3}        = $plevel3;
            }else{
                $hash{$KO}{plevel3}       .= "|".$plevel3;
            }

            if($hash{$KO}{plevel3_desc} eq $plevel3_desc || $hash{$KO}{plevel3_desc} eq ""){
                $hash{$KO}{plevel3_desc}   = $plevel3_desc;
            }else{
                $hash{$KO}{plevel3_desc}  .= "|".$plevel3_desc;
            }
        }
        if($module eq "undefined"){
            $hash{$KO}{mlevel1_desc}  .= "|";
            $hash{$KO}{mlevel2_desc}  .= "|";
            $hash{$KO}{mlevel3_desc}  .= "|";
            $hash{$KO}{module}        .= "|";
            $hash{$KO}{module_desc}   .= "|";
        }else{
            if($hash{$KO}{mlevel1_desc} eq $mlevel1_desc || $hash{$KO}{mlevel1_desc} eq ""){
                $hash{$KO}{mlevel1_desc}   = $mlevel1_desc;
            }else{
                $hash{$KO}{mlevel1_desc}  .= "|".$mlevel1_desc;
            }
            
            if($hash{$KO}{mlevel2_desc} eq $mlevel2_desc || $hash{$KO}{mlevel2_desc} eq ""){
                $hash{$KO}{mlevel2_desc}   = $mlevel2_desc;
            }else{
                $hash{$KO}{mlevel2_desc}  .= "|".$mlevel2_desc;
            }

            if($hash{$KO}{mlevel3_desc} eq $mlevel3_desc || $hash{$KO}{mlevel3_desc} eq ""){
                $hash{$KO}{mlevel3_desc}   = $mlevel3_desc;
            }else{
                $hash{$KO}{mlevel3_desc}  .= "|".$mlevel3_desc;
            }
            
            if($hash{$KO}{module} eq $module || $hash{$KO}{module} eq ""){
                $hash{$KO}{module}         = $module;
            }else{
                $hash{$KO}{module}        .= "|".$module;
            }
            
            if($hash{$KO}{module_desc} eq $module_desc || $hash{$KO}{module_desc} eq ""){
                $hash{$KO}{module_desc}    = $module_desc;
            }else{
                $hash{$KO}{module_desc}   .= "|".$module_desc;
            }
        }
    }else{
        $hash{$KO}{KO}            = $KO;
        $hash{$KO}{KO_desc}       = $KO_desc;
        if($plevel3 eq "undefined"){ 
            $hash{$KO}{plevel1} = "";
            $hash{$KO}{plevel1_desc} = "";
            $hash{$KO}{plevel2} = "";
            $hash{$KO}{plevel2_desc} = "";
            $hash{$KO}{plevel3} = "";
            $hash{$KO}{plevel3_desc} = "";
        }else{
            $hash{$KO}{plevel1}       = $plevel1;
            $hash{$KO}{plevel1_desc}  = $plevel1_desc;
            $hash{$KO}{plevel2}       = $plevel2;
            $hash{$KO}{plevel2_desc}  = $plevel2_desc;
            $hash{$KO}{plevel3}       = $plevel3;
            $hash{$KO}{plevel3_desc}  = $plevel3_desc;
        }

        if($module eq "undefined"){
            $hash{$KO}{mlevel1_desc}  = ""; 
            $hash{$KO}{mlevel2_desc}  = "";
            $hash{$KO}{mlevel3_desc}  = "";
            $hash{$KO}{module}        = "";
            $hash{$KO}{module_desc}   = "";
        }else{
            $hash{$KO}{mlevel1_desc}  = $mlevel1_desc;
            $hash{$KO}{mlevel2_desc}  = $mlevel2_desc;
            $hash{$KO}{mlevel3_desc}  = $mlevel3_desc;
            $hash{$KO}{module}        = $module;
            $hash{$KO}{module_desc}   = $module_desc;
        }
    }
}
close(IN);

print STDERR Dumper(\%hash);

print STDOUT "KO\tKO_desc\tplevel1\tplevel1_desc\tplevel2\tplevel2_desc\tpathway\tpathway_desc\tmlevel1_desc\tmlevel2_desc\tmlevel3_desc\tmodule\tmodule_desc\n";
foreach my $KO (keys (%hash)){
    print $hash{$KO}{KO}."\t";
    print $hash{$KO}{KO_desc}."\t";

    my $plevel1 = $hash{$KO}{plevel1};
    $plevel1 =~ s/^\|.*//;
    $plevel1 =~ s/\|\|//;
    print STDOUT $plevel1."\t";

    #print $hash{$KO}{plevel1_desc}."\t";
    my $plevel1_desc = $hash{$KO}{plevel1_desc};
    $plevel1_desc =~ s/^\|.*//;
    $plevel1_desc =~ s/\|\|//;
    print STDOUT $plevel1_desc."\t";
    
    my $plevel2 = $hash{$KO}{plevel2};
    $plevel2 =~ s/^\|.*//;
    $plevel2 =~ s/\|\|//;
    print STDOUT $plevel2."\t";
    
    my $plevel2_desc = $hash{$KO}{plevel2_desc};
    $plevel2_desc =~ s/^\|.*//;
    $plevel2_desc =~ s/\|\|//;
    print STDOUT $plevel2_desc."\t";
    
    my $plevel3 = $hash{$KO}{plevel3};
    $plevel3 =~ s/^\|.*//;
    $plevel3 =~ s/\|\|//;
    print STDOUT $plevel3."\t";
    
    my $plevel3_desc = $hash{$KO}{plevel3_desc};
    $plevel3_desc =~ s/^\|.*//;
    $plevel3_desc =~ s/\|\|//;
    print STDOUT $plevel3_desc."\t";
    
    my $mlevel1_desc = $hash{$KO}{mlevel1_desc};
    $mlevel1_desc =~ s/^\|.*//;
    $plevel1_desc =~ s/\|\|//;
    print STDOUT $mlevel1_desc."\t";
    
    my $mlevel2_desc = $hash{$KO}{mlevel2_desc};
    $mlevel2_desc =~ s/^\|.*//;
    $plevel2_desc =~ s/\|\|//;
    print STDOUT $mlevel2_desc."\t";
    
    my $mlevel3_desc = $hash{$KO}{mlevel3_desc};
    $mlevel3_desc =~ s/^\|.*//;
    $plevel3_desc =~ s/\|\|//;
    print STDOUT $mlevel3_desc."\t";
    
    my $module = $hash{$KO}{module};
    $module =~ s/^\|.*//;
    $module =~ s/\|\|//;
    print STDOUT $module."\t";
    
    my $module_desc = $hash{$KO}{module_desc};
    $module_desc =~ s/^\|.*//;
    $module_desc =~ s/\|\|//;
    print STDOUT $module_desc."\n";
}
exit;


