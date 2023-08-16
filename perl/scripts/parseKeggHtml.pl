#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use JSON::PP;# qw( decode_json );
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseKeggHtml.pl

PURPOSE:
1- Get all the KEGG module files :
https://www.genome.jp/module/M00001
...
https://www.genome.jp/module/M00950
and run the getKeggModulesKOLink.pl script to generate one single file linking each KO with their corresponding module id.

2- Then get the 2 following web pages from the KEGG web page :
https://www.genome.jp/kegg/pathway.html  :  hierarchy file.
https://www.genome.jp/brite/ko00002.json :  Module file (click on the download json link).


INPUT:
--infile_pathways <string>        : The pathway.html file obtained with wget.
--infile_modules <string>         : KEGG module file
--infile_link_KO_modules <string> : file resulting from 1) above.
				
OUTPUT:
<STDOUT> <string>                 : Conveniently tsv format table of each KO with their corresponding pathways and modules.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_pathways, $infile_link_KO_modules, $infile_modules);
my $verbose = 0;

GetOptions(
   'infile_pathways=s'        => \$infile_pathways,
   'infile_link_KO_modules=s' => \$infile_link_KO_modules,
   'infile_modules=s'         => \$infile_modules,
   'verbose'                  => \$verbose,
   'help'                     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
# link between KO and modules
my %hash_KO_modules;
my %hash_KO_modules_seen;
open(IN, "<".$infile_link_KO_modules) or die "Can't open $infile_link_KO_modules\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    #row[1] = KO, row[0] = module ID
    if(exists $hash_KO_modules{$row[1]}){
        $hash_KO_modules{$row[1]} .= ",".$row[0];
    }else{
        $hash_KO_modules{$row[1]} = $row[0];
        $hash_KO_modules_seen{$row[1]}{"seen"} = "undef";
    }
}
close(IN);
#print STDERR Dumper(\%hash_KO_modules);

# Modules
open(IN, "<".$infile_modules) or die "Can't open $infile_modules\n";
my $data = join '', <IN>;
my $hash = decode_json($data);
close(IN);

#%hash = undef(%hash);
my %hash = %$hash;
#print STDERR Dumper(\%hash);
my %hash_modules;
foreach my $level1 (@{ $hash{children} }) {
        
    my $str1 = %$level1{name}."\n"; # actual string value;
    chomp($str1);
    $str1 =~ s/^\s+//;
    #print STDOUT "L1: ".$str1."\n";
    my %hash2 = %$level1;

    ##
    foreach my $level2 (@{ $hash2{children} }) {
        my $str2 = %$level2{name}."\n"; # actual string value;
        chomp($str2);
        $str2 =~ s/^\s+//;
        #print STDOUT "L2: ".$str2."\n";
        my %hash3 = %$level2;

        ##
        foreach my $level3 (@{ $hash3{children} }) {
            my $str3 = %$level3{name}."\n"; # actual string value;
            chomp($str3);
            $str3 =~ s/^\s+//;
            #print STDOUT "L3: ".$str3."\n";
            my %hash4 = %$level3;
        
            foreach my $level4 (@{ $hash4{children} }) {
                my $str4 = %$level4{name}."\n"; # actual string value;
                chomp($str4);
                #print STDOUT "L4: ".$str4."\n";
                # split module id and module desc.
                my $module = "";
                my $module_desc = "";
                if($str4 =~ m/^(M\d{5}) (.*) \[.*\]/){
                    $module = $1;
                    $module_desc = $2;
                    $module_desc =~ s/^\s+//;
                }elsif($str4 =~ m/^(M\d{5}) (.*)$/){
                    $module = $1;
                    $module_desc = $2;
                    $module_desc =~ s/^\s+//;
                }
            
                $hash_modules{$module}{mlevel1_desc} = $str1;
                $hash_modules{$module}{mlevel2_desc} = $str2;
                $hash_modules{$module}{mlevel3_desc} = $str3;
                $hash_modules{$module}{module} = $module;
                $hash_modules{$module}{module_desc} = $module_desc;
            }
        }
    }
}
#print STDERR Dumper(%hash_modules);
# finally the pathways.
my %final_hash;
$hash=undef($hash);
open(IN, "<".$infile_pathways) or die "Can't open $infile_pathways\n";
while(<IN>){
    chomp;
    #my @row = split(/\t/, $_);
    if($_ =~ m/root\:/){
        #print STDERR "found root:\n";
        my $ref = $_;
        $ref =~ s/^\s+root\://;
        #$ref =~ s/<a href=\\"\/entry//g;
        #$ref =~ s/\\">//g;
        #$ref =~ s/<a href=\\"\/entry\/.*\\">//g;
        $hash = decode_json $ref;
    }
}
close(IN);

%hash = %$hash;

# Hash Hash of array of hashes of array of hashes... :D :D

# Pathways
foreach my $level1 (@{ $hash{children} }) {
    if(!defined( @{ %$level1{values} }[0] )){next;}
        
    my $str1 = @{ %$level1{values} }[0]."\n"; # actual string value;
    #print STDOUT "L1: ".$str1."\n";
    my %hash2 = %$level1;
    my $plevel1_id = "";
    my $plevel1_desc = "";
    if($str1 =~ m/^(\d+) (.*)/){
        $plevel1_id = $1;
        $plevel1_desc = $2;
    }


    ##
    foreach my $level3 (@{ $hash2{children} }) {
        if(!defined( @{ %$level3{values} }[0] )){next;}
        my $str2 = @{ %$level3{values} }[0]."\n"; # actual string value;
        #print STDOUT "L2: ".$str2."\n";
        my %hash3 = %$level3;
        my $plevel2_id = "";
        my $plevel2_desc = "";
        if($str2 =~ m/^(\d+) (.*)/){
            $plevel2_id = $1;
            $plevel2_desc = $2;
        }

        ##
        foreach my $level5 (@{ $hash3{children} }) {
            if(!defined( @{ %$level5{values} }[0] )){next;}
            my $str3 = @{ %$level5{values} }[0]."\n"; # actual string value;
            my %hash4 = %$level5;
            #print STDOUT "L3: ".$str3."\n";
            my $plevel3_id = "";
            my $plevel3_desc = "";
            if($str3 =~ m/^(\d+) (.*) \[.*\]/){
                $plevel3_id = "ko".$1;
                $plevel3_desc = $2;
            }
            
            ##
            foreach my $level7 (@{ $hash4{children} }) {
                if(!defined( @{ %$level7{values} }[0] )){next;}
                my $str4 = @{ %$level7{values} }[0]."\n"; # actual string value;
                my %hash5 = %$level7;
                #print STDOUT "L4: ".$str4."\n";
                my $KO = "";
                my $KO_desc = $str4;
                if($str4 =~ m/(K\d{5})/){
                    $KO = $1;
                }
                $KO_desc =~ s/\<a.*\/a\> //;
                $KO_desc =~ s/\[.*\]//;
                chomp($KO_desc);

                $final_hash{$KO."-".$plevel3_id}{plevel1_id} = $plevel1_id;
                $final_hash{$KO."-".$plevel3_id}{plevel1_desc} = $plevel1_desc;
                $final_hash{$KO."-".$plevel3_id}{plevel2_id} = $plevel2_id;
                $final_hash{$KO."-".$plevel3_id}{plevel2_desc} = $plevel2_desc;
                $final_hash{$KO."-".$plevel3_id}{plevel3_id} = $plevel3_id;
                $final_hash{$KO."-".$plevel3_id}{plevel3_desc} = $plevel3_desc;
                $final_hash{$KO."-".$plevel3_id}{KO_desc} = $KO_desc;
                $final_hash{$KO."-".$plevel3_id}{KO} = $KO;
                
                # Then add modules. Valid for modules that may or may not be present in pathways.
                # But what about modules that exist, but have no associated pathways? They won't be catched here. Need another loop after this one...
                if(exists $hash_KO_modules{$KO}){ # can be multiple modules.
                    $hash_KO_modules_seen{$KO}{seen} = "yes";
                    my $modules = $hash_KO_modules{$KO};
                    my @modules = split(/,/, $modules);
                    foreach my $module (@modules){
                        $final_hash{$KO."-".$plevel3_id}{module} = $module;
                        $final_hash{$KO."-".$plevel3_id}{module_desc} = $hash_modules{$module}{module_desc};
                        $final_hash{$KO."-".$plevel3_id}{mlevel1_desc} = $hash_modules{$module}{mlevel1_desc};
                        $final_hash{$KO."-".$plevel3_id}{mlevel2_desc} = $hash_modules{$module}{mlevel2_desc};
                        $final_hash{$KO."-".$plevel3_id}{mlevel3_desc} = $hash_modules{$module}{mlevel3_desc};
                    }
                }else{
                    $hash_KO_modules_seen{$KO}{seen} = "no";
                    $final_hash{$KO."-".$plevel3_id}{module} = "undefined"; 
                    $final_hash{$KO."-".$plevel3_id}{module_desc} = "undefined";
                    $final_hash{$KO."-".$plevel3_id}{mlevel1_desc} = "undefined";
                    $final_hash{$KO."-".$plevel3_id}{mlevel2_desc} = "undefined";
                    $final_hash{$KO."-".$plevel3_id}{mlevel3_desc} = "undefined";
                }
            }
        }
    }
}


#print STDERR Dumper(\%final_hash);
# Then print final hash;
print STDOUT "plevel1\tplevel1_desc\tplevel2\tplevel2_desc\tpathway\tpathway_desc\tKO\tKO_desc\tmlevel1_desc\tmlevel2_desc\tmlevel3_desc\tmodule\tmodule_desc\n";
foreach my $KO_plevel3 (keys %final_hash) {
    print STDOUT $final_hash{$KO_plevel3}{plevel1_id}."\t";
    print STDOUT $final_hash{$KO_plevel3}{plevel1_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{plevel2_id}."\t";
    print STDOUT $final_hash{$KO_plevel3}{plevel2_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{plevel3_id}."\t";
    print STDOUT $final_hash{$KO_plevel3}{plevel3_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{KO}."\t";
    print STDOUT $final_hash{$KO_plevel3}{KO_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{mlevel1_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{mlevel2_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{mlevel3_desc}."\t";
    print STDOUT $final_hash{$KO_plevel3}{module}."\t";
    print STDOUT $final_hash{$KO_plevel3}{module_desc}."\n";
    #if(!defined $final_hash{$KO_plevel3}{mlevel1_desc}){print STDERR "undef:".$KO_plevel3."\n";}
}

#Finally, loop back into the module hash to find modules with not pathways.
foreach my $KO (keys %hash_KO_modules_seen) {
    if($hash_KO_modules_seen{$KO}{seen} eq "undef"){
        print STDOUT "undefined\t"; 
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT $KO."\t";
        print STDOUT "undefined\t"; 
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT "undefined\t";
        print STDOUT "undefined\n";
    }
}
#print STDERR Dumper(\%hash_KO_modules_seen);
exit;
