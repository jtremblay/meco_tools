#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
convertTsneToJson.pl

PURPOSE:

INPUT:
--infile_tsne <string>     : tsne file <contig_id\txcoord\tycoord>.
--infile_length <string>   : length of contigs (txt \n file).
--infile_taxonomy <string> : infile taxonomy 
--infile_link <string>     : infile link. In file: binid\tcontig_id\tgene_id\n
--genus_select <string>    : Default:none. ex: Thalassolituus,Oleispira
--outfile_tsv <string>     : outfile tsv

OUTPUT:
STDOUT                     : same as gc table but in json format.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $infile_length, $infile_taxonomy, $infile_link, $genus_select, $outfile_tsv);
my $verbose = 0;

GetOptions(
   'infile_tsne=s'     => \$infile,
   'infile_length=s'   => \$infile_length,
   'infile_taxonomy=s' => \$infile_taxonomy,
   'infile_link=s'     => \$infile_link,
   'outfile_tsv=s'     => \$outfile_tsv,
   'genus_select=s'    => \$genus_select,
#   'infile_kegg=s'     => \$infile_kegg,
   'verbose'           => \$verbose,
   'help'              => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash_select;
if($genus_select){
    my @select = split(/,/, $genus_select);
    for my $el (@select){
        $hash_select{$el} = $el;
    }
}else{
    my @select = ();
}

my %hash;
my $total_el = 0;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){ 
   chomp;
   if($. == 1){
       next;
   }
   my @row = split(/\t/, $_);
   $hash{$row[0]}{x} = $row[1];
   $hash{$row[0]}{y} = $row[2];
   $total_el++;
   
}
close(IN);

open(IN, "<".$infile_length) or die "Can't open $infile_length\n";
while(<IN>){ 
   chomp;
   my @row = split(/\t/, $_);
   if(exists $hash{$row[0]}){
     $hash{$row[0]}{contig_length} = $row[1];
   }
}
close(IN);

open(IN, "<".$infile_link) or die "Can't open $infile_link\n";
while(<IN>){ 
   chomp;
   my @row = split(/\t/, $_);
   if(exists $hash{$row[1]}){
     $hash{$row[1]}{bin_id} = $row[0];
   }
}
close(IN);

open(IN, "<".$infile_taxonomy) or die "Can't open $infile_taxonomy\n";
while(<IN>){ 
   chomp;
   my @row = split(/\t/, $_);
    
   if($genus_select){
      my $found = 0;

      if(exists $hash_select{$row[7]}){
         if(exists $hash{$row[0]}){
            $hash{$row[0]}{kingdom} = $row[2];
            $hash{$row[0]}{phylum}  = $row[3];
            $hash{$row[0]}{class}   = $row[4];
            $hash{$row[0]}{order}   = $row[5];
            $hash{$row[0]}{family}  = $row[6];
            $hash{$row[0]}{genus}   = $row[7];
            $hash{$row[0]}{species} = $row[8];
         }
      }else{
        delete $hash{$row[0]};
      }
   }else{
       if(exists $hash{$row[0]}){
         $hash{$row[0]}{kingdom} = $row[2];
         $hash{$row[0]}{phylum}  = $row[3];
         $hash{$row[0]}{class}   = $row[4];
         $hash{$row[0]}{order}   = $row[5];
         $hash{$row[0]}{family}  = $row[6];
         $hash{$row[0]}{genus}   = $row[7];
         $hash{$row[0]}{species} = $row[8];
       }
   }
}
close(IN);

#open(IN, "<".$infile_kegg) or die "Can't open $infile_kegg\n";
#while(<IN>){ 
#   chomp;
#   my @row = split(/\t/, $_);
#   if(exists $hash{$row[0]}){
#     $hash{$row[0]}{KO} = $row[2];
#   }
#}
#close(IN);
if($outfile_tsv){
    open(OUT, ">".$outfile_tsv) or die "Can't open $outfile_tsv\n";
    print OUT "contig_id\ttax_kingdom\ttax_phylum\ttax_class\ttax_order\ttax_family\ttax_genus\ttax_species\tposx\tposy\tbin_id\tcontig_length\n";
}



print STDOUT "{\"nodes\":[{\n";
my $i = 1;
for my $key (keys %hash) {
    if(!defined($hash{$key})){ $i++;next;};
    if($hash{$key}{genus} eq "Zea"){$i++;next;}
	if($i == 1){
		print STDOUT "    \"contig_id\": \"".$key."\",\n";
		print STDOUT "    \"tax_kingdom\": \"".$hash{$key}{kingdom}."\",\n";
		print STDOUT "    \"tax_phylum\": \"".$hash{$key}{phylum}."\",\n";
		print STDOUT "    \"tax_class\": \"".$hash{$key}{class}."\",\n";
		print STDOUT "    \"tax_order\": \"".$hash{$key}{order}."\",\n";
		print STDOUT "    \"tax_family\": \"".$hash{$key}{family}."\",\n";
		print STDOUT "    \"tax_genus\": \"".$hash{$key}{genus}."\",\n";
		print STDOUT "    \"tax_species\": \"".$hash{$key}{species}."\",\n";
		print STDOUT "    \"KO\": [\"Tyrosine metabolism\",\"Benzoate degradation\"],\n";
		print STDOUT "    \"posx\": \"".$hash{$key}{x}."\",\n";
		print STDOUT "    \"posy\": \"".$hash{$key}{y}."\",\n";
		print STDOUT "    \"bin_id\": \"".$hash{$key}{bin_id}."\",\n";
		print STDOUT "    \"contig_length\": \"".$hash{$key}{bin_id}."\"\n";
	}elsif($i < $total_el){
    	print STDOUT "  }, {\n";
		print STDOUT "    \"contig_id\": \"".$key."\",\n";
		print STDOUT "    \"tax_kingdom\": \"".$hash{$key}{kingdom}."\",\n";
		print STDOUT "    \"tax_phylum\": \"".$hash{$key}{phylum}."\",\n";
		print STDOUT "    \"tax_class\": \"".$hash{$key}{class}."\",\n";
		print STDOUT "    \"tax_order\": \"".$hash{$key}{order}."\",\n";
		print STDOUT "    \"tax_family\": \"".$hash{$key}{family}."\",\n";
		print STDOUT "    \"tax_genus\": \"".$hash{$key}{genus}."\",\n";
		print STDOUT "    \"tax_species\": \"".$hash{$key}{species}."\",\n";
		print STDOUT "    \"KO\": [\"Tyrosine metabolism\",\"Benzoate degradation\"],\n";
		print STDOUT "    \"posx\": \"".$hash{$key}{x}."\",\n";
		print STDOUT "    \"posy\": \"".$hash{$key}{y}."\",\n";
		print STDOUT "    \"bin_id\": \"".$hash{$key}{bin_id}."\",\n";
		print STDOUT "    \"contig_length\": \"".$hash{$key}{contig_length}."\"\n";
	}elsif($i == $total_el){
    	print STDOUT "  }, {\n";
		print STDOUT "    \"contig_id\": \"".$key."\",\n";
		print STDOUT "    \"tax_kingdom\": \"".$hash{$key}{kingdom}."\",\n";
		print STDOUT "    \"tax_phylum\": \"".$hash{$key}{phylum}."\",\n";
		print STDOUT "    \"tax_class\": \"".$hash{$key}{class}."\",\n";
		print STDOUT "    \"tax_order\": \"".$hash{$key}{order}."\",\n";
		print STDOUT "    \"tax_family\": \"".$hash{$key}{family}."\",\n";
		print STDOUT "    \"tax_genus\": \"".$hash{$key}{genus}."\",\n";
		print STDOUT "    \"tax_species\": \"".$hash{$key}{species}."\",\n";
		print STDOUT "    \"KO\": [\"Tyrosine metabolism\",\"Benzoate degradation\"],\n";
		print STDOUT "    \"posx\": \"".$hash{$key}{x}."\",\n";
		print STDOUT "    \"posy\": \"".$hash{$key}{y}."\",\n";
		print STDOUT "    \"bin_id\": \"".$hash{$key}{bin_id}."\",\n";
		print STDOUT "    \"contig_length\": \"".$hash{$key}{contig_length}."\"\n";
		print STDOUT "  }\n";
	}
	$i++;

    if($outfile_tsv){
        print OUT $key."\t";
        print OUT $hash{$key}{kingdom}."\t";
        print OUT $hash{$key}{phylum}."\t";
        print OUT $hash{$key}{class}."\t";
        print OUT $hash{$key}{order}."\t";
        print OUT $hash{$key}{family}."\t";
        print OUT $hash{$key}{genus}."\t";
        print OUT $hash{$key}{species}."\t";
        print OUT $hash{$key}{x}."\t";
        print OUT $hash{$key}{y}."\t";
        print OUT $hash{$key}{bin_id}."\t";
        print OUT $hash{$key}{contig_length}."\n";
    }
}
print STDOUT "]}";
close(OUT);
exit;

#    var data = [{
#      "id": 1,
#      "first_name": "Beverlie",
#      "last_name": "Yeoman",
#      "email": "byeoman0@google.it",
#      "gender": "Female"
#    }, {
#      "id": 2,
#      "first_name": "Mabel",
#      "last_name": "Saurin",
#      "email": "msaurin1@facebook.com",
#      "gender": "Female"
#    }, {
#      "id": 3,
#      "first_name": "Magnum",
#      "last_name": "Moniker",
#      "email": "mmoniker2@geocities.com",
#      "gender": "Male"
#    }]
