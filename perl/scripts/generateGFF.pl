#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
generateGFF.pl

PURPOSE:

INPUT:
--infile_fasta <string>     : Sequence file
--infile_gff <string>       : gff file (i.e. output of prodigal)
--prefix <string>           : will add a prefix to contig names (i.e. meaningful project id)
--pfam <string              : pfam hmmscan results (.tblout)
--cog <string               : cog rpsblast results
--kog <string               : kog rpsblasts results
--kegg <string              : kegg blastpn results
--taxonomy <string>         : taxonomy results from both NCBI nt blastn and NCBI nr - diamond blastp nr/refseq
--ublast <string>           : ublast_nr results - i.e. can submit an empty file.
--taxonomy_type <string>    : 'consensus' or 'contigs'. Default='contigs'

OUTPUT:
--outfile_gff <string>      : Fully annotated gff. To use in genome viewers.
--outfile_fasta <string>    : Fasta file with compatible headers for genome viewers.
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
, Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_fasta, $infile_gff, $pfam, $cog, $kog, $kegg, $taxonomy, $ublast, $prefix, $outfile_fasta, $outfile_gff, $taxonomy_type);
my $verbose = 0;

GetOptions(
   'infile_fasta=s'     => \$infile_fasta,
   'infile_gff=s'       => \$infile_gff,
   'outfile_fasta=s'    => \$outfile_fasta,
   'outfile_gff=s'      => \$outfile_gff,
   'prefix=s'           => \$prefix,
   'pfam=s'             => \$pfam,
   'cog=s'              => \$cog,
   'kog=s'              => \$kog,
   'kegg=s'             => \$kegg,
   'ublast=s'           => \$ublast,
   'taxonomy=s'         => \$taxonomy,
   'taxonomy_type=s'    => \$taxonomy_type,
   'verbose' 	        => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile_fasta is missing\n" unless($infile_fasta);
die "--infile_gff is missing\n" unless($infile_gff);
die "--outfile_fasta is missing\n" unless($outfile_fasta);
die "--outfile_gff is missing\n" unless($outfile_gff);
die "--prefix is missing\n" unless($prefix);
die "--pfam is missing\n" unless($pfam);
die "--cog is missing\n" unless($cog);
die "--kog is missing\n" unless($kog);
die "--kegg is missing\n" unless($kegg);
die "--taxonomy is missing\n" unless($taxonomy);
#die "--taxonomy_type is missing\n" unless($taxonomy_type);
die "--ublast is missing\n" unless($ublast);

#die if($taxonomy_type ne "consensus" && $taxonomy_type ne "contigs");
$taxonomy_type = "contigs" unless($taxonomy_type);


## MAIN
# First, rename contigs (remove spurious characters in headers). For instance convert this: ">contig-0 34353 nucleotides" to this: ">contig-0"
my %hash;
my %hash_gene_to_contig;

open(OUT_FASTA, ">".$outfile_fasta) or die "Can't open $outfile_fasta\n";
open(OUT_GFF, ">".$outfile_gff) or die "Can't open $outfile_gff\n";

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ){
   my $header = $curr->header;
   $header = $1 if $header =~ m/(>\S+).*$/;
   print OUT_FASTA $header."\n".$curr->seq."\n";
}
close(OUT_FASTA);

my $megahit_flag = 0;
my $megahit_id = "";

# Parse metagenemark's gff output.
open(IN, "<".$infile_gff) or die "Can't open $infile_gff\n";
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   my $contig_id = $row[0];

   # Hack if assembly is done with ray or megahit
   if($contig_id =~ m/^contig-(\d+)/){
      $contig_id = $1; 
   }elsif($contig_id =~ m/^(k\d+_)(\d+) /){
      $contig_id = $2;
      $megahit_flag = 1;
      $megahit_id = $1;
   }elsif($contig_id =~ m/^i(scaffold\d+\|size\d+)/){
      $contig_id = $1;
      $megahit_flag = 0;
      $megahit_id = $1;
   }

   #1 to 7
   my $source     = $row[1];
   my $type       = $row[2];
   my $start      = $row[3];
   my $end        = $row[4];
   my $score      = $row[5];
   my $strand     = $row[6];
   my $phase      = $row[7];
   my $attributes = $row[8];

   my @attributes = split(/;/, $attributes);
   my $gene_id = $attributes[0];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   #my @gene_id_row = split(/;/, $gene_id);
   #$gene_id = $gene_id_row[0];

   $hash{$contig_id}{$gene_id}{SOURCE} = $source;
   $hash{$contig_id}{$gene_id}{TYPE}   = $type;
   $hash{$contig_id}{$gene_id}{START}  = $start;
   $hash{$contig_id}{$gene_id}{END}    = $end;
   $hash{$contig_id}{$gene_id}{SCORE}  = $score;
   $hash{$contig_id}{$gene_id}{STRAND} = $strand;
   $hash{$contig_id}{$gene_id}{PHASE}  = $phase;
      
   $hash_gene_to_contig{$gene_id} = $contig_id;
   
   #if($. > 1000){
   #   print STDERR Dumper(\%hash);
   #   print STDERR Dumper(\%hash_gene_to_contig);
   #   exit;
   #}
}
close(IN);
#print STDERR Dumper(\%hash);
#print STDERR Dumper(\%hash_gene_to_contig);

open(IN, "<".$pfam) or die "Can't open $pfam\n";
my %seen;
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/.*#.*/); #new since using hmmsearch instead of hmmscan
   next if($_ =~ m/^$/);
   my @row = split(/\s+/, $_);
   my $gene_id = $row[2];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   if(!exists $seen{$gene_id}){
      # get contig id first...
      my $contig_id = $hash_gene_to_contig{$gene_id};
      next if($contig_id eq "");
      #print STDERR "ContigID: $contig_id\n";

      $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME} = $row[0];
      $hash{$contig_id}{$gene_id}{PFAM_ACCESS} = $row[1];

      # Then get desc.
      my @desc = @row[18 .. (@row - 1)];
      my $desc = join(" ", @desc);
      $hash{$contig_id}{$gene_id}{PFAM_DESC} = $desc;

      $seen{$gene_id} = 1;
   }
   
   #if($. > 100){
   #   print STDERR Dumper(\%hash);
   #   print STDERR Dumper(\%hash_gene_to_contig);
   #   last;
   #}
}
close(IN); 
#print STDERR Dumper(\%hash);
#print STDERR Dumper(\%hash_gene_to_contig);

## Tigrfam deprecated - just add dummy value in a dummy hash.
#open(IN, "<".$tigrfam) or die "Can't open $tigrfam\n";
#undef %seen;
#while(<IN>){
#   chomp;
#   next if($_ =~ m/^#/);
#   next if($_ =~ m/^$/);
#   my @row = split(/\s+/, $_);
#   my $gene_id = $row[2];
#   $gene_id =~ s/gene_id=//;
#   $gene_id =~ s/gene_id_//;
#   if(!exists $seen{$gene_id}){
#      # get contig id first...
#      my $contig_id = $hash_gene_to_contig{$gene_id};
#      
#      # Parse field $row[18]
#      my $length = scalar(@row) - 1;
#      my @slice =  @row[18..$length];
#      my $desc_field = join(" ", @slice);
#      my @desc_field = split(/:/, $desc_field);
#      my $product_name = $desc_field[0];
#
#      $hash{$contig_id}{$gene_id}{TIGRFAM_NAME} = $row[0];
#      $hash{$contig_id}{$gene_id}{TIGRFAM_ACCESS} = $row[1];
#      $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME} = $product_name;
#      
#      # Then get desc.
#      my @desc = @row[18 .. (@row - 1)];
#      my $desc = join(" ", @desc);
#      $hash{$contig_id}{$gene_id}{TIGRFAM_DESC} = $desc;
#
#      $seen{$gene_id} = 1;
#   }
#   
#   #if($. > 100){
#   #    print STDERR Dumper(\%hash);
#   #    print STDERR Dumper(\%hash_gene_to_contig);
#   #   last;
#   #}
#}
#close(IN);
#print STDERR Dumper(\%hash);
#print STDERR Dumper(\%hash_gene_to_contig);
#$hash{"contig_id_000"}{"gene_id_000"}{TIGRFAM_NAME} = "NULL";
#$hash{"contig_id_000"}{"gene_id_000"}{TIGRFAM_ACCESS} = "NULL";
#$hash{"contig_id_000"}{"gene_id_000"}{TIGRFAM_PRODUCT_NAME} = "NULL";

open(IN, "<".$cog) or die "Can't open $cog\n";
undef %seen;
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   #0, 12, 13
   my $gene_id = $row[0];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   if(!exists $seen{$gene_id}){
      # get contig id first...
      my $contig_id = $hash_gene_to_contig{$gene_id};
      
      my @terms = split(/,/, $row[12]);
      $hash{$contig_id}{$gene_id}{COG_ACCESS} = $terms[0];
      my $term1 = $terms[1];
      $term1 =~ s/^\s+//;
      $hash{$contig_id}{$gene_id}{COG_NAME} = $term1;
      my @desc = @terms[2 .. (@terms - 1)];
      my $desc = join(",", @desc);
      $desc =~ s/^\s+//;
      my $cog1;
      my $cog2;
      if($desc =~ m/^(.*)\[/){
         $cog1 = $1;
         $cog1 =~ s/\s+$//; 
         $hash{$contig_id}{$gene_id}{COG_FUNCTION} = $cog1;
      }
      if($desc =~ m/\[(.*)\]/){
         $cog2 = $1; 
         $cog2 =~ s/\s+$//; 
         $hash{$contig_id}{$gene_id}{COG_CATEGORY} = $cog2;
      }

      $seen{$gene_id} = 1;
   }
   
   #if($. > 100){
   #print STDERR Dumper(\%hash);
   #    print STDERR Dumper(\%hash_gene_to_contig);
   #   last;
   #}
}
close(IN);
#print STDERR Dumper(\%hash);
#print STDERR Dumper(\%hash_gene_to_contig);

open(IN, "<".$kog) or die "Can't open $kog\n";
undef %seen;
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   #0, 12, 13
   my $gene_id = $row[0];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   if(!exists $seen{$gene_id}){
      # get contig id first...
      my $contig_id = $hash_gene_to_contig{$gene_id};

      my @terms = split(/,/, $row[12]);
      $hash{$contig_id}{$gene_id}{KOG_ACCESS} = $terms[0];
      my $term1 = $terms[1];
      $term1 =~ s/^\s+//;
      $hash{$contig_id}{$gene_id}{KOG_NAME} = $term1;
      my @desc = @terms[3 .. (@terms - 1)];
      my $desc = join(",", @desc);
      $desc =~ s/^\s+//;
      my $cog1;
      my $cog2;
      if($desc =~ m/^(.*)\[/){
         $cog1 = $1;
         $cog1 =~ s/\s+$//; 
         $hash{$contig_id}{$gene_id}{KOG_FUNCTION} = $cog1;
      }
      if($desc =~ m/\[(.*)\]/){
         $cog2 = $1; 
         $cog2 =~ s/\s+$//; 
         $hash{$contig_id}{$gene_id}{KOG_CATEGORY} = $cog2;
      }

      $seen{$gene_id} = 1;
   }
   
   #if($. > 100){
   #   print STDERR Dumper(\%hash);
   #   print STDERR Dumper(\%hash_gene_to_contig);
   #   last;
   #}
}
close(IN);

# Ublastp nr
open(IN, "<".$ublast) or die "Can't open $ublast\n";
undef %seen;
while(<IN>){
   chomp;
   if($. == 1){
      next;
   }
   my @row = split(/\t/, $_);
   #0, 12, 13
   my $gene_id = $row[0];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   if(!exists $seen{$gene_id}){
      # get contig id first...
      my $contig_id = $hash_gene_to_contig{$gene_id};

      $hash{$contig_id}{$gene_id}{UBLASTP_NR} = $row[1];
      $hash{$contig_id}{$gene_id}{UBLASTP_NR_DESC} = $row[12];

      $seen{$gene_id} = 1;
   }
}
close(IN);

# Contigs-based taxonomy...
my %hash_tax;
if($taxonomy_type eq "contigs"){
    open(IN, "<".$taxonomy) or die "Can't open $taxonomy\n";
    undef %seen;
    while(<IN>){
       chomp;
       if($. == 1){
          next;
       }
       my @row = split(/\t/, $_);
       #0, 12, 13
       my $contig_id = $row[0];
       if($megahit_flag == 1){ # Megahit
          if($contig_id =~ m/^(k\d+_)(\d+)/){
             $contig_id = $2;
          }
       }else{ #Ray
          if($contig_id =~ m/^contig-(\d+)/){
             $contig_id = $1;
          }elsif($contig_id =~ m/^(scaffold\d+\|size\d+)/){
             $contig_id = $1;
          }
       }
       
       if(!exists $seen{$contig_id}){
          # get contig id first...
          #my $contig_id = $hash_gene_to_contig{$gene_id};
    
          $hash_tax{$contig_id}{TAX_K} = $row[1];
          $hash_tax{$contig_id}{TAX_P} = $row[2];
          $hash_tax{$contig_id}{TAX_C} = $row[3];
          $hash_tax{$contig_id}{TAX_O} = $row[4];
          $hash_tax{$contig_id}{TAX_F} = $row[5];
          $hash_tax{$contig_id}{TAX_G} = $row[6];
          $hash_tax{$contig_id}{TAX_S} = $row[7];
    
          $seen{$contig_id} = 1;
       }
    }
    close(IN);
}else{
    #print STDERR Dumper(\%hash_tax);
    #print STDERR Dumper(\%hash_gene_to_contig);
    #exit;

    # Gene based taxonomy - consensus taxonomy from blastn NCBI nt and diamond blastp NCBI nr
    open(IN, "<".$taxonomy) or die "Can't open $taxonomy\n";
    undef %seen;
    while(<IN>){
       chomp;
       if($. == 1){
          next;
       }
       my @row = split(/\t/, $_);
       #0, 12, 13
       my $gene_id = $row[0];
       $gene_id =~ s/gene_id=//;
       $gene_id =~ s/gene_id_//;
       if(!exists $seen{$gene_id}){
          # get contig id first...
          my $contig_id = $hash_gene_to_contig{$gene_id};
    
          $hash_tax{$contig_id}{$gene_id}{TAX_K} = $row[2];
          $hash_tax{$contig_id}{$gene_id}{TAX_P} = $row[3];
          $hash_tax{$contig_id}{$gene_id}{TAX_C} = $row[4];
          $hash_tax{$contig_id}{$gene_id}{TAX_O} = $row[5];
          $hash_tax{$contig_id}{$gene_id}{TAX_F} = $row[6];
          $hash_tax{$contig_id}{$gene_id}{TAX_G} = $row[7];
          $hash_tax{$contig_id}{$gene_id}{TAX_S} = $row[8];
    
          $seen{$gene_id} = 1;
       }
    }
    close(IN);
    print STDERR Dumper(\%hash_tax);
    #print STDERR Dumper(\%hash_gene_to_contig);
}

open(IN, "<".$kegg) or die "Can't open $kegg\n";
undef %seen;
while(<IN>){
   chomp;
   next if($_ =~ m/^#/);
   next if($_ =~ m/.*#.*/);
   next if($_ =~ m/^Program/);
   next if($_ =~ m/^$/);
   my @row = split(/\t/, $_);
   #5-7
   my $gene_id = $row[0];
   $gene_id =~ s/gene_id=//;
   $gene_id =~ s/gene_id_//;
   if(!exists $seen{$gene_id}){
      # get contig id first...
      my $contig_id = $hash_gene_to_contig{$gene_id};
   
      # OLD ORDER DO NOT USE.
      #0-query   
      #1-kegg_gene_id  
      #2-KO_id 
      #3-NAME  
      #4-ENTRY 
      #5-DEFINITION  
      #6-PATHWAY_ID  
      #7-PATHWAY_DESC   
      #8-MODULE_ID   
      #9-MODULE_DESC

      $hash{$contig_id}{$gene_id}{KEGG_ENTRY} = $row[2];
      $hash{$contig_id}{$gene_id}{KEGG_NAME} = $row[3];
      my $kegg_desc = $row[5]; #definition.
      $kegg_desc =~ s/\#|;|\&|\"//g;
      $hash{$contig_id}{$gene_id}{KEGG_DEFINITION} = $kegg_desc;
      $hash{$contig_id}{$gene_id}{KEGG_PATHWAY} = $row[6];
      $hash{$contig_id}{$gene_id}{KEGG_PATHWAY_DESC} = $row[7];
      $hash{$contig_id}{$gene_id}{KEGG_MODULE} = $row[8];
      $hash{$contig_id}{$gene_id}{KEGG_MODULE_DESC} = $row[9];

      $seen{$gene_id} = 1;
   }
   #if($. > 100){
   #   print STDERR Dumper(\%hash);
   #   print STDERR Dumper(\%hash_gene_to_contig);
   #   last;
   #}
}
close(IN);
#print STDERR "dumper hash\n";
#print STDERR Dumper(\%hash);

# Once everything is annotated, print final GFF
for my $contig_id (sort {$a cmp $b} keys %hash) {
   for my $gene_id (sort {$a cmp $b} keys %{ $hash{$contig_id} }){

      # Chose which name to display follow JGI guidelines.
      # tigrgam -> COG -> Pfam -> Blast -> hypothetical
      my $product_name;
      if( (exists $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME}) && ($hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME} !~ m/TIGR\d+/) ){
         $product_name = $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME};

      }elsif(exists $hash{$contig_id}{$gene_id}{COG_NAME}){
         $product_name = $hash{$contig_id}{$gene_id}{COG_NAME};

      }elsif(exists $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME}){
          $product_name = $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME};

      #}elsif(exists $hash{$contig_id}{$gene_id}{BLASTN_PRODUCT_NAME}){ 

      }else{
         $product_name = "hypothetical";
      }

      #if($megahit_flag == 0){ print OUT_GFF "contig-".$contig_id."\t"; } else { print OUT_GFF $megahit_id.$contig_id."\t"; }
      if($megahit_flag == 0){ print OUT_GFF $contig_id."\t"; } else { print OUT_GFF $megahit_id.$contig_id."\t"; }
      print OUT_GFF $hash{$contig_id}{$gene_id}{SOURCE}."\t";
      print OUT_GFF $hash{$contig_id}{$gene_id}{TYPE}."\t";
      print OUT_GFF $hash{$contig_id}{$gene_id}{START}."\t";
      print OUT_GFF $hash{$contig_id}{$gene_id}{END}."\t";
      print OUT_GFF $hash{$contig_id}{$gene_id}{SCORE}."\t";
      print OUT_GFF $hash{$contig_id}{$gene_id}{STRAND}."\t";
      print OUT_GFF $hash{$contig_id}{$gene_id}{PHASE}."\t";
      print OUT_GFF "ID=gene_id_$gene_id; ";
      print OUT_GFF "gene_id=gene_id_$gene_id; ";
      print OUT_GFF "Name=$product_name; ";

      if(exists $hash{$contig_id}{$gene_id}{KEGG_ENTRY}){
         print OUT_GFF "Kegg_entry=".$hash{$contig_id}{$gene_id}{KEGG_ENTRY}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{KEGG_DEFINITION}){
         print OUT_GFF "Kegg_definition=".$hash{$contig_id}{$gene_id}{KEGG_DEFINITION}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{KEGG_MODULE}){
         print OUT_GFF "Kegg_module=".$hash{$contig_id}{$gene_id}{KEGG_MODULE}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{KEGG_PATHWAY}){
         print OUT_GFF "kegg_pathway=".$hash{$contig_id}{$gene_id}{KEGG_PATHWAY}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{PFAM_ACCESS}){
         print OUT_GFF "pfam_accession=".$hash{$contig_id}{$gene_id}{PFAM_ACCESS}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME}){
         print OUT_GFF "pfam_name=".$hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_ACCESS}){
         print OUT_GFF "tigrfam_accession=".$hash{$contig_id}{$gene_id}{TIGRFAM_ACCESS}."; "; 
      }
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_NAME}){
         print OUT_GFF "tigrfam_name=".$hash{$contig_id}{$gene_id}{TIGRFAM_NAME}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME}){
         print OUT_GFF "tigrfam_product_name=".$hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{COG_ACCESS}){
         print OUT_GFF "cog_accession=".$hash{$contig_id}{$gene_id}{COG_ACCESS}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{COG_NAME}){
         print OUT_GFF "cog_name=".$hash{$contig_id}{$gene_id}{COG_NAME}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{KOG_ACCESS}){
         print OUT_GFF "kog_accession=".$hash{$contig_id}{$gene_id}{KOG_ACCESS}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{KOG_NAME}){
         print OUT_GFF "kog_name=".$hash{$contig_id}{$gene_id}{KOG_NAME}."; ";
      }

      # If we ever need to revert to gene_id based taxonomy, use hash struture like this: $hash{$contig_id}{gene_id}{TAX_S}
      if(exists $hash_tax{$contig_id}{TAX_K}){
         print OUT_GFF "tax_kingdom=".$hash_tax{$contig_id}{TAX_K}."; ";
      }
      if(exists $hash_tax{$contig_id}{TAX_P}){
         print OUT_GFF "tax_phylum=".$hash_tax{$contig_id}{TAX_P}."; ";
      }
      if(exists $hash_tax{$contig_id}{TAX_C}){
         print OUT_GFF "tax_class=".$hash_tax{$contig_id}{TAX_C}."; ";
      }
      if(exists $hash_tax{$contig_id}{TAX_O}){
         print OUT_GFF "tax_order=".$hash_tax{$contig_id}{TAX_O}."; ";
      }
      if(exists $hash_tax{$contig_id}{TAX_F}){
         print OUT_GFF "tax_family=".$hash_tax{$contig_id}{TAX_F}."; ";
      }
      if(exists $hash_tax{$contig_id}{TAX_G}){
         print OUT_GFF "tax_genus=".$hash_tax{$contig_id}{TAX_G}."; ";
      }
      if(exists $hash_tax{$contig_id}{TAX_S}){
         print OUT_GFF "tax_species=".$hash_tax{$contig_id}{TAX_S}."; ";
      }

      if(exists $hash{$contig_id}{$gene_id}{UBLASTP_NR}){
         print OUT_GFF "ublastp_nr=".$hash{$contig_id}{$gene_id}{UBLASTP_NR}."; ";
      }
      if(exists $hash{$contig_id}{$gene_id}{UBLASTP_NR_DESC}){
         print OUT_GFF "ublastp_nr_desc=".$hash{$contig_id}{$gene_id}{UBLASTP_NR_DESC}."; ";
      }
      print OUT_GFF "note=$prefix; \n";
   }
}
close(OUT_GFF);

# Then print reference tsv table
print STDOUT "#contig_id\tgene_id\tproduct_name\tkegg_entry\tkegg_definition\tkegg_module\tkegg_module_desc\tkegg_pathway\tkegg_pathway_desc\t";
print STDOUT "pfam_access\tpfam_product_name\tpfam_desc\ttigrfam_access\ttigrfam_name\ttigrfam_product_name\ttigrfam_desc\t";
print STDOUT "cog_access\tcog_name\tcog_function\tcog_category\tkog_access\tkog_name\tkog_function\tkog_category\tublastp\tublastp_desc\ttax_kingdom\ttax_phylum\ttax_class\ttax_order\ttax_family\ttax_genus\ttax_specie\n";
for my $contig_id (sort {$a cmp $b} keys %hash) {
   for my $gene_id (sort {$a cmp $b} keys %{ $hash{$contig_id} }){

      # Chose which name to display follow JGI guidelines.
      # tigrgam -> COG -> Pfam -> Blast -> hypothetical
      # Actually... I found that the KEGG gene names where more meaningful. So it'll be: KEGG -> tigrgam -> COG -> Pfam -> Blast -> hypothetical instead...
      my $product_name;
      if( (exists $hash{$contig_id}{$gene_id}{KEGG_NAME}) && $hash{$contig_id}{$gene_id}{KEGG_NAME} !~ m/E\d+\.\d+\.\d+/){
         $product_name = $hash{$contig_id}{$gene_id}{KEGG_NAME};

      }elsif( (exists $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME}) && ($hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME} !~ m/TIGR\d+/) ){
         $product_name = $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME};

      }elsif(exists $hash{$contig_id}{$gene_id}{COG_NAME}){
         $product_name = $hash{$contig_id}{$gene_id}{COG_NAME};

      }elsif(exists $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME}){
          $product_name = $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME};

      #}elsif(exists $hash{$contig_id}{$gene_id}{BLASTN_PRODUCT_NAME}){ 

      }else{
         $product_name = "hypothetical";
      }

      if($megahit_flag == 0){ print STDOUT $contig_id."\t"; } else { print STDOUT $megahit_id.$contig_id."\t"; }
      #if($megahit_flag == 0){ print STDOUT "contig-".$contig_id."\t"; } else { print STDOUT $megahit_id.$contig_id."\t"; }
      print STDOUT "gene_id_$gene_id\t";
      print STDOUT "$product_name\t";
      if(exists $hash{$contig_id}{$gene_id}{KEGG_ENTRY}){           print STDOUT $hash{$contig_id}{$gene_id}{KEGG_ENTRY}."\t";}           else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KEGG_DEFINITION}){      print STDOUT $hash{$contig_id}{$gene_id}{KEGG_DEFINITION}."\t";}      else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KEGG_MODULE}){          print STDOUT $hash{$contig_id}{$gene_id}{KEGG_MODULE}."\t";}          else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KEGG_MODULE_DESC}){     print STDOUT $hash{$contig_id}{$gene_id}{KEGG_MODULE_DESC}."\t";}     else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KEGG_PATHWAY}){         print STDOUT $hash{$contig_id}{$gene_id}{KEGG_PATHWAY}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KEGG_PATHWAY_DESC}){    print STDOUT $hash{$contig_id}{$gene_id}{KEGG_PATHWAY_DESC}."\t";}    else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{PFAM_ACCESS}){          print STDOUT $hash{$contig_id}{$gene_id}{PFAM_ACCESS}."\t";}          else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME}){    print STDOUT $hash{$contig_id}{$gene_id}{PFAM_PRODUCT_NAME}."\t";}    else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{PFAM_DESC}){            print STDOUT $hash{$contig_id}{$gene_id}{PFAM_DESC}."\t";}            else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_ACCESS}){       print STDOUT $hash{$contig_id}{$gene_id}{TIGRFAM_ACCESS}."\t";}       else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_NAME}){         print STDOUT $hash{$contig_id}{$gene_id}{TIGRFAM_NAME}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME}){ print STDOUT $hash{$contig_id}{$gene_id}{TIGRFAM_PRODUCT_NAME}."\t";} else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{TIGRFAM_DESC}){         print STDOUT $hash{$contig_id}{$gene_id}{TIGRFAM_DESC}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{COG_ACCESS}){           print STDOUT $hash{$contig_id}{$gene_id}{COG_ACCESS}."\t";}           else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{COG_NAME}){             print STDOUT $hash{$contig_id}{$gene_id}{COG_NAME}."\t";}             else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{COG_FUNCTION}){         print STDOUT $hash{$contig_id}{$gene_id}{COG_FUNCTION}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{COG_CATEGORY}){         print STDOUT $hash{$contig_id}{$gene_id}{COG_CATEGORY}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KOG_ACCESS}){           print STDOUT $hash{$contig_id}{$gene_id}{KOG_ACCESS}."\t";}           else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KOG_NAME}){             print STDOUT $hash{$contig_id}{$gene_id}{KOG_NAME}."\t"; }            else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KOG_FUNCTION}){         print STDOUT $hash{$contig_id}{$gene_id}{KOG_FUNCTION}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{KOG_CATEGORY}){         print STDOUT $hash{$contig_id}{$gene_id}{KOG_CATEGORY}."\t";}         else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{UBLASTP_NR}){           print STDOUT $hash{$contig_id}{$gene_id}{UBLASTP_NR}."\t";}           else{ print STDOUT "NULL\t";}
      if(exists $hash{$contig_id}{$gene_id}{UBLASTP_NR_DESC}){      print STDOUT $hash{$contig_id}{$gene_id}{UBLASTP_NR_DESC}."\t";}      else{ print STDOUT "NULL\t";}
      # Check if there is only one gene on contig (1 contig and 1 gene). If so and if taxonomy values of contigs are all to set to NULL, assign the gene taxonomy to contig taxonomy
      # TODO only output consensus taxonomy on the fly (directly in this script).
      #my $size = scalar keys %{ $hash_tax{$contig_id} };
      #print STDERR "size1: $size         ";
      #$size = $size - 7;
      #print STDERR "size2: $size\n";
      if($taxonomy_type eq "consensus"){
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_K}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_K}."\t";}            else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_P}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_P}."\t";}            else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_C}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_C}."\t";}            else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_O}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_O}."\t";}            else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_F}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_F}."\t";}            else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_G}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_G}."\t";}            else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{$gene_id}{TAX_S}){            print STDOUT $hash_tax{$contig_id}{$gene_id}{TAX_S}."";}              else{ print STDOUT "NULL";}
      }else{
          if(exists $hash_tax{$contig_id}{TAX_K}){                      print STDOUT $hash_tax{$contig_id}{TAX_K}."\t";}                      else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{TAX_P}){                      print STDOUT $hash_tax{$contig_id}{TAX_P}."\t";}                      else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{TAX_C}){                      print STDOUT $hash_tax{$contig_id}{TAX_C}."\t";}                      else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{TAX_O}){                      print STDOUT $hash_tax{$contig_id}{TAX_O}."\t";}                      else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{TAX_F}){                      print STDOUT $hash_tax{$contig_id}{TAX_F}."\t";}                      else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{TAX_G}){                      print STDOUT $hash_tax{$contig_id}{TAX_G}."\t";}                      else{ print STDOUT "NULL\t";}
          if(exists $hash_tax{$contig_id}{TAX_S}){                      print STDOUT $hash_tax{$contig_id}{TAX_S}."";}                        else{ print STDOUT "NULL";}
      }
      print STDOUT "\n";
   }
}
