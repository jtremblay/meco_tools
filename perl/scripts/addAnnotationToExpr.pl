#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
addAnnotationToExpr.pl

PURPOSE:

INPUT:
--infile <string> : annotation tabular file.
--indir <string>  : indir where your DDA directory structure is.
--prefix <string> : prefix. will search for files having this prefix.
				
OUTPUT:
Nothing, will find output of edger and generate a file with _annotated.tsv extension.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $indir);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'indir=s'   => \$indir,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--infile missing\n" unless($infile);
die "--indir missing\n" unless($indir);

my %hash_edger;
my %hash_expr;

## SUBS
sub eachFile{
   my $filename = $_;
   my $fullpath = $File::Find::name;
   #remember that File::Find changes your CWD, 
   #so you can call open with just $_

   if (-e $filename) { 
      
      if($filename =~ m/_edger\.tsv$/){
         my $path = dirname($fullpath);
         $hash_edger{$path}{fc} = $filename;
      }
      if($filename =~ m/_normalized_significant\.tsv$/){
         my $path = dirname($fullpath);
         $hash_edger{$path}{expr} = $filename;
      }
   }
}

## MAIN

# Populate hash with file names.
find (\&eachFile, $indir);

# Populate annotation hash
my %hash;
my $number_of_el;
my @header_ann;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   if($_ =~ m/^#/){
      my @row = split(/\t/, $_);
      shift(@row); shift(@row);
      @header_ann = @row;
      next;
   }
   my @row = split(/\t/, $_);
   my $gene_id = $row[1];
   shift(@row); shift(@row);
   my $row = join("\t", @row);
   $hash{$gene_id} = $row;
   $number_of_el = scalar(@row);
}
close(IN);

# then loop through all hash_edger and hash_expr and write outfile (which is basically files in hash_edger with gene annotations.
for my $key (%hash_edger){
   if(defined $hash_edger{$key}{expr} and defined $hash_edger{$key}{fc}){
      my $fcFile = $key."/".$hash_edger{$key}{fc};
      my $exprFile = $key."/".$hash_edger{$key}{expr};
      print STDERR "Processing ".$fcFile."\n";
      print STDERR "Processing ".$exprFile."\n";
      my $outfile = $fcFile;
      #if($outfile =~ m/_normalized_significant\.tsv/){
      #   $outfile =~ s/_normalized_significant\.tsv/_normalized_significant_annotated\.tsv/;
      #}else{
      #   die "file suffix/ext does not match 'edger.tsv'\n";
      #}
      if($outfile =~ m/_edger\.tsv/){
         $outfile =~ s/_edger\.tsv/_edger_annotated\.tsv/;
      }else{
         die "file suffix/ext does not match 'edger.tsv'\n";
      }
      open(OUT, ">".$outfile) or die "Can't open file $outfile\n";

      
      my %hash_fc;
      my %hash_expr;
      my @header;

      open(IN, "<".$fcFile) or die "Can't open $fcFile\n";
      while(<IN>){
         chomp;
         if($. == 1){
            @header = split(/\t/, $_);
            shift(@header);
            next;
         }
         my @row = split(/\t/, $_);
         my $gene_id = shift(@row);
         my $values = join("\t", @row);
         $hash_fc{$gene_id} = $values;
      }
      close(IN);
      
      open(IN, "<".$exprFile) or die "Can't open $exprFile\n";
      while(<IN>){
         chomp;
         if($. == 1){
            my @row = split(/\t/, $_);
            shift(@row);
            push @header, @row;
            push @header, @header_ann;
            next;
         }
         my @row = split(/\t/, $_);
         my $gene_id = shift(@row);
         my $values = join("\t", @row);
         $hash_expr{$gene_id} = $values;
      }
      close(IN);
      
      print OUT "#gene_id\t"; print OUT join("\t", @header); print OUT "\n";
      for my $key2 (keys %hash_expr){
         print OUT $key2."\t";
         print OUT $hash_fc{$key2}."\t";
         print OUT $hash_expr{$key2}."\t";
         if(exists $hash{$key2}){
            print OUT $hash{$key2}."\n";
         }else{
            print OUT "NULL";
            for(my $k=0;$k<=$number_of_el;$k++){ print "\tNULL"; }
         }
      }
      close(OUT);
   }
}

   
