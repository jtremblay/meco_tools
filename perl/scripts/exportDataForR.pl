#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use File::Spec;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
exportDataForR.pl

PURPOSE:
Generate a tar ball containing all key files needed for 
downstream analysis done with R.

INPUT:
--indir <string>      : Directory
--type <string>       : either 'amplicons' or 'shotgun'
--gz                  : if flag set, will look for .gz. much of the files are
                        gzipped if cleanup.pl is run before.
--reads_type <string> : Default = reads_12. Can be reads_1 or reads_12
--prefix <string>     : Defaut = "". add prefix to archive name.
--bins <string>       : Default = "true". If true, will generate symlinks for bins results.
--table_prefix        : feature_table or otu_table

OUTPUT:

NOTES:
Just cd to a root metagenome project directory and execute like this: 
exportDataForR.pl --indir ./

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $type, $gz, $reads_type, $prefix, $bins, $table_prefix);
my $verbose = 0;

GetOptions(
   'indir=s'        => \$indir,
   'type=s'         => \$type,
   'gz'             => \$gz,
   'bins=s'         => \$bins,
   'reads_type=s'   => \$reads_type,
   'prefix=s'       => \$prefix,
   'table_prefix=s' => \$table_prefix,
   'verbose'        => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
sub findNonRecursiveParsedBins(){
	my $currIndir = "./binning/metabat2/parsed_bins";
	my $currIndir2 = $indir;
	my $extension = "fa";
    my @filenames = glob($currIndir."/*.".$extension);
    foreach my $filename (@filenames){
        print STDERR $filename."\n";    
        if (-e $filename) {
            my $basename = basename($filename);
            next if($basename  eq "all_bins.fa");
            #print STDERR substr($basename, -(length($extension)+1))."\n";
    
           if(substr($basename, -(length($extension)+1)) eq ".".$extension){
                #print STDERR "Parsing ".$filename." into single fasta file archive...\n";
                #my $outfile = "$outdir/$basename";
                #print STDERR $outfile."\n";
                #my $symlink_exists = eval { symlink("../../../binning/metabat/parsed_bins/$basename", "$currIndir2/export/parsed_bins/fasta/$basename"); 1 };
                system("ln -f -s -r -t $currIndir2/export/bins_metabat2/parsed_bins/fasta/ binning/metabat2/parsed_bins/$basename");
            }
        }
    }   
}

sub findNonRecursiveBinsMetabat(){
	my $currIndir = "./binning/metabat2/";
	my $currIndir2 = $indir;
	my $extension = "fa";
    my @filenames = glob($currIndir."/*.".$extension);
    foreach my $filename (@filenames){
        print STDERR $filename."\n";    
        if (-e $filename) {
            my $basename = basename($filename);
            next if($basename  eq "all_bins.fa");
            #print STDERR substr($basename, -(length($extension)+1))."\n";
    
           if(substr($basename, -(length($extension)+1)) eq ".".$extension){
                #print STDERR "Parsing ".$filename." into single fasta file archive...\n";
                #my $outfile = "$outdir/$basename";
                #print STDERR $outfile."\n";
                #my $symlink_exists = eval { symlink("../../../binning/metabat2/$basename", "$currIndir2/export/bins_metabat2/fasta/$basename"); 1 };
                system("ln -f -s -r -t $currIndir2/export/bins_metabat2/fasta/ binning/metabat2/$basename");
            }
        }
    }   
}

sub findNonRecursiveBinsMaxbin(){
	my $currIndir = "./binning/maxbin2/";
	my $currIndir2 = $indir;
	my $extension = "fa";
    my @filenames = glob($currIndir."/*.".$extension);
    foreach my $filename (@filenames){
        print STDERR $filename."\n";    
        if (-e $filename) {
            my $basename = basename($filename);
            next if($basename  eq "all_bins.fa");
            #print STDERR substr($basename, -(length($extension)+1))."\n";
    
           if(substr($basename, -(length($extension)+1)) eq ".".$extension){
                #print STDERR "Parsing ".$filename." into single fasta file archive...\n";
                #my $outfile = "$outdir/$basename";
                #print STDERR $outfile."\n";
                #my $symlink_exists = eval { symlink("../../../binning/maxbin2/$basename", "$currIndir2/export/bins_maxbin2/fasta/$basename"); 1 };
                system("ln -f -s -r -t $currIndir2/export/bins_maxbin2/fasta/ binning/maxbin2/$basename $basename");
            }
        }
    }   
}

$bins = "true" unless($bins);
if(!defined($reads_type)){
    $reads_type = "reads_12";
}else{
    die "--reads_type has to be 'reads_12' or 'reads_1'\n" unless($reads_type eq "reads_12" || $reads_type eq "reads_1" || $reads_type eq "reads_2");
}
if($prefix){
    $prefix = "_".$prefix;
}else{
    $prefix = "";
}

if(!defined($table_prefix)){
    $table_prefix = "feature_table";
}

my $symlink_exists;

my $indir2 = File::Spec->abs2rel($indir);
if($type eq "shotgun"){

    # First create symlinks.
    $indir = abs_path($indir);
    $indir = File::Spec->abs2rel($indir);
    system("mkdir -p $indir/export/");
    system("mkdir -p $indir/export/consensus/relative");
    system("mkdir -p $indir/export/consensus/absolute");
    #system("mkdir -p $indir/export/parsed_bins/fasta");
    system("mkdir -p $indir/export/genes");
    system("mkdir -p $indir/export/contigs");
    if($bins eq "true"){
        system("test -f $indir/binning/metabat2/link.tsv  &&                         ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/binning/metabat2/link.tsv");
        system("test -f $indir/binning/metabat2/raw_bins.tsv &&                      ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/binning/metabat2/raw_bins.tsv");                     
        system("test -f $indir/binning/metabat2/out_checkm.txt &&                    ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/binning/metabat2/out_checkm.txt");                   
        system("test -f $indir/binning/metabat2/parsed_bins.tsv &&                   ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/binning/metabat2/parsed_bins.tsv");                  
        system("test -f $indir/binning/metabat2/parsed_bins/out_checkm.txt &&        ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/binning/metabat2/parsed_bins/out_checkm.txt        out_checkm_parsed_bins.txt");       
        system("test -f $indir/annotations/metabat2/bins/feature_table.tsv &&        ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/annotations/metabat2/bins/feature_table.tsv");       
        system("test -f $indir/annotations/metabat2/parsed_bins/feature_table.tsv && ln -f -s -r -t $indir2/export/bins_metabat2/       $indir/annotations/metabat2/parsed_bins/feature_table.tsv feature_table_parsed_bins.tsv");
                                                                                                            
        system("test -f $indir/binning/maxbin2/link.tsv &&                           ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/binning/maxbin2/link.tsv");                          
        system("test -f $indir/binning/maxbin2/raw_bins.tsv &&                       ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/binning/maxbin2/raw_bins.tsv");                      
        system("test -f $indir/binning/maxbin2/out_checkm.txt &&                     ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/binning/maxbin2/out_checkm.txt");                    
        system("test -f $indir/binning/maxbin2/parsed_bins.tsv &&                    ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/binning/maxbin2/parsed_bins.tsv");                   
        system("test -f $indir/binning/maxbin2/parsed_bins/out_checkm.txt &&         ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/binning/maxbin2/parsed_bins/out_checkm.txt         out_checkm_parsed_bins.txt");        
        system("test -f $indir/annotations/maxbin2/bins/feature_table.tsv &&         ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/annotations/maxbin2/bins/feature_table.tsv");        
        system("test -f $indir/annotations/maxbin2/parsed_bins/feature_table.tsv &&  ln -f -s -r -t $indir2/export/bins_maxbin2/        $indir/annotations/maxbin2/parsed_bins/feature_table.tsv  feature_table_parsed_bins.tsv"); 
    
        system("mkdir -p $indir/export/bins_metabat2/fasta");
        system("mkdir -p $indir/export/bins_maxbin2/fasta");

        findNonRecursiveBinsMaxbin();
        findNonRecursiveBinsMetabat();
        
   }
   if($gz){
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_bins/".$table_prefix.".tsv.gz",                                 "$indir2/export/bins/".$table_prefix.".tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_parsed_bins/".$table_prefix.".tsv.gz",                          "$indir2/export/parsed_bins/".$table_prefix.".tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/".$table_prefix."_normalized_bacteriaArchaea.tsv.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/taxonomy.tsv.gz",                                     "$indir2/export/consensus/taxonomy.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/annotations.tsv.gz",                                                     "$indir2/export/annotations.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/kegg_matrix_modules.tsv.gz",                                             "$indir2/export/kegg_matrix_modules.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/kegg_modules.tsv.gz",                                                    "$indir2/export/kegg_modules.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/gene_abundance/merged_gene_abundance_cpm.tsv.gz",                              "$indir2/export/genes/merged_gene_abundance_cpm.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/binning/metabat/parsed_bins/link.tsv.gz",                                      "$indir2/export/bins/link.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/mapping_file.tsv",                                                             "$indir2/export/mapping_file.tsv"); 1 };
      $symlink_exists = eval { symlink("$indir/contig_abundance/qc_mapping_stats.tsv",                                       "$indir2/export/qc_mapping_stats.tsv"); 1 };
      $symlink_exists = eval { symlink("$indir/assembly/ray_assembly_stats.txt",                                              "$indir2/export/ray_assembly_stats.txt"); 1 };
      #$symlink_exists = eval { symlink("$indir/betadiv/bray_curtis_".$table_prefix."_normalized.txt.gz",                             "$indir2/export/betadiv/bray_curtis_".$table_prefix."_normalized.txt.gz"); 1 };
      #$symlink_exists = eval { symlink("$indir/betadiv/bray_curtis_".$table_prefix."_normalized_coords.tsv.gz",                      "$indir2/export/betadiv/bray_curtis_".$table_prefix."_normalized_coords.tsv.gz"); 1 };
      #$symlink_exists = eval { symlink("$indir/betadiv/bray_curtis_".$table_prefix."_normalized_bacteriaArchaea.txt.gz",             "$indir2/export/betadiv/bray_curtis_".$table_prefix."_normalized_bacteriaArchaea.txt.gz"); 1 };
      #$symlink_exists = eval { symlink("$indir/betadiv/bray_curtis_".$table_prefix."_normalized_bacteriaArchaea_coords.tsv.gz",      "$indir2/export/betadiv/bray_curtis_".$table_prefix."_normalized_bacteriaArchaea_coords.tsv.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/beta_div",                                                                   "$indir2/export/beta_div"); 1 };
      $symlink_exists = eval { symlink("$indir/alpha_div",                                                                   "$indir2/export/alpha_div"); 1 };
      #$symlink_exists = eval { symlink("$indir/alphadiv/bray_curtis_".$table_prefix."_normalized_bacteriaArchaea/collated",           "$indir2/export/alphadiv/bray_curtis_".$table_prefix."_normalized_bacteriaArchaea"); 1 };
      #$symlink_exists = eval { symlink("$indir/alphadiv/bray_curtis_".$table_prefix."_normalized/collated",                           "$indir2/export/alphadiv/bray_curtis_".$table_prefix."_normalized"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L1.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L1.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L2.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L2.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L3.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L3.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L4.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L4.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L5.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L5.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L6.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L6.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L7.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_bacteriaArchaea_L7.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L1.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L1.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L2.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L2.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L3.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L3.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L4.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L4.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L5.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L5.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L6.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L6.txt.gz"); 1 };
      $symlink_exists = eval { symlink("$indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L7.txt.gz", "$indir2/export/consensus/".$table_prefix."_normalized_L7.txt.gz"); 1 };
   }else{
       system("test -f $indir/annotations/taxonomy_consensus/".$table_prefix."_normalized_bacteriaArchaea.tsv && ln -f -s -r -t $indir2/export/consensus/       $indir/annotations/taxonomy_consensus/".$table_prefix."_normalized_bacteriaArchaea.tsv"   );
       system("test -f $indir/annotations/taxonomy_consensus/".$table_prefix."_bacteriaArchaea.tsv && ln -f -s -r -t            $indir2/export/consensus/       $indir/annotations/taxonomy_consensus/".$table_prefix."_bacteriaArchaea.tsv"   );
       system("test -f $indir/annotations/taxonomy_consensus/".$table_prefix."_normalized.tsv && ln -f -s -r -t                 $indir2/export/consensus/       $indir/annotations/taxonomy_consensus/".$table_prefix."_normalized.tsv"   );
       system("test -f $indir/annotations/taxonomy_consensus/".$table_prefix.".tsv && ln -f -s -r -t                            $indir2/export/consensus/       $indir/annotations/taxonomy_consensus/".$table_prefix.".tsv"   );
       system("test -f $indir/annotations/taxonomy_consensus/taxonomy.tsv && ln -f -s -r -t                                     $indir2/export/consensus/       $indir/annotations/taxonomy_consensus/taxonomy.tsv"   );
       system("test -f $indir/annotations/annotations.tsv && ln -f -s -r -t                                                     $indir2/export/                 $indir/annotations/annotations.tsv"   );
       system("test -f $indir/gene_abundance/merged_gene_abundance_cpm.tsv && ln -f -s -r -t                                    $indir2/export/genes/           $indir/gene_abundance/merged_gene_abundance_cpm.tsv"   );
       system("test -f $indir/gene_abundance/merged_gene_abundance.tsv && ln -f -s -r -t                                        $indir2/export/genes/           $indir/gene_abundance/merged_gene_abundance.tsv"   );
       system("test -f $indir/contig_abundance/merged_contig_abundance_cpm.tsv && ln -f -s -r -t                                $indir2/export/contigs/         $indir/contig_abundance/merged_contig_abundance_cpm.tsv"   );
       system("test -f $indir/contig_abundance/merged_contig_abundance.tsv && ln -f -s -r -t                                    $indir2/export/contigs/         $indir/contig_abundance/merged_contig_abundance.tsv"   );
       system("test -f $indir/contigs_abundance/merged_contigs_abundance_cpm.tsv && ln -f -s -r -t                              $indir2/export/contigs/         $indir/contigs_abundance/merged_contigs_abundance_cpm.tsv"   );
       system("test -f $indir/contigs_abundance/merged_contigs_abundance.tsv && ln -f -s -r -t                                  $indir2/export/contigs/         $indir/contigs_abundance/merged_contigs_abundance.tsv"   );
       system("test -f $indir/mapping_file.tsv && ln -f -s -r -t                                                                $indir2/export/                 $indir/mapping_file.tsv"   );
       system("test -f $indir/contig_abundance/qc_mapping_stats.tsv && ln -f -s -r -t                                           $indir2/export/contigs/         $indir/contig_abundance/qc_mapping_stats.tsv"   );
       system("test -f $indir/contigs_abundance/qc_mapping_stats.tsv && ln -f -s -r -t                                          $indir2/export/contigs/         $indir/contigs_abundance/qc_mapping_stats.tsv"   );
       system("test -f $indir/assembly/assembly_stats.txt && ln -f -s -r -t                                                     $indir2/export/contigs/         $indir/assembly/assembly_stats.txt"   );
       system("test -f $indir/assembly/Contigs.fasta && ln -f -s -r -t                                                          $indir2/export/contigs/         $indir/assembly/Contigs.fasta"   );
       system("test -f $indir/gene_prediction/Contigs_renamed.fna && ln -f -s -r -t                                             $indir2/export/genes/           $indir/gene_prediction/Contigs_renamed.fna && mv $indir2/export/genes/Contigs_renamed.fna $indir2/export/genes/genes.fna");
       system("test -f $indir/gene_prediction/Contigs_renamed.faa && ln -f -s -r -t                                             $indir2/export/genes/           $indir/gene_prediction/Contigs_renamed.faa && mv $indir2/export/genes/Contigs_renamed.faa $indir2/export/genes/genes.faa");
       system("test -f $indir/assembly/Contigs_tf5mer_clr_bhtsne.tsv && ln -f -s -r -t                                          $indir2/export/contigs/         $indir/assembly/Contigs_tf5mer_clr_bhtsne.tsv"   );
       system("test -e $indir/beta_div && ln -f -s -r -t                                                                        $indir2/export/                 $indir/beta_div"   );
       system("test -e $indir/alpha_div && ln -f -s -r -t                                                                       $indir2/export/                 $indir/alpha_div"   );
       system("test -e $indir/assembly/quast && ln -f -s -r -t                                                                  $indir2/export/contigs/         $indir/assembly/quast"   );
      
       ####### ABSOLUTE
       # Normalized
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L1.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L2.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L3.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L4.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L5.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L6.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L7.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_normalized_bacteriaArchaea_L7.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L1.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L1.tsv");                           
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L2.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L2.tsv");                           
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L3.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L3.tsv");                           
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L4.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L4.tsv");                           
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L5.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L5.tsv");                           
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L6.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L6.tsv");                           
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L7.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_normalized_L7.tsv");                           
      
       # Raw (non-normalized)
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L1.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L2.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L3.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L4.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L5.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L6.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L7.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/absolute/".$table_prefix."_bacteriaArchaea_L7.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L1.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L2.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L3.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L4.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L5.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L6.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L7.tsv && ln -f -s -r -t $indir2/export/consensus/absolute/ $indir/annotations/taxonomy_consensus/all/absolute/".$table_prefix."_L7.tsv");

       # RELATIVE
       # Normalized 
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L1.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L2.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L3.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L4.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L5.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L6.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L7.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_normalized_bacteriaArchaea_L7.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L1.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L2.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L3.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L4.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L5.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L6.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L7.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_normalized_L7.tsv");
      
       # Raw (non-Normalized)
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L1.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L2.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L3.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L4.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L5.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L6.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L7.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/bacteriaArchaea/relative/".$table_prefix."_bacteriaArchaea_L7.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L1.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L1.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L2.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L2.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L3.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L3.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L4.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L4.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L5.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L5.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L6.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L6.tsv");
       system("test -f $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L7.tsv && ln -f -s -r -t $indir2/export/consensus/relative/ $indir/annotations/taxonomy_consensus/all/relative/".$table_prefix."_L7.tsv");
   }
   #my @array = (
   #   "$indir2/export/bins/otu_table_bacteriaArchaea.tsv",
   #   "$indir2/export/contigs/otu_table_normalized_bacteriaArchaea.tsv",
   #   "$indir2/export/contigs/taxonomy.tsv",
   #   "$indir2/export/annotations.tsv",
   #   "$indir2/export/genes/merged_gene_abundance_cpm.tsv",
   #   "$indir2/export/bins/link.tsv",
   #   "$indir2/export/mapping_file.tsv",
   #   "$indir2/export/contigs/otu_table_normalized_bacteriaArchaea_L6.txt",
   #   "$indir2/export/contigs/otu_table_normalized_bacteriaArchaea_L7.txt"
   #);
   print STDERR "Archiving directory with command : tar -zchvf ./export$prefix.tar.gz ./export/\n";
   system("tar -zchvf ./export$prefix.tar.gz ./export/");

}elsif($type eq "amplicons"){
   # First create symlinks.
   my $indir2 = $indir;
   $indir = abs_path($indir);
   $indir = File::Spec->abs2rel($indir);
   system("mkdir -p $indir/export/");
   system("mkdir -p $indir/export/betadiv");
   system("mkdir -p $indir/export/blast");
   system("mkdir -p $indir/export/alphadiv");
   system("mkdir -p $indir/export/taxonomy");
   system("mkdir -p $indir/export/$table_prefix"."s");
   #system("mkdir -p $indir/export/otu_tables");
   #system("mkdir -p $indir/export/picrust");
   system("mkdir -p $indir/export/qscores");
   system("mkdir -p $indir/export/tree");
   system("mkdir -p $indir/export/rdp");
   if(-e "$indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L7.tsv"){
       #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L7.tsv",       "$indir2/export/taxonomy/".$table_prefix."_final_normalized_L7.tsv"); 1 };
      system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L7.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L7.tsv");
   }
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L6.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L6.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L5.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L5.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L4.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L4.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L3.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L3.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L2.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L2.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L1.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_normalized_L1.tsv");
   
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L6.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L6.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L5.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L5.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L4.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L4.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L3.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L3.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L2.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L2.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L1.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/absolute/".$table_prefix."_final_rarefied_L1.tsv");
   
   if(-e "$indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L7.tsv"){
      system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L7.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_rel_L7.tsv");
   }
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L6.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L6.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L5.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L5.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L4.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L4.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L3.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L3.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L2.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L2.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L1.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_normalized_L1.tsv");
   
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L6.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L6.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L5.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L5.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L4.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L4.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L3.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L3.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L2.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L2.tsv");
   system("test -f $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L1.tsv && ln -f -s -r -t $indir2/export/taxonomy/ $indir/amplicontagger/$reads_type/tax_summary/relative/".$table_prefix."_final_rarefied_L1.tsv");
   
   system("test -f $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_normalized.tsv         && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_normalized.tsv");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_normalized_coords.tsv  && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_normalized_coords.tsv");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_normalized.tsv              && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_normalized.tsv");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_normalized_coords.tsv       && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_normalized_coords.tsv");

   system("test -f $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_rarefied.tsv           && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_rarefied.tsv");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_rarefied_coords.tsv    && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/weighted_unifrac_".$table_prefix."_final_rarefied_coords.tsv");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_rarefied.tsv                && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_rarefied.tsv");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_rarefied_coords.tsv         && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/bray_curtis_".$table_prefix."_final_rarefied_coords.tsv");

   system("test -f $indir/amplicontagger/$reads_type/beta_div/plots/3d_weighted_unifrac                            && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/3d_weighted_unifrac");
   system("test -f $indir/amplicontagger/$reads_type/beta_div/plots/3d_bray_curtis                                 && ln -f -s -r -t $indir2/export/betadiv/ $indir/amplicontagger/$reads_type/beta_div/3d_bray_curtis");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_normalized.tsv     && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_normalized.tsv");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_normalized.biom    && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_normalized.biom");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_rarefied.tsv       && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_rarefied.tsv");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_rarefied.biom      && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_final_rarefied.biom");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered.tsv             && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered.tsv");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_normalized.tsv  && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_normalized.tsv");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered.biom            && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered.biom");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_normalized.biom && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_normalized.biom");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_rarefied.tsv    && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_rarefied.tsv");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_rarefied.biom   && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_filtered_rarefied.biom");
   system("test -f $indir/amplicontagger/$reads_type/".$table_prefix."s/features.fasta                             && ln -f -s -r -t $indir2/export/".$table_prefix."s/ $indir/amplicontagger/$reads_type/".$table_prefix."s/features.fasta");

   system("test -f $indir/mapping_file.tsv                                         && ln -f -s -r -t $indir2/export/          $indir/mapping_file.tsv");
   system("test -f $indir/amplicontagger/$reads_type/alpha_div/alpha_richness.tsv  && ln -f -s -r -t $indir2/export/alphadiv/ $indir/amplicontagger/$reads_type/alpha_div/alpha_richness.tsv");
   system("test -f $indir/amplicontagger/$reads_type/alpha_div/alpha_shannon.tsv   && ln -f -s -r -t $indir2/export/alphadiv/ $indir/amplicontagger/$reads_type/alpha_div/alpha_shannon.tsv");
   system("test -f $indir/amplicontagger/$reads_type/alpha_div/alpha_simpson.tsv   && ln -f -s -r -t $indir2/export/alphadiv/ $indir/amplicontagger/$reads_type/alpha_div/alpha_simpson.tsv");
   system("test -f $indir/amplicontagger/$reads_type/alpha_div/alpha_chao1.tsv     && ln -f -s -r -t $indir2/export/alphadiv/ $indir/amplicontagger/$reads_type/alpha_div/alpha_chao1.tsv");
   system("test -f $indir/amplicontagger/$reads_type/countReport.tsv               && ln -f -s -r -t $indir2/export/          $indir/amplicontagger/$reads_type/countReport.tsv");
   system("test -f $indir/amplicontagger/$reads_type/blast/blast.out.besthit       && ln -f -s -r -t $indir2/export/blast/    $indir/amplicontagger/$reads_type/blast/blast.out.besthit");
   system("test -f $indir/amplicontagger/$reads_type/blast/blast.out               && ln -f -s -r -t $indir2/export/blast/    $indir/amplicontagger/$reads_type/blast/blast.out");
   system("test -f $indir/amplicontagger/$reads_type/rdp/rdp.tsv                   && ln -f -s -r -t $indir2/export/rdp/      $indir/amplicontagger/$reads_type/rdp/rdp.tsv");
   system("test -f $indir/amplicontagger/$reads_type/fasttree/tree.fasttree        && ln -f -s -r -t $indir2/export/tree/     $indir/amplicontagger/$reads_type/fasttree/tree.fasttree");

   # Picrust
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes.tsv",                                  "$indir2/export/picrust/".$table_prefix."_greengenes.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_precalculated.tsv",                    "$indir2/export/picrust/".$table_prefix."_greengenes_precalculated.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups_nbcn_picrust_KO_sum.tsv",       "$indir2/export/picrust/".$table_prefix."_greengenes_nodups_nbcn_picrust_KO_sum.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups_nbcn_picrust_KO_sum.biom",      "$indir2/export/picrust/".$table_prefix."_greengenes_nodups_nbcn_picrust_KO_sum.biom"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups_nbcn_picrust_KO_occurence.tsv", "$indir2/export/picrust/".$table_prefix."_greengenes_nodups_nbcn_picrust_KO_occurence.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups_nbcn.tsv",                      "$indir2/export/picrust/".$table_prefix."_greengenes_nodups_nbcn.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups_nbcn.biom",                     "$indir2/export/picrust/".$table_prefix."_greengenes_nodups_nbcn.biom"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups.biom",                          "$indir2/export/picrust/".$table_prefix."_greengenes_nodups.biom"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/".$table_prefix."_greengenes_nodups.tsv",                           "$indir2/export/picrust/".$table_prefix."_greengenes_nodups.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/".$table_prefix."s/link.tsv",                                                          "$indir2/export/picrust/link.tsv"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/qual_stats_unassembled.pdf",                                                           "$indir2/export/qscores/qual_stats_unassembled.pdf"); 1 };
   #$symlink_exists = eval { symlink("$indir/amplicontagger/$reads_type/qscores/qscores_assembled_QCpassed.pdf",                                               "$indir2/export/qscores/qscores_assembled_QCpassed.pdf"); 1 };

   system("tar -zchvf ./export$prefix.tar.gz ./export/");
}else{
   die("--type has to be either 'amplicons' or 'shotgun'\n");
}

