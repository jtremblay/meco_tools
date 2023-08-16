#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Env qw(TMPDIR);
use File::Temp;

my $usage=<<'ENDHERE';
NAME:
compareCategories.pl

PURPOSE:

INPUT:
--infile <string>       : distance matrix (i.e. for instance, a betadiversity distance matrix)
--mapping_file <string> : Qiime formatted mapping file

OUTPUT:
--outdir <string>  : Where text files containing all statistical results
                     from compare_categories.py (Qiime) will be written

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $outdir, $mapping_file);
my $verbose = 0;

GetOptions(
   'infile=s' 	     => \$infile,
   'outdir=s'       => \$outdir,
   'mapping_file=s' => \$mapping_file,
   'verbose' 	     => \$verbose,
   'help'           => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $tmpdir = File::Temp->newdir(
    "tmpdir-compareCategories-XXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

open(MAP,"<".$mapping_file) or die "Can't open mapping file\n";
my $header;
while(<MAP>){
   chomp;
   $header = $_;
   last;
}
close(MAP);

my @factors = split(/\t/, $header);
shift(@factors);

#compare_categories.py --method adonis -i unweighted_unifrac_dm.txt -m Fasting_Map.txt -c Treatment -o adonis_out -n 999
my @cmds;
foreach(@factors){
   my $variable = $_;
   print STDERR $variable."\n" if($verbose);
   #print STDERR "memtime compare_categories.py --method adonis -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/adonis_$variable/\n";
   system("memtime compare_categories.py --method adonis -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/adonis_$variable/");
   system("memtime compare_categories.py --method anosim -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/anosim_$variable/");
   system("memtime compare_categories.py --method mrpp -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/mrpp_$variable/");
   system("memtime compare_categories.py --method permanova -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/permanova_$variable/");
   system("memtime compare_categories.py --method permdisp -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/permdisp_$variable/");
   system("memtime compare_categories.py --method dbrda -i $infile -m $mapping_file -c $variable -n 999 -o $tmpdir/dbrda_$variable/");

   system("echo  > $outdir/$variable.txt");
   system("echo '=================ADONIS=======================\n' >> $outdir/$variable.txt");
   system("cat  $tmpdir/adonis_$variable/*.txt >> $outdir/$variable.txt"); 
   system("echo '=================ANOSIM========================\n' >> $outdir/$variable.txt");
   system("cat  $tmpdir/anosim_$variable/*.txt >> $outdir/$variable.txt"); 
   system("echo '=================MRPP=========================\n' >> $outdir/$variable.txt");
   system("cat  $tmpdir/mrpp_$variable/*.txt >> $outdir/$variable.txt"); 
   system("echo '=================PERMANOVA====================\n' >> $outdir/$variable.txt");
   system("cat  $tmpdir/permanova_$variable/*.txt >> $outdir/$variable.txt"); 
   system("echo '=================PERMDISP=====================\n' >> $outdir/$variable.txt");
   system("cat  $tmpdir/permdisp_$variable/*.txt >> $outdir/$variable.txt"); 
   system("echo '=================DBRDA========================\n' >> $outdir/$variable.txt");
   system("cat  $tmpdir/dbrda_$variable/*.txt >> $outdir/$variable.txt"); 
}

