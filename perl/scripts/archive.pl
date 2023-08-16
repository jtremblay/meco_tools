#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use File::Find;
use Cwd qw(getcwd);

my $usage=<<'ENDHERE';
NAME:
archive.pl

PURPOSE:
To archive project in tarballs. Use in the context of nearline archiving.
Necessary in the context that we have file size limit on our account
and also a limit on the total number of files.
Should run cleanup.pl first.

INPUT:
--indir <string> : Root directory of project
--type <string>  : Type of data. 'shotgun' or 'amplicons'
				
OUTPUT:
STDOUT           : Commands to delete files. Can be sbacthed (for SLURM)...

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $type);
my $verbose = 0;

GetOptions(
   'indir=s' 	=> \$indir,
   'type=s'     => \$type,
   'verbose' 	=> \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

my %hash;

print STDOUT "#!/bin/bash\n";
print STDOUT "#SBATCH --time=10:00:00\n";
print STDOUT "#SBATCH --nodes=1\n";
print STDOUT "#SBATCH --account=rrg-cgreer\n";
print STDOUT "#SBATCH -n 1\n";
print STDOUT "#SBATCH --mem=6000\n";
print STDOUT "#SBATCH -o ./stdout_archive.txt\n";
print STDOUT "#SBATCH -e ./stderr_archive.txt\n\n";

if($type eq "shotgun"){
    print STDOUT "tar -chvf raw_reads.tar $indir/raw_reads\n";
    print STDOUT "tar -chvf qced_reads.tar $indir/qced_reads\n";
    print STDOUT "tar -zchvf assembly.tar.gz $indir/assembly\n";
    if(-d "$indir/gene_annotation"){
        print STDOUT "tar -chvf gene_annotation.tar $indir/gene_annotation\n";
    }elsif(-d "$indir/annotations"){
        print STDOUT "tar -chvf annotations.tar $indir/annotations\n";
    }

    if(-d "$indir/misc"){
        print STDOUT "tar -zchvf misc.tar.gz $indir/misc\n";
    }
    print STDOUT "tar -zchvf gene_prediction.tar.gz $indir/gene_prediction\n";
    print STDOUT "tar -zchvf gene_abundance.tar.gz $indir/gene_abundance\n";
    if(-d "$indir/contigs_abundance"){
        print STDOUT "tar -zchvf contigs_abundance.tar.gz $indir/contigs_abundance\n";
    }elsif(-d "$indir/contig_abundance"){
        print STDOUT "tar -zchvf contig_abundance.tar.gz $indir/contig_abundance\n";
    }
    print STDOUT "tar -zchvf config_files.tar.gz $indir/*.sh $indir/*.tsv $indir/*.ini $indir/*.txt\n";

}elsif($type eq "amplicons"){

    my $prefix = ""; 
    if(-d "$indir/rrnatagger"){
        $prefix = "rrnatagger";
    }elsif(-d "$indir/amplicontagger"){
        $prefix = "amplicontagger";
    }else{
        die "Can't find rrnatagger or amplicontagger dir\n";
    }
   
    my $cmd = "tar -zchvf $prefix.tar.gz";
    $cmd .= " --exclude='*ncontam.fastq'";
    $cmd .= " --exclude='*contam.fastq'";
    $cmd .= " --exclude='*ncontam_phix.fastq'";
    $cmd .= " --exclude='*ncontam_nphix.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_1.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_2.fastq'";
    $cmd .= " --exclude='*ncontam_nphix.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_paired.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_trimmed_1.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_trimmed_2.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_unpaired_1.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_unpaired_2.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_trimmed.extendedFrags.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_trimmed.notCombined_1.fastq'";
    $cmd .= " --exclude='*ncontam_nphix_trimmed.notCombined_2.fastq'";
    $cmd .= " --exclude='*obs/demultiplexed'";
    $cmd .= " --exclude='*obs/dada2/filtered'";
    $cmd .= " --exclude='*_complete.tar.gz'";
    $cmd .= " -C $indir/ ./$prefix";

    print STDOUT $cmd."\n";
  
    $cmd  = "if [ -d '$indir/raw_reads' ]; then\n";
    $cmd .= "    echo '$indir/raw_reads does exist.'\n";
    $cmd .= "    tar -chvf raw_reads.tar -C $indir/ ./raw_reads\n";
    $cmd .= "fi\n";
    print STDOUT $cmd."\n";

    $cmd = "find $indir -name '*_complete.tar.gz' -exec cp {} ./ \\;";
    print STDOUT $cmd."\n";
    
    $cmd = "find $indir -name '*export*.tar.gz' -exec cp {} ./ \\;";
    print STDOUT $cmd."\n";
    
    $cmd = "find $indir -maxdepth 1 -name 'mapping_file*' -exec cp {} ./ \\;";
    print STDOUT $cmd."\n";
    
    $cmd = "find $indir -maxdepth 1 -name 'barcodes*' -exec cp {} ./ \\;";
    print STDOUT $cmd."\n";
    
    $cmd = "find $indir -maxdepth 1 -name '*.ini' -exec cp {} ./ \\;";
    print STDOUT $cmd."\n";
}

exit;
