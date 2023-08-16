#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use File::Find;
use Cwd qw(getcwd);

my $usage=<<'ENDHERE';
NAME:
cleanup.pl

PURPOSE:
To remove unuseful data from a project directory.
Necessary in the context that we have file size limit on our account
and also a limit on the total number of files.

INPUT:
--indir <string> : Root directory of project
--type <string>  : Type of data. 'shotgun' or 'amplicons'
--gz             : Specify flag if key files are to be gzipped.
--job_output     : Specify flag if job_output folder is to be gzipped. 
                   The original job_output folder will also be deleted.
				
OUTPUT:
STDOUT           : Commands to delete files. Can be sbacthed (for SLURM)...

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $type, $gz, $job_output);
my $verbose = 0;

GetOptions(
   'indir=s' 	=> \$indir,
   'type=s'     => \$type,
   'gz'         => \$gz,
   'job_output' => \$job_output,
   'verbose' 	=> \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

my %hash;

# Remove the following from qced_reads
$indir = abs_path($indir);
my @to_match = (
    ".*.ncontam.fastq.gz",
    ".*.contam.fastq.gz",
    ".*.trim.interleaved.fastq.gz",
    ".*.trim.single2.fastq.gz",
    ".*.trim.pair2.fastq.gz",
    ".*.trim.single1.fastq.gz",
    ".*.trim.pair1.fastq.gz",
    ".*.ncontam_R1.fastq.gz",
    ".*.ncontam_R2.fastq.gz",
    ".*.contam_R1.fastq",
    ".*.contam_R2.fastq",
    ".*.ncontam_unpaired_R1.fastq.gz",
    ".*.ncontam_unpaired_R2.fastq.gz",
    ".*.ncontam_paired_R1.fastq.gz",
    ".*.ncontam_paired_R2.fastq.gz"
);

my @to_match2 = (
    "hmmscan_chunk_.*",
    "blastn_chunk_.*",
    "blastp_chunk_.*",
    "rpsblast_chunk_.*",
    "ublastp_chunk_.*"
);

my @to_gzip = (
    "hmmscan_tigrfam_pfam_pfamtblout.tsv",
    "hmmscan_pfam_pfam_pfamtblout.tsv",
    "hmmscan_tigrfam_domtblout.tsv",
    "Contigs.fasta",
    "hmmscan_pfam_domtblout.tsv",
    "Contigs.gff",
    "hmmscan_tigrfam_tblout.tsv",
    "annotations.tsv",
    "hmmscan_pfam_tblout.tsv",
    "rpsblast_cog.tsv",
    "blastp_kegg_parsed.tsv",
    "blastn_nt.tsv",
    "blastn_nt_besthit.tsv",
    "ublastp_nr_annotated.tsv",
    "ublastp_nr.tsv",
    "rpsblast_kog.tsv",
    "blastn_nt_contigs.tsv",
    "blastp_kegg.tsv",
    "blastn_nt_contigs_besthit.tsv",
    "gc_table.tsv",
    "merged_contigs_abundance.*",
    "merged_gene_abundance.*",
    "otu_table.*",
    "link.tsv",
    "taxonomy.tsv",
    "bray_curtis_otu_table.*" 
);

my @amplicons_to_delete = (
    "alpha_rarefaction_.*",
    "rarefaction_.*",
    "ncontam_nphix_trimmed.extendedFrags.fastq",
    "ncontam_nphix_trimmed.notCombined_1.fastq",
    "ncontam_nphix_trimmed.notCombined_2.fastq",
    "ncontam.fastq",
    "ncontam_nphix.fastq",
    "ncontam_phix.fastq",
    "ncontam_nphix_1.fastq",
    "ncontam_nphix_2.fastq",
    "ncontam_nphix.fastq",
    "ncontam_nphix_trimmed_1.fastq",
    "ncontam_nphix_trimmed_2.fastq",
    "ncontam_nphix_trimmed.extendedFrags_QCfailed.fastq",
    "ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq",
    "ncontam_nphix_unpaired_1.fastq",
    "ncontam_nphix_unpaired_2.fastq"
);

my @amplicons_to_gz = (
    "derep1_099_derep2.fasta",
    "derep1_099.fasta",
    "derep1.fasta",
    "qscores_1.tsv",
    "qscores_2.tsv",
    "ncontam_nphix_paired.fastq"
);

#BiologicalAbundances
#Scaffolds.fasta


sub eachFile{
    my $filename = $_;
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_

    if (-e $filename) { 
    
        #if(substr($filename, -9) eq ".fastq.gz" || substr($filename, -6) eq ".fastq"){

        if($type eq "shotgun"){
            foreach my $to_match (@to_match){
                if($filename =~ m/$to_match/){
                    #print STDOUT "rm ".$fullpath."\n";
                    $hash{$fullpath} = "rm ".$fullpath."\n";
                }  
            }    
            foreach my $to_match (@to_match2){
                if($filename =~ m/$to_match/){
                    #print STDOUT "rm ".$fullpath."\n";
                    $hash{$fullpath} = "rm ".$fullpath."\n";
                }  
            }

            if($gz){
                foreach my $to_gzip (@to_gzip){
                    if($filename =~ m/$to_gzip/){
                        #print STDOUT "rm ".$fullpath."\n";
                        $hash{$fullpath} = "pigz -p 4 ".$fullpath."\n";
                    }  
                }   
             } 
        }elsif($type eq "amplicons"){
            foreach my $to_match (@amplicons_to_delete){
                if($filename =~ m/$to_match/){
                    #print STDOUT "rm ".$fullpath."\n";
                    $hash{$fullpath} = "rm ".$fullpath."\n";
                }  
            }    
            if($gz){
                foreach my $to_gzip (@amplicons_to_gz){
                    if($filename =~ m/$to_gzip/){
                        #print STDOUT "rm ".$fullpath."\n";
                        $hash{$fullpath} = "pigz -p 4 ".$fullpath."\n";
                    }  
                }
            }    
        }
    }
}

## MAIN

if($type eq "amplicons"){
    my $indir_amplicons = $indir."/";
    find (\&eachFile, $indir_amplicons);

}elsif($type eq "shotgun"){
    my $indir_qc = $indir."/qced_reads/";
    my $indir_annotation = $indir."/annotations/";
    my $indir_contigs_abundance = $indir."/contig_abundance/";
    my $indir_gene_abundance = $indir."/gene_abundance/";
    find (\&eachFile, $indir_qc);
    find (\&eachFile, $indir_annotation);
    find (\&eachFile, $indir_contigs_abundance);
    find (\&eachFile, $indir_gene_abundance);
}else{
    die "--type amplicons or --type shotgun only...\n";
}


# Print to STDOUT
#print STDOUT "#!/bin/bash\n";
#print STDOUT "#PBS -l nodes=1:ppn=1\n";
#print STDOUT "#PBS -q metaq\n";
#print STDOUT "#PBS -l walltime=12:00:00\n";
#print STDOUT "#PBS -A rrg-cgreer\n";
#print STDOUT "#PBS -o outputfile\n";
#print STDOUT "#PBS -e errorfile\n";
#print STDOUT "#PBS -V\n";
#print STDOUT "#PBS -N cleanup\n";
#print STDOUT "#PBS -d ".getcwd."\n\n";

print STDOUT "#!/bin/bash\n";
print STDOUT "#SBATCH --time=10:00:00\n";
print STDOUT "#SBATCH --nodes=1\n";
print STDOUT "#SBATCH --account=rrg-cgreer\n";
print STDOUT "#SBATCH -n 1\n";
print STDOUT "#SBATCH --mem=6000\n";
print STDOUT "#SBATCH -o ./stdout_cleanup.txt\n";
print STDOUT "#SBATCH -e ./stderr_cleanup.txt\n\n";

if($type eq "shotgun"){
    print STDOUT "find $indir/ -name '*_chunk_*' -exec rm {} \\;\n";
    print STDOUT "rm -rf $indir/binning/metabat2/parsed_bins/out_checkm/bins";
    print STDOUT "rm -rf $indir/binning/metabat2/parsed_bins/out_checkm/storage";
    print STDOUT "rm -rf $indir/binning/metabat2/out_checkm/bins";
    print STDOUT "rm -rf $indir/binning/metabat2/out_checkm/storage";

}

for my $key (keys %hash){
    print STDOUT $hash{$key};
}

if($job_output){
    print STDOUT "tar -zcvf $indir/job_output.tar.gz $indir/job_output && rm -r $indir/job_output\n";
}
exit;
