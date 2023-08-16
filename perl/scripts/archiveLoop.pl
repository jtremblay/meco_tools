#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use File::Find;
use Cwd qw(getcwd);
use File::Spec;
use Data::Dumper;
use Cwd;

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

## MAIN
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

my %hash;
my %hash_del;
my $extension = ".ini";
my $cwd = getcwd;

sub eachFile{
    my $filename = $_; 
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_
    #if(-d $filename and $filename =~ m/plots/) {
    #    print STDERR $fullpath."\n";
    #    next;
    #} 
    
    if (-e $filename) { 
        #if(substr($filename, -(length($extension)+1)) eq $extension ){
        if($filename =~ m/RRNATagger.*\.ini$/ || $filename =~ m/AmpliconTagger*\.ini$/){
            my($vol, $dir, $file) = File::Spec->splitpath($fullpath);
            print STDERR $fullpath."\n";
            print STDERR $filename."\n";
            print STDERR $dir."\n";


            if(-d "$dir/amplicontagger"){
                $hash{$dir} = "amplicontagger"
            }elsif(-d "$dir/rrnatagger"){
                $hash{$dir} = "rrnatagger"
            }
        }
        foreach my $to_match (@amplicons_to_delete){
            if($filename =~ m/$to_match/){
                #print STDOUT "rm ".$fullpath."\n";
                $hash_del{$fullpath} = "rm -f ".$fullpath."\n";
            }  
        }
    }
}

sub eachFileSMG{
    my $filename = $_; 
    my $fullpath = $File::Find::name;
    #remember that File::Find changes your CWD, 
    #so you can call open with just $_
    
    
    if (-e $filename) { 
        #if(substr($filename, -(length($extension)+1)) eq $extension ){
        if($filename =~ m/Shotgun*\.ini$/ || $filename =~ m/metagenomics\.ini$/){
            my($vol, $dir, $file) = File::Spec->splitpath($fullpath);
            print STDERR $fullpath."\n";
            print STDERR $filename."\n";
            print STDERR $dir."\n";

            #if(-d "$dir/amplicontagger"){
            #    $hash{$dir} = "amplicontagger"
            #}elsif(-d "$dir/rrnatagger"){
            #    $hash{$dir} = "rrnatagger"
            #}
            $hash{$dir} = $dir 
        }
    }
}

print STDOUT "#!/bin/bash\n";
print STDOUT "#SBATCH --time=10:00:00\n";
print STDOUT "#SBATCH --nodes=1\n";
print STDOUT "#SBATCH --account=rrg-cgreer\n";
print STDOUT "#SBATCH -n 1\n";
print STDOUT "#SBATCH --mem=6000\n";
print STDOUT "#SBATCH -o ./stdout_archive.txt\n";
print STDOUT "#SBATCH -e ./stderr_archive.txt\n\n";

if($type eq "shotgun"){
    find(\&eachFileSMG, $indir);
    print STDERR Dumper(\%hash);
    my $last_indir_term = (split '\/', $indir)[-1];

    for my $key (keys %hash){
        my $prefix = $hash{$key};
        my $curr_indir = $key;
        my ($curr_indir2) = $curr_indir =~ m/($last_indir_term.*)$/;
        my ($curr_indir3) = $curr_indir =~ m/$last_indir_term\/(.*)$/;

        print STDOUT "#### Start $prefix #####\n";
        my $cmd = "cd $cwd && mkdir -p $curr_indir2 && cd $curr_indir2";
        print STDOUT $cmd."\n";

        print STDOUT "tar -chvf raw_reads.tar $curr_indir/raw_reads\n";
        print STDOUT "tar -chvf qced_reads.tar $curr_indir/qced_reads\n";
        print STDOUT "tar -zchvf assembly.tar.gz $curr_indir/assembly\n";
        if(-d "$curr_indir/gene_annotation"){
            print STDOUT "tar -chvf gene_annotation.tar $curr_indir/gene_annotation\n";
        }elsif(-d "$curr_indir/annotations"){
            print STDOUT "tar -chvf annotations.tar $curr_indir/annotations\n";
        }

        if(-d "$indir/misc"){
            print STDOUT "tar -zchvf misc.tar.gz $curr_indir/misc\n";
        }
        print STDOUT "tar -zchvf gene_prediction.tar.gz $curr_indir/gene_prediction\n";
        print STDOUT "tar -zchvf gene_abundance.tar.gz $curr_indir/gene_abundance\n";
        if(-d "$curr_indir/contigs_abundance"){
            print STDOUT "tar -zchvf contigs_abundance.tar.gz $curr_indir/contigs_abundance\n";
        }elsif(-d "$curr_indir/contig_abundance"){
            print STDOUT "tar -zchvf contig_abundance.tar.gz $curr_indir/contig_abundance\n";
        }
        #print STDOUT "tar -zchvf config_files.tar.gz $curr_indir/*.sh $curr_indir/*.tsv $curr_indir/*.ini $curr_indir/*.txt\n";
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name 'mapping_file*' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name 'workflow*' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name '*.ini' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";       
        print STDOUT "#### End $prefix #####\n\n\n";
    }
}elsif($type eq "amplicons"){

    # find directories were there is an .ini file + a rrnatagger or amplicontagger directory.
    find(\&eachFile, $indir);
    print STDERR Dumper(\%hash);
    #exit;
    my $last_indir_term = (split '\/', $indir)[-1];


    #my $prefix = ""; 
    #if(-d "$indir/rrnatagger"){
    #    $prefix = "rrnatagger";
    #}elsif(-d "$indir/amplicontagger"){
    #    $prefix = "amplicontagger";
    #}else{
    #    die "Can't find rrnatagger or amplicontagger dir\n";
    #}
    for my $key (keys %hash){
        my $prefix = $hash{$key};
        my $curr_indir = $key;
        my ($curr_indir2) = $curr_indir =~ m/($last_indir_term.*)$/;
        my ($curr_indir3) = $curr_indir =~ m/$last_indir_term\/(.*)$/;

        print STDOUT "#### Start $prefix #####\n";
        for my $key (keys %hash_del){
            print STDOUT $hash_del{$key};
        }

        my $cmd = "cd $cwd && mkdir -p $curr_indir2 && cd $curr_indir2";
        print STDOUT $cmd."\n";

        
        $cmd  = "if [ -d '$indir/$curr_indir3/amplicontagger/reads_12/feature_tables/rarefactions' ]; then\n";
        $cmd .= "    echo '$indir/$curr_indir3/amplicontagger/reads_12/feature_tables/rarefactions does exist.'\n";
        $cmd .= "    rm -rf $indir/$curr_indir3/amplicontagger/reads_12/feature_tables/rarefactions\n";
        $cmd .= "fi\n";
        print STDOUT $cmd."\n";
        
        $cmd  = "if [ -d '$indir/$curr_indir3/amplicontagger/reads_1/feature_tables/rarefactions' ]; then\n";
        $cmd .= "    echo '$indir/$curr_indir3/amplicontagger/reads_1/feature_tables/rarefactions does exist.'\n";
        $cmd .= "    rm -rf $indir/$curr_indir3/amplicontagger/reads_1/feature_tables/rarefactions\n";
        $cmd .= "fi\n";
        print STDOUT $cmd."\n";
        
        $cmd  = "if [ -d '$indir/$curr_indir3/amplicontagger/reads_12/otu_tables/rarefactions' ]; then\n";
        $cmd .= "    echo '$indir/$curr_indir3/amplicontagger/reads_12/otu_tables/rarefactions does exist.'\n";
        $cmd .= "    rm -rf $indir/$curr_indir3/amplicontagger/reads_12/otu_tables/rarefactions\n";
        $cmd .= "fi\n";
        print STDOUT $cmd."\n";
      
        $cmd  = "if [ -d '$indir/$curr_indir3/raw_reads' ]; then\n";
        $cmd .= "    echo '$indir/$curr_indir3/raw_reads does exist.'\n";
        $cmd .= "    tar -chvf raw_reads.tar --exclude='*preprocessed.fastq.gz' -C $indir/$curr_indir3/ ./raw_reads\n";
        $cmd .= "fi\n";
        print STDOUT $cmd."\n";

        $cmd = "find $indir/$curr_indir3 -name '*_complete.tar.gz' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "find $indir/$curr_indir3 -name '*export*.tar.gz' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name 'mapping_file*' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name 'barcodes*' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name '*.ini' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";       
        
        $cmd = "find $indir/$curr_indir3 -maxdepth 1 -name 'workflow*' -exec cp {} ./ \\;";
        print STDOUT $cmd."\n";
        
        $cmd = "tar -zchvf $prefix.tar.gz";
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
        $cmd .= " --exclude='*ncontam_nphix_1_trimmed_QCpassed.fastq'";
        $cmd .= " --exclude='*ncontam_nphix_1_trimmednoN.fastq'";
        $cmd .= " --exclude='*ncontam_nphix_2_trimmed_QCpassed.fastq'";
        $cmd .= " --exclude='*ncontam_nphix_2_trimmed.fastq'";
        $cmd .= " --exclude='*ncontam_nphix_1_trimmed.fastq'";
        $cmd .= " --exclude='*ncontam_nphix_2_trimmednoN.fastq'";
        $cmd .= " --exclude='*obs/demultiplexed'";
        $cmd .= " --exclude='*obs/dada2/filtered'";
        $cmd .= " --exclude='*_complete.tar.gz'";
        $cmd .= " --exclude='*/plots/tmp*'";
        $cmd .= " -C $curr_indir/ ./$prefix";
        print STDOUT $cmd."\n";
        print STDOUT "#### End $prefix #####\n\n\n";
    } 
}
print STDOUT "cd $cwd";

exit;
