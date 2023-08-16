#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use File::Temp;
use POSIX qw(mkfifo);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
dereplicateGz.pl

PURPOSE:

INPUT:
--fasta <string>    : Sequence file (fasta format).
OR
--fastq <string>    : Sequence file (fastq.gz format). Can be multiple, separated by a <,>
--sort              : Set if sequences are to be sorted and 
                      not dereplicated
--minsize <int>     : Only sequences >= minsize <int> will be
                      kept. So sequences < minsize are discarded.
                      For instance, if --minsize=2, only cluster 
                      that have a total of at least 2 will be kept.
        
OUTPUT:
STDOUT <string>     : Sequence file 100% clustered.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Microbiomes and Genomics
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $fasta, $fastq, $sort, $minsize, $minsize_per_sample);
my $verbose = 0;

GetOptions(
  'fasta=s'              => \$fasta,
  'fastq=s'              => \$fastq,
  'sort'                 => \$sort,
  'minsize=i'            => \$minsize,
  'minsize_per_sample=i' => \$minsize_per_sample,
  'verbose'              => \$verbose,
  'help'                 => \$help
);
if ($help) { print $usage; exit; }
 
## Validate
die "--fasta OR --fastq required.\n" if(!defined($fastq) && !defined($fasta));

my $tmpdir = File::Temp->newdir(
    "tmpDirDereplicateGzXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0
);


####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################

my %hash = ();
my %hash_barcodes = ();
my %hash_headers = ();

if($fastq){
    my $j = 0;
    my @fastqs = split(/,/, $fastq);
    for my $infile (@fastqs){
        my $basename = $infile;
        $basename =~ s{.*/}{};         # removes path  
        #$basename =~ s{\.[^.]+$}{};    # removes extension
        #$basename =~ s{\.[^.]+$}{};    # removes extension
        #TODO add more robust managing/check of file extension.
        $basename =~ s/_unmapped\.fastq\.gz//;
        $basename =~ s/_unmapped2b\.fastq\.gz//;
        $basename =~ s/_mapped\.fastq\.gz//;

        # make pipe
        my $pipe = "$tmpdir/reads.pipe";
        system("mkfifo $pipe");
        # gunzip to pipes
        print STDERR "[DEBUG] Integrating $infile\n";
        system("gunzip -c ".$infile." > $pipe &");
        
        my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fasta file, $pipe\n");
        while( my $curr = $ref_fastq_db->next_seq() ) {
            if(exists $hash{$curr->seq}){
                $hash{$curr->seq}++;
                if(exists $hash_barcodes{$curr->seq}{$basename} ){
                    $hash_barcodes{$curr->seq}{$basename}++;
                }else{
                    $hash_barcodes{$curr->seq}{$basename}=1;
                }
            }else{
                $hash{$curr->seq} = 1;
                $hash_barcodes{$curr->seq}{$basename}=1;
            }
      
            print STDERR "processing seq ".$j."\n" if((($j % 100000) == 0) && $verbose );
            $j++;
        }
        system("rm $pipe");
    }
}

#print STDERR Dumper(\%hash);
#print STDERR Dumper(\%hash_barcodes);

my $counter = 1;

$minsize = 0 unless($minsize);
#foreach my $key (sort {$hash{$b} <=> $hash{$a} } keys %hash){
#    my $line = ">".$hash_headers{$key}.";size=".$hash{$key}."\n".$key."\n";
#    $line =~ s/;;/;/;
#    $line =~ s/>+/>/;
#    print STDOUT $line if($hash{$key} >= $minsize);
#}
foreach my $key (sort {$hash{$b} <=> $hash{$a} } keys %hash){
    print STDOUT ">".$counter.";";
    for my $barcode (keys %{ $hash_barcodes{$key} }){  
        print STDOUT "#".$barcode."=".$hash_barcodes{$key}{$barcode}.";";
    }
    print STDOUT "size=".$hash{$key}."\n";
    print STDOUT $key."\n";
    $counter++;
}


## REMOVE TEMP FILES
sub END{
  local $?;
  system("rm ".$tmpdir." -rf");
}

exit;
