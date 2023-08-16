#!/usr/bin/env perl

use strict;
no strict 'refs';
use warnings;

use Env qw/TMPDIR/;
use File::Which;
use String::Approx 'aslice';
use List::Util qw(sum);
use threads;
use threads::shared;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use Data::Dumper;
use FileCache;
use File::Basename;
use File::Temp;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
barcodes.pl

PURPOSE:
Split/demultiplex a fastq file by its barcode( i.e. index) sequences.

INPUT:
--infile <string>   : Sequence file.
--num_threads <int> : Number of threads.
--barcodes <string> : barcodes in fasta format.
--trim_length <int> : Reads will be trimmed to that length. Shorter
                      reads will be discarded. Optional. Only works
                      with the option --outdir (not --outfile).

OUTPUT:
--log <string>      : Log file having barcodes count.
--outfile <string>  : Sequence file in fastq having all good barcodes.
OR
--outdir <string>   : DIR where there will be one fastq file per barcode.
--suffix <string>   : Optional. Useful to specify if reads are R1 or R2.


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $log, $outdir, $barcodes, $num_threads, $trim_length, $suffix);

my $verbose = 0;

GetOptions(
    'infile=s'      => \$infile,
    'outfile=s'     => \$outfile,
    'outdir=s'      => \$outdir,
    'barcodes=s'    => \$barcodes,
    'log=s'         => \$log,
    'trim_length=i' => \$trim_length,
    'suffix=s'      => \$suffix,
    'num_threads=i' => \$num_threads,
    'verbose'       => \$verbose,
    'help'          => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--infile might be empty or wrong file path?\n") if((!-e $infile) and (!-s $infile));
die("--outfile OR --outdir arg required\n") unless($outfile or $outdir);
#die("--num_threads arg required\n") unless($num_threads);
$num_threads = $num_threads unless($num_threads);
die("Please specify a \$TMPDIR environment variable in your environment.\n") unless($TMPDIR);
die("Please specify a --log <string> arg/file.\n") unless($log);
my $phred = 33;

my $tmpdir = File::Temp->newdir(
    "tmpDirAmpliconTaggerBarcodesXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0
);

$suffix = "" unless($suffix);

my $outfile_failed;
if($outdir){
  $outfile_failed = $outdir."/UNKOWN.fastq";
  system("mkdir -p ".$outdir);
}elsif($outfile){
  my($filename, $directories) = fileparse($outfile);
  $outfile_failed = $directories."/UNKNOWN_barcodes.fastq";
}else{
  die "Something went wrong...\n";
}

## First, put barcodes into a hash.
my %hash_fh = ();
my %hash = ();
my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  my $header = $curr->header;
  $header =~ s/>//;
  $header =~ s/\//-/; # / will interfere in the moment of writing into file.
  $header = $header.$suffix;
  
  die "Barcode ".$curr->seq." is present in duplicate in barcode fasta file\n" if(exists $hash{$curr->seq});
  $hash{$curr->seq} = [$header, 0];
  
  # If outdir, open file handlers for writing and push them into hash structure.
  if($outdir){
    my $FASTQ_BARCODE = $outdir."/".$header.".fastq";
    push @{$hash{$curr->seq}}, $FASTQ_BARCODE;
    
    # Open it once to initialize the file and overwrite potential old file.
    open my $fh, '>', $FASTQ_BARCODE or die "Can't open ".$FASTQ_BARCODE." for writing.\n";
    close($fh);

    # Open all file handles at once
    open $fh, '>>', $hash{$curr->seq}->[2] or die $!;
    $hash_fh{$curr->seq} = $fh;
  }
}

#print Dumper(\%hash);
open my $FASTQ, '<', $infile or die $!;  
open my $FASTQ_OUT, '>', $outfile or die $! if($outfile);
#open my $FASTQ_OUT_FAILED, '>', $outfile_failed or die $!;

## If a no-zero start posns, find the start of the next full record.
my $size =  -s $infile or die "$! : $infile".". Cannot calculate size of file.\n";

## process records until the end of this threads section.
while( tell( $FASTQ ) < $size ) {
    my @lines = map scalar( <$FASTQ> ), 1 .. 4;
    chomp @lines;

    ## Validate header
    my $validate = new Iterator::ValidateFastq($lines[0], $lines[3], $phred);
    my $base = "@".$validate->base;
    my $barcode = $validate->barcode;
    my $pair = $validate->pair;
    my $header = $lines[0];
    my $seq = $lines[1];
    my $qual = $validate->qual;
    #my $qual = $lines[3];
    
    ## Perform barcode binning here.
    if(defined ($pair) and defined ($barcode)){
      if(exists $hash{$barcode}){
        $hash{$barcode}[1]++;
        if($outfile){
          print $FASTQ_OUT $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");    
        }elsif($outdir){
          #open my $fh, '>>', $hash{$barcode}->[2] or die $!;
          #print $fh $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");    
          #close($fh);
          print {$hash_fh{$barcode}} $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");    

            
        }else{
          die "Somehting went wrong...\n";
        }
      }else{
          #print $FASTQ_OUT_FAILED $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");    
      }
    
    }elsif(defined($pair)and !defined($barcode)){
      die "Fastq sequence file does not contain barcodes\n";
    
    }elsif(!defined($pair) and defined($barcode)){
        if(exists $hash{$barcode}){
            $hash{$barcode}[1]++;
            if($outfile){
                print $FASTQ_OUT $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");  
            }elsif($outdir){
                if($trim_length){
                    if(length($seq) >= $trim_length){
                        my $trimed_seq = substr $seq, 0, $trim_length;
                        my $trimed_qual = substr $qual, 0, $trim_length;
                        #open my $fh, '>>', $hash{$barcode}->[2] or die $!;
                        #print  $fh $base."#".$barcode."\n".$trimed_seq."\n+\n".$trimed_qual."\n" or die("NFS failure\n");  
                        #close($fh);
                        print {$hash_fh{$barcode}} $base."#".$barcode."\n".$trimed_seq."\n+\n".$trimed_qual."\n" or die("NFS failure\n");  
                    }#else do not print.
                }else{
                    #open my $fh, '>>', $hash{$barcode}->[2] or die $!;
                    #print  $fh $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");  
                    #close($fh);
                    print {$hash_fh{$barcode}} $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");  
            
                }
            }else{
                die "Something went wrong...\n";
            }
        }else{
            #print $FASTQ_OUT_FAILED $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");  
        }
    }else{
        die "Fastq sequence file does not contain barcodes and pair info\n";
    }   
}
close($FASTQ);
close($FASTQ_OUT) if($outfile);
#close($FASTQ_OUT_FAILED);

close $_ foreach values %hash_fh;

#print STDERR Dumper(\%hash);

#Print values in final log file.
my $total = 0;
while ( my ($key, $value) = each(%hash) ) {
    $total = $total + $hash{$key}[1];
}

open(LOG, ">".$log);
print LOG "#name\tsequence\tcount\tperc\n";
while ( my ($key, $value) = each(%hash) ) { 
    my $perc = sprintf '%.2f', $hash{$key}[1]/$total * 100;
    print LOG $hash{$key}[0]."\t".$key."\t".$hash{$key}[1]."\t".$perc."%\n";
}
close(LOG);

## REMOVE TEMP FILES
sub END{
  local $?;
  system("rm ".$tmpdir." -rf");
}

exit;
