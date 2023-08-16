#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use Cwd 'abs_path';

my $usage=<<'ENDHERE';
NAME:
preprocessIonTorrent.pl

PURPOSE:
To remove barcodes sequences in 5' region of reads and
place barcode sequence in reads's headers + cut reads
to a specific length.

INPUT:
--fastq <string>     : Fastq file.
OR --indir <string   : Indirectory where is located your fastq.

--barcodes <string>  : barcodes in fasta format.
--minLength <int>    : Reads shorter than this value will
                       be discarded.
--5primeSearch <int> : Length of nucl to search for primer
                       in the 5' region.
--cut                : cut all reads to --minLength <int>

OUTPUT:
STDOUT               : trimmed sequences
--log <string>       : log/stats
--outFailed <string> : Reads that did not either made the
                       length cutoff or did not have any
                       barcode sequence associated with it.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $fastq, $indir, $barcodes, $failed, $log, $minLength, $_5primeSearch, $cut);
my $verbose = 0;

GetOptions(
   'fastq=s'         => \$fastq,
   'indir=s'         => \$indir,
   'barcodes=s'      => \$barcodes,
   'failed=s'        => \$failed,
   'log=s'           => \$log,
   'minLength=s'     => \$minLength,
   '5primeSearch=i'  => \$_5primeSearch,
   'cut'             => \$cut,
   'verbose'         => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
if(!$fastq && !$indir){ die "provide at least --fastq OR --indir";}
die ("--barcodes missing\n") unless $barcodes;
die ("--failed missing\n") unless $failed;
die ("--log missing\n") unless $log;
#die ("--minLength missing\n") unless $minLength;
die ("--5primeSearch missing\n") unless $_5primeSearch;

## SUB
sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) {  
    if($filename !~ m/preprocess/){
      if(substr($filename, -9) eq ".fastq.gz"){
         $fastq = $fullpath;
       }elsif(substr($filename, -6) eq ".fastq"){
         $fastq = $fullpath;
         print STDERR $fastq."\n";
       }
     }
   }
}

my %hash;

open(FAILED, ">$failed") or die "Can't open $failed\n";
open(LOG, ">$log") or die "Can't open $log\n";

# Parse barcodes
my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
   my $seq = $curr->seq();
   $hash{$seq} = 0;
   #print STDERR $seq."\n";
}

# Parse fastq if indir, take fastq in indirectory (limitation: must be only one fastq in indir...)
if($indir){
   $indir = abs_path($indir);
   find (\&eachFile, $indir);
}
print STDERR "Preprocessing: ".$fastq."\n";

if(substr($fastq, -10) =~ m/\.fastq\.gz/){
   my $fastqgz = $fastq;
   $fastq =~ s/\.gz//;
   system("gunzip -c $fastqgz > $fastq");
   $? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly uncompressed ".$fastqgz." into .fastq file...\n";
}

if(!$minLength){
	## evaluate what length cutoff to use. 150, 175, 200 or 220
	my $length_100 = 0;
	my $length_125 = 0;
	my $length_150 = 0;
	my $length_175 = 0;
	my $length_200 = 0;
	my $length_225 = 0;
	my $length_250 = 0;
	my $length_300 = 0;
	my $length_350 = 0;
	my $length_400 = 0;
	my $total_seq = 0;
	
	my $ref_fastq_db = Iterator::FastqDb->new($fastq) or die("Unable to open Fastq file, $fastq\n");
	while(my $curr = $ref_fastq_db->next_seq() ){
	   $length_100++ if(length($curr->seq) >= 100); 
	   $length_125++ if(length($curr->seq) >= 125); 
	   $length_150++ if(length($curr->seq) >= 150); 
	   $length_175++ if(length($curr->seq) >= 175); 
	   $length_200++ if(length($curr->seq) >= 200); 
	   $length_225++ if(length($curr->seq) >= 225); 
	   $length_250++ if(length($curr->seq) >= 250); 
	   $length_300++ if(length($curr->seq) >= 300); 
	   $length_350++ if(length($curr->seq) >= 350); 
	   $length_400++ if(length($curr->seq) >= 400); 
	   $total_seq++;
	}
	
	$minLength = 100 if(($length_100 / $total_seq)*100 >= 60);
	$minLength = 125 if(($length_125 / $total_seq)*100 >= 60);
	$minLength = 150 if(($length_150 / $total_seq)*100 >= 60);
	$minLength = 175 if(($length_175 / $total_seq)*100 >= 60);
	$minLength = 200 if(($length_200 / $total_seq)*100 >= 60);
	$minLength = 225 if(($length_225 / $total_seq)*100 >= 60);
	$minLength = 250 if(($length_250 / $total_seq)*100 >= 60);
	$minLength = 300 if(($length_300 / $total_seq)*100 >= 60);
	$minLength = 350 if(($length_350 / $total_seq)*100 >= 60);
	$minLength = 400 if(($length_400 / $total_seq)*100 >= 60);
	
	if(!defined($minLength)){
	   print STDERR "optimal cutoff length lower than 100 bases! Will use 75 bases as cutoff...\n";
	   $minLength = 75;
	}
   print STDERR "Optimal minLength cutoff was found to be $minLength\n";
}

print STDERR "Trimming reads at $minLength bases\n";

my $success = 0;
my $nosuccess = 0;
my $ref_fastq_db = Iterator::FastqDb->new($fastq) or die("Unable to open Fastq file, $fastq\n");
while(my $curr = $ref_fastq_db->next_seq() ){
   
   # Discard reads shorter than the barcode search space.
   if(length($curr->seq) < $_5primeSearch){
      print FAILED $curr->output;
      $nosuccess++;
      next;
   }

   my $found = 0;
   for my $key (keys %hash){
      #print STDERR "barcode ".$key."\n";
      my $region = substr($curr->seq, 0, $_5primeSearch);

      if($region =~ m/($key)/g){
         my $pos = pos($region);
         
         if(length($curr->seq) - $pos < $minLength){
            print FAILED $curr->output;
            $nosuccess++;
         }else{
            $hash{$key}++;
            if($cut){
                print STDOUT $curr->header."#$key\n".substr($curr->seq, $pos, $minLength)."\n+\n".substr($curr->qual, $pos, $minLength)."\n";
            }else{
                print STDOUT $curr->header."#$key\n".substr($curr->seq, $pos)."\n+\n".substr($curr->qual, $pos)."\n";
            }
         }
         $success++;
         $found = 1;
      }
      last if($found == 1);
   }
   if($found == 0){
      print FAILED $curr->output;
      $nosuccess++;
   }
}
close(FAILED);

for my $key (keys %hash){
    print LOG "$key\t".commify($hash{$key})."\n";
}
print LOG "Preprocessing passed reads\t$success\n";
print LOG "Preprocessing failed reads\t$nosuccess\n";
close(LOG);

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

exit;
