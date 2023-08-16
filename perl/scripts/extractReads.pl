#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::Pipe;
use Parallel::ForkManager;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use Data::Dumper;
use Cache::FastMmap;

my $usage=<<'ENDHERE';
NAME:
extractReads.pl

PURPOSE:

INPUT:
--infiles <string>    : bam files
--clusters <string>   : clusters
--flag <string>       : bam/sam file flag for instance: 0x040 or 0x080

OUTPUT:
--outdir <string>     : outdir

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infiles, $outdir, $clusters, $num_threads, $flag);
my $verbose = 0;

GetOptions(
   'infiles=s' 	   => \$infiles,
   'outdir=s' 	      => \$outdir,
   'clusters=s'      => \$clusters,
   'flag=s'          => \$flag,
   'num_threads=i'   => \$num_threads,
   'verbose' 	      => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

die "--infiles missing\n" unless($infiles);
die "--outdir missing\n" unless($outdir);
die "--clusters missing\n" unless($clusters);
die "--flag missing\n" unless($flag);
#die "--num_threads missing\n" unless($num_threads);

$num_threads = 1 unless($num_threads);

my $prefix = "";
if($flag eq "0x040"){
   $prefix = "R1";
}elsif($flag eq "0x080"){
   $prefix = "R2";
}

## MAIN
my @infiles = split(/,/, $infiles);

my %hash;
my %hash_fh1;
my %hash_unpaired;
open(IN, "<".$clusters) or die "Can't open $clusters\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $cluster = $row[0];
   my $ref_id = $row[1];
   
   #T17G_contig-13
   my $contigId;
   if($ref_id =~ m/^(\S+_contig-\d+)\..*/){
      $contigId = $1;
   }else{
      die "Did not matched regex...\n";
   }
   my $finalId = $cluster."=".$contigId;

   $hash{$finalId} = $cluster;

   #if(!exists $hash_fh1{$cluster}){
   #   my $output = $outdir."/".$cluster.".".$prefix.".fastq.gz";
   #   my $z = new IO::Compress::Gzip $output or die "gzip failed: $GzipError\n";
   #   $hash_fh1{$cluster} = $z;
   #}
   
   #if(!exists $hash_unpaired{$cluster}){
   #   my $output =  $outdir."/".$cluster.".".$prefix.".unpaired.fastq.gz";
   #   my $z = new IO::Compress::Gzip $output or die "gzip failed: $GzipError\n";
   #   $hash_unpaired{$cluster} = $z;
   #}
      
}
close(IN);

#print STDERR Dumper(\%hash);
#exit;

my $pm = new Parallel::ForkManager($num_threads);

my %catHash;
my %catHashUnpaired;

## Shared memory...
my $cachePaired = Cache::FastMmap->new();
my $cacheUnpaired = Cache::FastMmap->new();

my $i=0;
foreach my $bam (@infiles){
   $i++;
   my $pid = $pm->start and next;

   print STDERR "Processing $bam with loop number: $i\n";

   # open file handles for paired and unpaired data.
   
   # reads 1
   my $pipe = IO::Pipe->new();
   $pipe->reader("samtools view -f $flag $bam");
   while(<$pipe>){
      my @row  = split(/\t/, $_);
      my $read_id = $row[0];
      my $ref_id = $row[2];
      my $seq  = $row[9];
      my $qual = $row[10];
      my $currFlag = $row[1]; # if unpaired reads. Flag will be either: 4, 16 or 0.
      
      if(exists $hash{$ref_id}){
         
          my $curr_ref = $hash{$ref_id};
         
          if(!exists $hash_fh1{$curr_ref.$prefix.$i}){
            my $output = $outdir."/".$curr_ref.".".$prefix."_".$i."_.fastq.gz";
            my $z = new IO::Compress::Gzip $output or die "gzip failed: $GzipError\n";
            #my $z = 1;
            $hash_fh1{$curr_ref.$prefix.$i} = $z;
		      my $file_string = $cachePaired->get($curr_ref.".".$prefix);
		      $cachePaired->set($curr_ref.".".$prefix, $file_string .= $output.";");
            #$catHash{$curr_ref.".".$prefix}{$i} = $output;
         }
   
         if(!exists $hash_unpaired{$curr_ref.$prefix.$i}){
            my $output =  $outdir."/".$curr_ref.".".$prefix.".unpaired_".$i."_.fastq.gz";
            my $z = new IO::Compress::Gzip $output or die "gzip failed: $GzipError\n";
            #my $z=1;
            $hash_unpaired{$curr_ref.$prefix.$i} = $z;
		      my $file_string = $cacheUnpaired->get($curr_ref.".".$prefix);
		      $cacheUnpaired->set($curr_ref.".".$prefix, $file_string .= $output.";");
            #$catHashUnpaired{$curr_ref.".".$prefix}{$i} = $output;
         }


         # Then based on the flag arg. decide if reads 1 or 2 is being printed. This is for the /1 or /2 field in fastq header.
         if($flag eq "0x040"){
            if($currFlag == 99 || $currFlag == 83) { # properly paired mapped for reads 1
               #print {$hash_unpaired{$curr_ref}} "@".$read_id."/1\n".$seq."\n+\n".$qual."\n"; 
               $hash_fh1{$curr_ref.$prefix.$i}->print("@".$read_id."/1\n".$seq."\n+\n".$qual."\n");
            }else{
               #print {$hash_fh1{$curr_ref}} "@".$read_id."/1\n".$seq."\n+\n".$qual."\n"; 
               $hash_unpaired{$curr_ref.$prefix.$i}->print("@".$read_id."/1\n".$seq."\n+\n".$qual."\n");
            }
        }elsif($flag eq "0x080"){
            if($currFlag == 147 || $currFlag == 163) { #properly paired mapped for reads 2
               #print {$hash_unpaired{$curr_ref}} "@".$read_id."/1\n".$seq."\n+\n".$qual."\n"; 
               $hash_fh1{$curr_ref.$prefix.$i}->print("@".$read_id."/2\n".$seq."\n+\n".$qual."\n");
            }else{
               #print {$hash_fh1{$curr_ref}} "@".$read_id."/1\n".$seq."\n+\n".$qual."\n"; 
               $hash_unpaired{$curr_ref.$prefix.$i}->print("@".$read_id."/2\n".$seq."\n+\n".$qual."\n");
            }
         }
      }
   }
   $pm->finish;
}
$pm->wait_all_children;


## PAIRED READS
$pm = new Parallel::ForkManager($num_threads);

my @array = $cachePaired->get_keys(0);
@array = sort(@array);
foreach(@array){
   my $pid = $pm->start and next;
   
   my $cmd = "cat ";
   my $el = $cachePaired->get($_);
   my @el = split(/;/, $el);
   my $prefix = $el[0];
   $prefix =~ m/(canopy\d{5}\.\S+)_\d+_\.fastq\.gz/;
   $prefix = $1;
   @el = sort(@el);
   $cmd .= join(" ", @el);
   $cmd .= " > ".$outdir."/".$prefix.".fastq.gz";
   system($cmd);

   # Then cleanup
   foreach my $currFile (@el){
      system("rm $currFile");
   }
   $pm->finish;

}
$pm->wait_all_children;

## UNPAIRED READS
$pm = new Parallel::ForkManager($num_threads);
@array = $cacheUnpaired->get_keys(0);
@array = sort(@array);

foreach(@array){
   my $pid = $pm->start and next;
   
   my $cmd = "cat ";
   my $el = $cacheUnpaired->get($_);
   my @el = split(/;/, $el);
   my $prefix = $el[0];
   $prefix =~ m/(canopy\d{5}\.\S+\.unpaired)_\d+_\.fastq\.gz/;
   $prefix = $1;
   @el = sort(@el);
   $cmd .= join(" ", @el);
   $cmd .= " > ".$outdir."/".$prefix.".fastq.gz";
   system($cmd);

   # Then cleanup
   foreach my $currFile (@el){
      system("rm $currFile");
   }
   $pm->finish;

}
$pm->wait_all_children;

exit;
