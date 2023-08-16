#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use Bio::Seq;
use Bio::Tools::IUPAC;
use String::Approx 'aslice';
use List::Util qw(sum);
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
tagsQC.pl

PURPOSE:
Filter amplicon reads.

INPUT:
--infile <string>          : Sequence file. Assumes it is fwd reads for single-end reads data or
                             merged/assembled amplicon fastq.
--infile_R2 <string>       : Optional. Rev reads.  if the --infile_R2 <string> arg. is provided, 
                             assumed --infile <string> is for fwd reads and only output reads that 
                             satisfy criteria for both ends (R1 and R2). 
--primer_5_prime <string>  : Fasta file containing ONE primer sequence
--primer_3_prime <string>  : Fasta file containing ONE primer sequence
--length_5_prime <int>     : Length of the 5' region to search for primer
--length_3_prime <int>     : Length of the 3' region to search for primer 
--primer_mismatch <int>    : % of mismatch allowed when searching for primer.
                             Default: 15.
--cut_first <int>          : Cut the first <int> bases
--cut_last <int>           : Cut last <int> bases
--qscore_1 <int>           : Average Qscore.
--N <int>                  : Max number of N tolerated
--qscore_2 <int>           : Low quality threshold
--lq_threshold <int>       : Number of tolerated bases below low quality 
                             threshold. Sequences having more than <lq_threshold>
                             bases below <qscore_2> will be rejected.
--max_length <int>         : Sequences having more than <int> bases length
                             will be rejected.
--min_length <int>         : Sequences having less than <int> bases length
                             will be rejected. 
--qual <int>               : Phred scale. Either 33 or 64.
--noqc                     : Set if no QC is to be done. For instance if primers
                             are to be removed, but no QC to be done.
--remove_1st_N             : Trim first base of read in 5' region if it is a N
                             (Optional).
--num_threads <int>        : Number of threads to use.
--reject_unmatched         : Set this flag if sequences not matching
                             at least one primer seq. To use in combination with
                             either --primer_5_prime or --primer_3_prime
--correct_orientation      : PacBio data reads can sometimes be in reverse orientation 
                             or reverse complement orientation.
                             If that's the case, they should be re-orientated in the 
                             fwd --> rev orientation. Option only available in combination
                             with --primer_5_prime and --primer_3_prime options.

OUTPUT:
--outfile_failed <string>     : File containing reads that failed to pass
                                QC.
--outfile_failed_R2 <string>  : File containing reads that failed to pass
                                QC. Needed if --infile_R2 <string> is specified.

OPTIONAL:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $debug, $infile, $outfile, $outfile_failed, $num_threads,
    $primer_5_prime, $primer_3_prime, $length_5_prime, $length_3_prime,
    $primer_mismatch, $cut_first, $cut_last, $qscore_1, $N, $lq_threshold,
    $qscore_2, $max_length, $min_length, $qual, $noqc, $reject_unmatched,
    $remove_N, $infile_R2, $outfile_R2, $outfile_failed_R2, $correct_orientation, $log
);
my $verbose = 0;

GetOptions(
  'log=s'                      => \$log,
  'infile=s'                   => \$infile,
  'outfile=s'                  => \$outfile,
  'infile_R2=s'                => \$infile_R2,
  'outfile_R2=s'               => \$outfile_R2,
  'num_threads=i'              => \$num_threads,
  'infile=s'                   => \$infile,
  'primer_5_prime=s'           => \$primer_5_prime,
  'primer_3_prime=s'           => \$primer_3_prime,
  'length_5_prime=i'           => \$length_5_prime,
  'length_3_prime=i'           => \$length_3_prime,
  'primer_mismatch=i'          => \$primer_mismatch,
  'remove_1st_N'               => \$remove_N,
  'cut_last=i'                 => \$cut_last,
  'cut_first=i'                => \$cut_first,
  'qscore_1=i'                 => \$qscore_1,
  'N=i'                        => \$N,
  'lq_threshold=i'             => \$lq_threshold,
  'qscore_2=i'                 => \$qscore_2,
  'max_length=i'               => \$max_length,
  'min_length=i'               => \$min_length,
  'qual=i'                     => \$qual,
  'outfile_failed=s'           => \$outfile_failed,
  'outfile_failed_R2=s'        => \$outfile_failed_R2,
  'reject_unmatched'           => \$reject_unmatched,
  'noqc'                       => \$noqc,
  'correct_orientation'        => \$correct_orientation,
  'debug'                      => \$debug,
  'verbose'                    => \$verbose,
  'help'                       => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--infile $infile might be empty or wrong file path?\n") if((!-e $infile) and (!-s $infile));
die("--outfile arg required\n") unless($outfile);
die("--num_threads arg required\n") unless($num_threads);
die("If --infile_R2 is provided, --outfile_R2 and --outfile_failed_R2 also has to provided\n") if($infile_R2 && !$outfile_R2);
die("If --infile_R2 is provided, --outfile_R2 and --outfile_failed_R2 also has to provided\n") if(!$infile_R2 && $outfile_R2);
if($primer_mismatch < 0 or $primer_mismatch > 51){
    die("--primer_mismatch int value between 0 and 50 required\n");
}
if($max_length){
    die("--max_length Please enter a value <= 5000\n") if($max_length > 5000);
}
die("--outfile outfile required\n") unless $outfile;

if( ($primer_5_prime and !$length_5_prime) or (!$primer_5_prime and $length_5_prime) ){
    die "--primer_5_prime AND --length_5_prime must be provided\n"
}
if( ($primer_3_prime and !$length_3_prime) or (!$primer_3_prime and $length_3_prime) ){
    die "--primer_3_prime AND --length_3_prime must be provided\n"
}
die("--log <string> required\n") unless $log;
my $left_primer = 0;
my $right_primer = 0;
if($primer_mismatch > -1 && $primer_5_prime && $length_5_prime){
    $left_primer = 1;
}
if($primer_mismatch > -1 && $primer_3_prime && $length_3_prime){
    $right_primer = 1;
}
if($noqc){
    $noqc = 1;
}else{
    $noqc = 0;
}
if($noqc == 0){
    $qscore_1 = 30 unless $qscore_1;
    $N = 3 unless $N;
    $lq_threshold = 3 unless $lq_threshold;
    $qscore_2 = 15 unless $qscore_2;
    #die "please enter either 33 or 64 for --qual arg\n" unless $qual == 33 or $qual == 64;
}

if(!defined($reject_unmatched)){ $reject_unmatched = 0; }else{$reject_unmatched = 1;}
print STDERR "reject_unmatched value: ".$reject_unmatched."\n" if($debug);
$cut_first = 0 unless($cut_first);
$cut_last = 0 unless($cut_last);

## Initialize hash for statistics
my %hash_stats;

## PREPARE PRIMERS
my @left_primer;
my @right_primer;

## Initialize log file.
open(LOG, ">".$log) or die "Can't open $log\n";

## GENERATE PRIMER SEQS
if($primer_5_prime){
    my $left_primer_seq;
    my $ref_fasta_db = new Iterator::FastaDb($primer_5_prime) or die("Unable to open Fasta file, $primer_5_prime\n");
    my $i=0;
    while( my $ref_seq = $ref_fasta_db->next_seq() ) {
        $left_primer_seq = uc($ref_seq->seq());
        $i++;
    }
    die "Please provide only one sequence in fasta file\n" if($i > 1);
    @left_primer = generate_primer_combination($left_primer_seq);
}

if($infile_R2){
   if($primer_3_prime){
       my $right_primer_seq;
       my $ref_fasta_db = new Iterator::FastaDb($primer_3_prime) or die("Unable to open Fasta file, $primer_3_prime\n");
       my $i=0;
       while( my $ref_seq = $ref_fasta_db->next_seq() ) {
           $right_primer_seq = uc($ref_seq->seq()); # do not rc if scanning from 5'
           #$right_primer_seq =~ tr/ACGT/TGCA/; #Just convert to complement since the reads will be reversed
           $i++;
       }
       die "Please provide only one sequence in fasta file\n" if($i > 1);
       @right_primer = generate_primer_combination($right_primer_seq);
   }
}else{
   if($primer_3_prime){
       my $right_primer_seq;
       my $ref_fasta_db = new Iterator::FastaDb($primer_3_prime) or die("Unable to open Fasta file, $primer_3_prime\n");
       my $i=0;
       while( my $ref_seq = $ref_fasta_db->next_seq() ) {
           $right_primer_seq = uc($ref_seq->seq());
           $right_primer_seq =~ tr/ACGT/TGCA/; #Just convert to complement since the reads will be reversed
           $i++;
       }
       die "Please provide only one sequence in fasta file\n" if($i > 1);
       @right_primer = generate_primer_combination($right_primer_seq);
   }
}

####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################

print $outfile."\n" if($verbose);
open my $FASTQ, '<', $infile or die $!;
open my $FASTQ_OUT, '>', $outfile or die $!;
open my $FASTQ_OUT_FAILED, '>', $outfile_failed or die $!;
my $FASTQ2;
my $FASTQ2_OUT;
my $FASTQ2_OUT_FAILED;

## process records until the end of this threads section.
if($infile_R2){
   open $FASTQ2, '<', $infile_R2 or die $!;
   open $FASTQ2_OUT, '>', $outfile_R2 or die $!;
   open $FASTQ2_OUT_FAILED, '>', $outfile_failed_R2 or die $!;
   
   my $fastq_db = new Iterator::FastqDb($infile) or die("Unable to open Fasta file, $infile\n");
   my $fastq_db2 = new Iterator::FastqDb($infile_R2) or die("Unable to open Fasta file, $infile_R2\n");
   while(my $seq = $fastq_db->next_seq() ) {
      $hash_stats{Input_R1}++;
      my $seq2 = $fastq_db2->next_seq();
      $hash_stats{Input_R2}++;

      my $base = $seq->base;  
      my $barcode = $seq->barcode;  
      my $pair = $seq->pair;
      my $header = $seq->header;
      my $qual = $seq->qual;
      $seq = $seq->seq;
      
      my $base2 = $seq2->base;  
      my $barcode2 = $seq2->barcode;  
      my $pair2 = $seq2->pair;
      my $header2 = $seq2->header;
      my $qual2 = $seq2->qual;
      $seq2 = $seq2->seq;
  
      # Sanity check before filtering:
      if($base ne $base2){
          print STDERR "[DEBUG] base R1: $base\n";
          print STDERR "[DEBUG] base R2: $base2\n";
          die "Fastq entries do not match between pairs (before filtering)\n";
      }

      my ($R1_p, $R1_f, $filtered_base1) = filter($base, $header, $seq,  $qual,  $FASTQ_OUT,  $FASTQ_OUT_FAILED, "fwd", "do_left_only");
      my ($R2_p, $R2_f, $filtered_base2) = filter($base2, $header2,$seq2, $qual2, $FASTQ2_OUT, $FASTQ2_OUT_FAILED, "rev", "do_right_only");
      if($R1_p && $R2_p){
         print $FASTQ_OUT $R1_p;
         print $FASTQ2_OUT $R2_p;
         $hash_stats{Passed_pair}++;
      }else{
         $hash_stats{Failed_pair}++;
         if($R1_f){
            print $FASTQ_OUT_FAILED $R1_f;
            $hash_stats{failed_R1}++;
	     }
         if($R2_f){
            print $FASTQ2_OUT_FAILED $R2_f;
            $hash_stats{failed_R2}++;
         }
      }
      # Sanity check after filtering:
      if($filtered_base1 ne $filtered_base2){
          print STDERR "[DEBUG] base R1: $filtered_base1\n";
          print STDERR "[DEBUG] base R2: $filtered_base2\n";
          die "Fastq entries do not match between pairs (after filtering)\n";
      }
   }
}else{
   my $fastq_db = new Iterator::FastqDb($infile) or die("Unable to open Fasta file, $infile\n");
   while(my $seq = $fastq_db->next_seq() ) {
      $hash_stats{Input_R1}++;
      my $base = $seq->base;  
      my $barcode = $seq->barcode;  
      my $pair = $seq->pair;
      my $header = $seq->header;
      my $qual = $seq->qual;
      my $seq = $seq->seq;
   
      my ($R1_p, $R1_f, $filtered_base) = filter($base, $header,$seq, $qual, $FASTQ_OUT, $FASTQ_OUT_FAILED, "fwd", "do_both");
      if($R1_p){
         print $FASTQ_OUT $R1_p;
         $hash_stats{Passed_R1}++;
      }else{
         print $FASTQ_OUT_FAILED $R1_f;
         $hash_stats{Failed_R1}++;
      }
   }
}
close $FASTQ_OUT;
close $FASTQ_OUT_FAILED;

for my $key (sort keys %hash_stats){
	print LOG $key.": ".$hash_stats{$key}."\n";
}
close(LOG);

## FILTER SEQUENCE
## input : sequence_header_string, dna_string, quality_string, out_filehandler, out_failed_filehandler.
## output: null
sub filter{
   my($base, $header, $newseq, $newqual, $OUT, $OUT_F, $orientation, $what_to_do) = @_;
   my $rc_newseq = reverse_complement_IUPAC($newseq);
   my $r_newqual = reverse($newqual);
   my %qscore_hash;
   setQscoreHash(\%qscore_hash);
   my $passed_read_entry = 0; 
   my $failed_read_entry = 0;

   my $curr_length_5_prime;
   my $curr_length_3_prime;
   my @curr_left_primer;
   my @curr_right_primer;

   my $found_left_primer = 0; 
   my $found_right_primer = 0; 
   my $found_at_least_one_primer = 0;
   
   #REMOVE LEFT PRIMER
   if($orientation eq "rev"){
      #$curr_length_5_prime = $length_3_prime;
      #$curr_length_3_prime = $length_5_prime;
      #@curr_left_primer    = @right_primer;
      #@curr_right_primer   = @left_primer;
      $curr_length_5_prime = $length_5_prime;
      $curr_length_3_prime = $length_3_prime;
      @curr_right_primer   = @right_primer;
      @curr_left_primer    = @left_primer;

   }elsif($orientation eq "fwd"){
      $curr_length_5_prime = $length_5_prime;
      $curr_length_3_prime = $length_3_prime;
      @curr_left_primer    = @left_primer;
      @curr_right_primer   = @right_primer;
   }
   
   my $matched_reverse_orientation_left_primer = 0;
   if($left_primer == 1 && ($what_to_do eq "do_both" || $what_to_do eq "do_left_only")){
      my $string = substr($newseq, 0, $curr_length_5_prime);
      my $matched = 0;
      my $lowest_length = 10000; #Arbitrary unrealistic high value.
      #my $max_length = 0;

	  print STDERR "Attempting to match fwd primers...\n" if($debug);
      foreach(@curr_left_primer){
         print STDERR "Current FWD Primer:\t".$_."\n" if($debug);
         print STDERR "primer_mismatch:".$primer_mismatch."\n" if($debug);
         
         my @index = aslice($_, ["i ".$primer_mismatch."%"] , $string);
         if(defined $index[0][0]){
            my $length_from_0 = $index[0][0] + $index[0][1];
            print STDERR "Found fwd->rev Index:".$index[0][0]." + ".$index[0][1]."\n" if($debug);
            $lowest_length = $length_from_0 if($length_from_0 < $lowest_length);
            $matched = 1;
         }elsif($correct_orientation){ # If sequence is in reverse orientation (rev-->fwd), rc complement it and look again for primers
            undef(@index);
            my $rc_string = substr($rc_newseq, 0, $curr_length_5_prime);
            @index = aslice($_, ["i ".$primer_mismatch."%"] , $rc_string);
            if(defined $index[0][0]){
               my $length_from_0 = $index[0][0] + $index[0][1];
               print STDERR "Found rev->fwd Index:".$index[0][0]." + ".$index[0][1]."\n" if($debug);
               $lowest_length = $length_from_0 if($length_from_0 < $lowest_length);
               $matched_reverse_orientation_left_primer = 1;
            }
            # if no primer match, can mean that sequence is too short to match primer. Should try to match reverse primer
            # if reverse primer match, keep the sequence, because we corrected the orientation. If it fails, reject it...
         }
      }

	  if($matched == 1 and $matched_reverse_orientation_left_primer == 1){
	     print STDERR "Found fwd primer in both 5' and 3' region... Something is wrong...\n";

      }elsif($matched == 1 and $lowest_length < 10000){
         print STDERR "Matched FWD primer in forward orientation: ".substr($newseq, 0, $lowest_length)."\n" if($debug);
         $newseq = substr($newseq, $lowest_length);
         $newqual = substr($newqual,$lowest_length);
         $found_left_primer = 1;
         $hash_stats{Matched_left_primer}++ ;
      }elsif($matched_reverse_orientation_left_primer == 1 and $lowest_length < 10000){
         print STDERR "Matched FWD primer in reverse orientation: ".substr($rc_newseq, 0, $lowest_length)."\n" if($debug);
         $newseq = substr($rc_newseq, $lowest_length);
         $newqual = substr($r_newqual,$lowest_length);
         $found_left_primer = 1;
         $hash_stats{reverse_complement_reads_5p_primer}++;
         $hash_stats{Matched_left_primer}++ ;
      }else{
         print STDERR "Failed to match FWD primer\n" if($debug);
      }
   }

   #REMOVE RIGHT PRIMER
   if($right_primer == 1 and ($what_to_do eq "do_both" || $what_to_do eq "do_right_only")){
	  print STDERR "Attempting to match rev primers...\n" if($debug);
      # if no primer match, can mean that sequence is too short to match primer. Should try to match reverse primer
      # if reverse primer match, keep the sequence, because we corrected the orientation. If it fails, reject it...
      # Reverse both seq/qual
      # Note: if found left primer in reverse orientation, newseq and newqual are already reversed-complemented and ready to be processed 'normally' in this section...
      if($infile and $infile_R2 and $orientation eq "rev"){
          print STDERR "infile and infile R2 and orientation=rev\n" if($debug);
          $newseq = $newseq;
          $newqual = $newqual;
          print STDERR "newseq: ".$newseq."\n" if($debug); 
      }else{
          print STDERR "infile (R1) only\n" if($debug);
          $newseq = reverse($newseq);
          $newqual = reverse($newqual);
      }
      my $newseq_rc = reverse_complement_IUPAC($newseq);
      my $newqual_r = reverse($newqual);
      my $string = substr($newseq, 0, $curr_length_3_prime);
      my $string_rc = substr($newseq_rc, 0, $curr_length_3_prime);
      my $matched = 0;
      my $matched_cc = 0;
      my $matched_c = 0;
      my $lowest_length = 10000; #Arbitrary unrealistic high value.
      
      foreach(@curr_right_primer){
      
         print STDERR "REV Primer:\t".$_."\n" if($debug);
         print STDERR "string matching against:\t".$string."\n" if($debug);
         my @index = aslice($_, ["i ".$primer_mismatch."%"] , $string);
   
         if(defined $index[0][0]){
            my $length_from_0 = $index[0][0] + $index[0][1];
            $lowest_length = $length_from_0 if($length_from_0 < $lowest_length);
            $matched = 1;
            $matched_c = 1;
            $found_right_primer = 1;
            print STDERR "found REV Primer:\t".$_."\n" if($debug);
            $hash_stats{Matched_right_primer}++ ;

         }elsif($correct_orientation and $matched == 0){ # If sequence is in reverse orientation (rev-->fwd) and no right primer have been found, 
                                                         # rc complement the primer and look again... Since the primer is already complemented, just reverse it.
            print STDERR "string matching against rc:\t".$string_rc."\n" if($debug);
            undef(@index);
            #my $rc_string = reverse_complement_IUPAC($string);
            @index = aslice($_, ["i ".$primer_mismatch."%"] , $string_rc);
            if(defined $index[0][0]){
               my $length_from_0 = $index[0][0] + $index[0][1];
               print STDERR "Match rev(seq) against REV primer Index:".$index[0][0]." + ".$index[0][1]."\n" if($debug);
               $lowest_length = $length_from_0 if($length_from_0 < $lowest_length);
               #$max_length = $length_from_0 if($length_from_0 > $max_length);
               $matched = 1;
               $matched_cc = 1;
               $found_right_primer = 1;
            }
         }
      }
      print STDERR "found_left_primer=$found_left_primer   found_right_primer=$found_right_primer\n" if($debug);
      if($reject_unmatched == 1){ # if reject primers = true; means reject if at least one unmatched primer unmatched primers...
          if($what_to_do eq "do_both" || $what_to_do eq "do_left_only"){
             if(($found_left_primer == 1 and $found_right_primer == 1) and $lowest_length < 10000){
                print STDERR "Matched both FWD and REV primers: ".substr($newseq, 0, $lowest_length)."\n" if($debug);
                if($matched_cc){
                   $newseq = substr($newseq_rc, $lowest_length);
                   $newqual = substr($newqual_r ,$lowest_length);
                }else{
                   $newseq = substr($newseq, $lowest_length);
                   $newqual = substr($newqual,$lowest_length);
                }
                $found_at_least_one_primer = 1;
                $hash_stats{Matched_both_primers}++ ;

             }elsif(($found_left_primer == 0 and $found_right_primer == 1) and $lowest_length < 10000){
                print STDERR "Matched REV primer only: ".substr($newseq, 0, $lowest_length)."\n" if($debug);
                $found_at_least_one_primer = 1;
                if($what_to_do eq "both"){
                    $hash_stats{Matched_right_primer_only}++ ;
                }else{
                    $hash_stats{Matched_right_primer}++ ;
                }

             }elsif(($found_left_primer == 1 and $found_right_primer == 0) and $lowest_length < 10000){
                $found_at_least_one_primer = 1;
                $hash_stats{Matched_left_primer_only}++ ;
     
             }else{
                print STDERR "Failed to match both FWD and REV primers\n" if($debug);
                if($reject_unmatched){
                   $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
                   $hash_stats{Failed_reads_due_to_unmatched_3p_and_5p_primer}++ ;
                   $found_at_least_one_primer = 0;
                }
             }
             print STDERR "Right primer Header:\t".$header."\n" if $debug;
          
             $newseq = reverse($newseq);
             $newqual = reverse($newqual);
          }elsif($what_to_do eq "do_right_only"){
               if($matched == 1){
                   $newseq = substr($newseq, $lowest_length);
                   $newqual = substr($newqual,$lowest_length);
               }
          }
      }else{ # if reject_unmatched is false
         if($matched == 1 and $lowest_length < 10000){
            print STDERR "Matched REV primer: ".substr($newseq, 0, $lowest_length)."\n" if($debug);
            $newseq = substr($newseq, $lowest_length);
            $newqual = substr($newqual,$lowest_length);
            $hash_stats{complement_reads_3p_primer}++ if($matched_c);
            $hash_stats{reverse_complement_reads_3p_primer}++ if($matched_cc);
   
         }else{
            print STDERR "Failed to match REV primer\n" if($debug);
            $hash_stats{Read_unmatched_3p_primer_but_not_rejected}++;
         }
         print STDERR "Right primer header:\t".$header."\n" if $debug;
      
         $newseq = reverse($newseq);
         $newqual = reverse($newqual);
      } 
   }#end if right primer.
   

   #Cut sequences before going into quality threshold.
   if($cut_first or $cut_last){
      my $length = (length($newseq)) - $cut_last - $cut_first;
      die "First nucleotides to cut must be higher or equal than/to read length.\n" if($cut_first >= length($newseq));
      $newseq  = substr($newseq, $cut_first, $length);
      $newqual = substr($newqual, $cut_first, $length);
   }
  
   #Remove 1st base if it is a N (typically seen in 5'reads).
   if($remove_N){
      if(substr($newseq, 0, 1) eq "N"){
         $newseq = substr($newseq, 1);  
         $newqual = substr($newqual, 1);  
      }
    
      # Also remove last base if last base is a N.  
      if(substr($newseq, -1) eq "N"){
         $newseq = substr($newseq, 0, -1);  
         $newqual = substr($newqual, 0, -1);  
      }  
   }

   my $passed_max_length = 0;
   my $passed_min_length = 0;
   my $passed_qual = 0;
   #QUALITY CONTROL (REMOVE READS WITH TOO MUCH Ns AND TOO MUCH NT HAVING LOW QUAL SCORE)
   if($noqc == 0){
      my @qual = unpack("C*", $newqual);
      my $prob_sum = 0;
      foreach(@qual){

         # For PacBio CCS data it seems some reads can have really low quality values (below 33...)
         if($_ < $qual){
            $_ = $qual;
         }

         $prob_sum = $prob_sum + $qscore_hash{$_ - $qual};
         die "Q score value does not exists in reference hash...".$_ - $qual."\n" if(!exists $qscore_hash{$_ - $qual});
      }
    
      if(@qual <= 0){
         $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
         print STDERR "Sequence $header has length of 0 bases. It you supplied the cut_last or cut_last paramters, it is possible that it made the sequence too short. Try to cut shorter...\n" if($debug);
         #next;
      }else{ 

         my $average = $prob_sum/@qual;
         print STDERR "\nAverage:\t".$average."\n" if($debug);
       
         my $Q_count = 0;
         foreach(@qual){
            $Q_count++ if ($_ - $qual) < $qscore_2;
         }
         print STDERR "Q_count:\t".$Q_count."\n" if $debug;
   
         my @seq = split(//, $newseq);
         my $N_count = 0;
         foreach(@seq){
            $N_count++ if uc($_) eq "N";
         }
         #print $N_count."\n" if $debug;
   
         if($average > $qscore_hash{$qscore_1},
            and $N_count < $N,
            and $Q_count < $lq_threshold ){
            
            $passed_qual = 1; 
            print STDERR "Average Q score is:\t".$average."\n" if $debug;
            #Max length
            if(defined($max_length) && $max_length != 0){
               if( length($newseq) <= $max_length ){
                   print STDERR $header."max_length passed\n" if $debug;
                   $passed_max_length = 1;
               }else{
                  $hash_stats{Failed_reads_due_to_read_gt_max_length}++;
                  $passed_max_length = 0;
               }
            }else{
                $passed_max_length = 1;
            }
            
            #Min length
            if($min_length != 0){
               if( length($newseq) >= $min_length ){
                  $passed_min_length = 1;  
               }else{
                  $hash_stats{Failed_reads_due_to_read_lt_min_length}++;
                  $passed_min_length = 0;  
               }
            }else{
               $passed_min_length = 1;
            } 
            #Reject read if it does not encounter filtering parameters
         }else{
            $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
            $hash_stats{Failed_reads_due_to_read_of_low_qual}++;
         }
      }
   # If no QC is to be done.
   }else{
      $passed_qual = 1; 
      if($cut_first or $cut_last){
         my $length = (length($newseq)) - $cut_last - $cut_first;
         $passed_read_entry = $header."\n".substr($newseq, $cut_first, $length)."\n+\n".substr($newqual, $cut_first, $length)."\n";
      }else{
         print $OUT $header."\n".$newseq."\n+\n".$newqual."\n";
         $passed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
      }
   }
   #if($reject_unmatched == 1 and $found_at_least_one_primer == 1){
   if($reject_unmatched == 1){ #means if reject_unmatched is true.
       if($what_to_do eq "do_both"){
           if($found_left_primer == 0 or $found_right_primer == 0){ #reject if one primer pair is missing
              $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
              $passed_read_entry = 0;
           }elsif($found_left_primer == 1 and $found_right_primer == 1){ #accept only if both primer pairs were found.
              $failed_read_entry = 0;
              $passed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
           }
        }elsif($what_to_do eq "do_left_only"){
           if($found_left_primer == 0){
              $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
              $passed_read_entry = 0;
           }elsif($found_left_primer == 1){ #accept only if both primer pairs were found.
              $failed_read_entry = 0;
              $passed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
           }
       }elsif($what_to_do eq "do_right_only"){
           if($found_right_primer == 0){
              $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
              $passed_read_entry = 0;
           }elsif($found_right_primer == 1){
              $failed_read_entry = 0;
              $passed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
           }
       }
   }elsif($reject_unmatched == 0){
      if($passed_max_length == 1 and $passed_min_length == 1 and $passed_qual == 1){
         $passed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
         $failed_read_entry = 0;	
      }else{
         $passed_read_entry = 0;	
         $failed_read_entry = $header."\n".$newseq."\n+\n".$newqual."\n";
      }
   }
   return ($passed_read_entry, $failed_read_entry, $base);
}

## GENERATE ALL POSSIBLE PRIMER COMBINATION
## input: one dna sequence that can contain ambiguous bases.
## output: an array of dna sequences with no ambiguous bases.
sub generate_primer_combination{
    my $primer_seq = $_[0];
    my @array = ();

    #For some obscur reason, Bio::Seq would always accept sequence as provided...
    $primer_seq =~ m/(\S+)/;
    $primer_seq = $1;

    my $ambiseq = Bio::Seq->new(-seq => $primer_seq, -alphabet => 'dna');
    my $stream  = Bio::Tools::IUPAC->new(-seq => $ambiseq);

    while (my $uniqueseq = $stream->next_seq()) {
        push(@array, $uniqueseq->seq());
    }
    return @array;
}

## SET QSCORE HASH
## Will populate a hash ref passed in paramter.
## input: a hash reference
## output: null;
sub setQscoreHash{
    my($hash_ref) = @_;

  $hash_ref->{0}=1;
  $hash_ref->{1}=1;
  $hash_ref->{2}=2;
  $hash_ref->{3}=3;
  $hash_ref->{4}=4;
  $hash_ref->{5}=5;
  $hash_ref->{6}=6;
  $hash_ref->{7}=7;
  $hash_ref->{8}=8;
  $hash_ref->{9}=9;
  $hash_ref->{10}=10;
  $hash_ref->{11}=11;
  $hash_ref->{12}=12;
  $hash_ref->{13}=13;
  $hash_ref->{14}=14;
  $hash_ref->{15}=15;
  $hash_ref->{16}=16;
  $hash_ref->{17}=17;
  $hash_ref->{18}=18;
  $hash_ref->{19}=19;
  $hash_ref->{20}=20;
  $hash_ref->{21}=21;
  $hash_ref->{22}=22;
  $hash_ref->{23}=23;
  $hash_ref->{24}=24;
  $hash_ref->{25}=25;
  $hash_ref->{26}=26;
  $hash_ref->{27}=27;
  $hash_ref->{28}=28;
  $hash_ref->{29}=29;
  $hash_ref->{30}=30;
  $hash_ref->{31}=31;
  $hash_ref->{32}=32;
  $hash_ref->{33}=33;
  $hash_ref->{34}=34;
  $hash_ref->{35}=35;
  $hash_ref->{36}=36;
  $hash_ref->{37}=37;
  $hash_ref->{38}=38;
  $hash_ref->{39}=39;
  $hash_ref->{40}=40;
  $hash_ref->{41}=41;
  $hash_ref->{42}=42;
  $hash_ref->{43}=43;
  $hash_ref->{44}=44;
  $hash_ref->{45}=45;
  $hash_ref->{46}=46;
  $hash_ref->{47}=47;
  $hash_ref->{48}=48;
  $hash_ref->{49}=49;
  $hash_ref->{50}=50;
  $hash_ref->{51}=51;
  $hash_ref->{52}=52;
  $hash_ref->{53}=53;
  $hash_ref->{54}=54;
  $hash_ref->{55}=55;
  $hash_ref->{56}=56;
  $hash_ref->{57}=57;
  $hash_ref->{58}=58;
  $hash_ref->{59}=59;
  $hash_ref->{60}=60;
  $hash_ref->{61}=61;
  $hash_ref->{62}=62;
  $hash_ref->{63}=63;
  $hash_ref->{64}=64;
  $hash_ref->{65}=65;
  $hash_ref->{66}=66;
  $hash_ref->{67}=67;
  $hash_ref->{68}=68;
  $hash_ref->{69}=69;
  $hash_ref->{70}=70;
  $hash_ref->{71}=71;
  $hash_ref->{72}=72;
  $hash_ref->{73}=73;
  $hash_ref->{74}=74;
  $hash_ref->{75}=75;
  $hash_ref->{76}=76;
  $hash_ref->{77}=77;
  $hash_ref->{78}=78;
  $hash_ref->{79}=79;
  $hash_ref->{80}=80;
  $hash_ref->{81}=81;
  $hash_ref->{82}=82;
  $hash_ref->{83}=83;
  $hash_ref->{84}=84;
  $hash_ref->{85}=85;
  $hash_ref->{86}=86;
  $hash_ref->{87}=87;
  $hash_ref->{88}=88;
  $hash_ref->{89}=89;
  $hash_ref->{90}=90;
  $hash_ref->{91}=91;
  $hash_ref->{92}=92;
  $hash_ref->{93}=93;
  $hash_ref->{94}=94;
  $hash_ref->{95}=95;
  $hash_ref->{96}=96;
  $hash_ref->{97}=97;
  $hash_ref->{98}=98;
  $hash_ref->{99}=99;
  $hash_ref->{100}=100;
  $hash_ref->{101}=101;
  $hash_ref->{102}=102;
  $hash_ref->{103}=103;
  $hash_ref->{104}=104;
  $hash_ref->{105}=105;
  $hash_ref->{106}=106;
  $hash_ref->{107}=107;
  $hash_ref->{108}=108;
  $hash_ref->{109}=109;
  $hash_ref->{110}=110;
  $hash_ref->{111}=111;
  $hash_ref->{112}=112;
  $hash_ref->{113}=113;
  $hash_ref->{114}=114;
  $hash_ref->{115}=115;
  $hash_ref->{116}=116;
  $hash_ref->{117}=117;
  $hash_ref->{118}=118;
  $hash_ref->{119}=119;
  $hash_ref->{120}=120;
  $hash_ref->{121}=121;
  $hash_ref->{122}=122;
  $hash_ref->{123}=123;
  $hash_ref->{124}=124;

}

sub reverse_complement_IUPAC {
   my $dna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $revcomp;
}

sub complement_IUPAC {
   my $dna = shift;


   # complement the reversed DNA sequence
   my $comp = $dna;
   $dna =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $comp;
}

