#!/usr/bin/env perl
=head1 DESCRIPTION
  Purpose: eliminate duplicate sequences
  Required: linux, awk, sort, uniq, perl
  The script starts by finding uniq sequences. This step is pretty fast. 
  The second step is checking reverse complement and subsequence, 
  the O(n x n) algorithm will take a loooooong time to finish (-all option). 
  A heuristic algorithm (default) is used to split contigs into 10000 sequence per pile and 
  check sequences' rc and sub for each pile. So, the time is O(pile). 
=head2 The States of Work
  initalized on 20110315
  add -rc and -sub option, @ 201104011
  add -all, @ 20110414
  add -cov @ 20110420
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Temp;
use Env qw/TMPDIR/;

$|=1;
my $file;
my $output;
my $rc;
my $sub;
my $all;
my $cov=1;
GetOptions('o=s'    => \$output,
           'i=s'    => \$file,
           'rc'     => \$rc,
           'sub'    => \$sub,
           'cov=f'  => \$cov,
           'all'    => \$all,
           'help|?' => sub{Usage()}
           );  
           
&Usage unless ($output && $file);

my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$file", qr/\.[^.]*/ );    
my ( $o_file_name, $o_path, $o_suffix ) = fileparse( "$output", qr/\.[^.]*/ );

# format into tab-delimited
my $table_seq = &fasta_table($file,$cov);

# find unique
#my @seqs=`awk -F"\t" '{print \$2}' $table_seq | sort -r | uniq`; 
#my $tmpsort = "${o_path}sort" . time();
#my $tmpsort = "/projectb/scratch/bfoster/tmp/" . basename("$output") . ".derepSort." . time();
my $tmpsort = File::Temp->newdir(
    "dereplicateContigsXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0
);
my $tmpsport .= basename("$output") . ".derepSort." . time();

while(-d $tmpsort){
    $tmpsort .= ".1";
}
system("mkdir $tmpsort");
my @seqs=`awk -F"\t" '{print length(\$2),\$2}' $table_seq | sort -nrk1,2 -T $tmpsort | awk '{print \$2}' | uniq`; 
system("rm -rf $tmpsort");

my $dup;
my $num_per_pile = 10000;

$num_per_pile = 10000000000 if ($all);

$dup=&check_rc_sub_heuristic(@seqs) if ($rc or $sub);

# clean sequences print out 
open (OUT,">$output") or die "can't open $output for writing\n";
my $n='Sequence0000000001';
my $remain_num = 0;
for my $i (0..$#seqs)
{
   if (! $dup->{$i}){
     my $out_seq = $seqs[$i];
     $out_seq =~ s/ //g;
     $out_seq =~ s/(.{100})/$1\n/g;
     chomp $out_seq;
     print OUT ">$n\n$out_seq\n";
     $n++;
     $remain_num++;
   }
}
close OUT;

print " Remaining sequence number: ",$remain_num,"\n";
unlink $table_seq;

&print_run_time;

sub check_rc_sub {
 # slow algorithm
  my @seqs = @_;
  my %dup;
  my $j;
  my $short_seq;
  my $rc_long_seq;
  my $i;
  print " Unique sequence number: " , scalar (@seqs), "\n";
  for $i (0..$#seqs)
  {
#     print "\r Checking subsequence or reverse compelement ", ($i + 1);
     next if ($dup{$i});
     my $seq = $seqs[$i];
     chomp $seq;
     $rc_long_seq = &ReverseComplement($seq) if ($rc);
     
     for $j (($i+1)..$#seqs)
     {
         
         next if ($dup{$j});
         $short_seq = $seqs[$j];
         chomp $short_seq;
         if ($sub){    
             $dup{$j} =1 if ($seq =~ /$short_seq/);
         }
         if ($rc){
             #chomp $rc_long_seq;
             $dup{$j} =1 if ($rc_long_seq eq "$short_seq" && !$dup{$j});
             $dup{$j} =1 if ($rc_long_seq =~ /$short_seq/ && $sub && !$dup{$j});
         }      
     }
     
  }
  print "\n";
  return \%dup;
}

sub check_rc_sub_heuristic {
 # fast algorithm (heristic algorithm) but 100% check all vs all comparison.
  my @seqs = @_;
  my %dup;
  my $i;
  my $uniq_sequence_num=scalar (@seqs);
  print " Uniq seqeunce number: " ,  $uniq_sequence_num, "\n";
  ## heuristic checking algorithm (10000 sequences per pile)
  my $piles_num = int($uniq_sequence_num/$num_per_pile);
  my $residule = $uniq_sequence_num % $num_per_pile;
  if ($piles_num>0){
      for my $pile(1..$piles_num)
      {
        my @keeps=();
        my $pile_start=($pile-1)*$num_per_pile;
        my $pile_end=$pile*$num_per_pile-1;
        for $i ($pile_start..$pile_end)
        {
 #          print "\r Checking subsequence or reverse compelement ", ($i + 1);
           my $seq = $seqs[$i];
           chomp $seq;
           my $rc_seq = &ReverseComplement($seq) if ($rc);
           #my $seq_len = length($seq);
           foreach my $j (@keeps)
           {
               last if ($i==0);
               next if ($dup{$j});
               my $long_seq = $seqs[$j];
               chomp $long_seq;
               if ($sub){    
                   $dup{$i} = 1 if ($long_seq =~ /$seq/);
               }
               last if ($dup{$i});
               if ($rc){            
                   $dup{$i} = 1 if ($long_seq eq "$rc_seq" );
                   $dup{$i} = 1 if ($long_seq =~ /$rc_seq/ && $sub);
               }   
               last if ($dup{$i});
           }
           push @keeps, $i;
        }
      }
  } # end if ($piles_num>0)
  ## checking residules
  my @keeps=();
  my $pile_start=$piles_num*$num_per_pile;
  my $pile_end=$piles_num*$num_per_pile+$residule-1;
  for $i ($pile_start..$pile_end)
  {
#       print "\r Checking subsequence or reverse compelement ", ($i + 1);
       my $seq = $seqs[$i];
       chomp $seq;
       my $rc_seq = &ReverseComplement($seq) if ($rc);
       #my $seq_len = length($seq);
       foreach my $j (@keeps)
       {
           last if ($i==0);
           next if ($dup{$j});
           my $long_seq = $seqs[$j];
           chomp $long_seq;
           if ($sub){    
               $dup{$i} = 1 if ($long_seq =~ /$seq/);
           }
           last if ($dup{$i});
           if ($rc){            
               $dup{$i} = 1 if ($long_seq eq "$rc_seq" );
               $dup{$i} = 1 if ($long_seq =~ /$rc_seq/ && $sub);
           }   
           last if ($dup{$i});
       }
       push @keeps, $i;
  }
  
  
  print "\n";
  return \%dup;
}

sub fasta_table{
   my ($file, $cutoff)=@_;
   my $seq;
   my $id;
   my $seq_number = 0;
   my $above_coverage_num = 0;
   open (my $fh,"< $file") or die "can't open $file for reading\n";
   my $tmp = "$o_path/tmp$$.tbl";
   open (OUTTMP, ">$tmp") or die "can't open $tmp for writing\n";
   while(<$fh>)
   {
      chomp;
      if(/^>(.*)/) 
      {  
        if ($seq){
             if(&coverageTest($id,$cutoff) == 0)
             {
                print OUTTMP "$id\t$seq\n";
                $above_coverage_num++;
             }   
        }
        $id =$1;
    	$seq ="";
    	$seq_number++;
      }
      else 
      {
        $seq .= $_;
      }
   }
       if ($seq){
           if(&coverageTest($id,$cutoff) == 0)
           {
              print OUTTMP "$id\t$seq\n";
              $above_coverage_num++;
           }
       }
    close $fh;
    close OUTTMP;
    print " Total sequence number: $seq_number\n";
    print " More than $cutoff coverage sequence: $above_coverage_num\n";
    return ($tmp);
}

sub ReverseComplement{
	my $dna = $_[0];
    my $ReverseCompSeq = reverse ($dna);
	$ReverseCompSeq =~ tr/atgcrywsmkATGCRYWSMK/tacgyrswkmTACGYRSWKM/;
	return($ReverseCompSeq);
}


sub Usage 
{
  print "perl $0 [-rc -sub -all -cov] -o <output> -i <fasta>\n";
  print "     -i    input fasta file\n";
  print "     -o    output fasta file\n";
  print "     -rc   check reverse complement sequence\n";
  print "     -sub  check 100% identity subsequence\n"; 
  print "     -all  [work with -rc or -sub]\n";
  print "           use all by all algorithm intead of heuristic algorithm\n";
  print "           it will be much clean but take a long time if input a lot of sequences\n";
  print "     -cov  contig coverage cutoff for SOAPdenovo and Velvet (default: 1.0)\n";
  print "     -h    print usage\n";
  exit;
}

sub print_run_time {
  # Print runtime #
  my $run_time = time() - $^T;
  printf("\n Total running time: %02d:%02d:%02d\n\n", int($run_time / 3600), int(($run_time % 3600) / 60), 
  int($run_time % 60));
}

 
sub coverageTest{
# by mscholz
# pass the coverage cutoff return 0;
                my ($test,$cutoff) = @_;
#SOAPdenovo coverage report:
                if($test =~ /cvg/){
                                my @temp =split /_/, $test;
                                if($temp[1] >= $cutoff){return 0;}
                                else {return 1;}
                }
#velvet coverage:
                elsif($test =~ /cov/){
                                my @temp = split /_/, $test;
                                if($temp[-1] >= $cutoff){return 0;}
                                else{return 1;}
                }
#return ALL 454 as OK
                elsif($test =~ /numreads/){
                                return 0;
                }
#return weird contigs names as normal
                else{
                                return 0;             
                }
}

=head1 COPYRIGHT

Copyright (c) 2011 LANS, LLC All rights reserved

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Chien-Chi Lo <chienchi at lanl.gov>

=cut
