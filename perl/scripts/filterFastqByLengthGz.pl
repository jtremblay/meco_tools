#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;

my $usage=<<'ENDHERE';
NAME:
filterFastqByLengthGz.pl

PURPOSE:
To keep reads lower or equal than --max_length <int> 

INPUT:
--infile <string>           : Sequence file
--max_length <int>          : maximum tolerated read length
--min_length <int>          : minimum tolerated read length. Default = 0
--remove_Ns                 : Flag where we remove Ns. Default = false

OUTPUT:
STDOUT                      : reads with length higher or equal to <min_length> AND lower or equal than <max_length>
--outfile_gtcutoff <string> : reads with length higher than <max_length>. Optional.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $max_length, $min_length, $outfile_gtcutoff, $remove_Ns);
my $verbose = 0;

GetOptions(
   'infile=s'           => \$infile,
   'max_length=i'       => \$max_length,
   'min_length=i'       => \$min_length,
   'outfile_gtcutoff=s' => \$outfile_gtcutoff,
   'remove_Ns'          => \$remove_Ns,
   'verbose'            => \$verbose,
   'help'               => \$help
);
if ($help) { print $usage; exit; }

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDir-filterFastqByLengthGz-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

# make pipe
my $pipe = "$tmpdir/reads.pipe";
system("mkfifo $pipe");
# gunzip to pipes
print STDERR "[DEBUG] Integrating $infile\n";
system("gunzip -c ".$infile." > $pipe &");

## MAIN

#my $size =  -s $infile or die "$! : $infile".". Cannot calculate size of file.\n";
open my $FASTQ, '<', $pipe or die $!;  
if($outfile_gtcutoff){
    open(OUTFILE, ">".$outfile_gtcutoff) or die "Can't create new file ".$outfile_gtcutoff."\n";
}

#while( <$FASTQ> ) {
#    my @lines = map scalar( <$FASTQ> ), 1 .. 4;
#    chomp @lines;
#    
#    my $header = $lines[0];  
#    my $seq = $lines[1];
#    my $qual = $lines[3];
#    if(length($seq) <= $max_length){
#        print STDOUT $header."\n".$seq."\n+\n".$qual."\n";
#    }
#}

my $Ns = 0;
my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fasta file, $pipe\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
    
    
    if($remove_Ns){
        if($curr->seq =~ m/N/) {
            $Ns++;
        }else{
            if(length($curr->seq) <= $max_length && length($curr->seq) >= $min_length){
                print STDOUT $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
            }else{
                if($outfile_gtcutoff){
                    print OUTFILE $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
                }
            }
        }
     }else{
        if(length($curr->seq) <= $max_length && length($curr->seq) >= $min_length){
            print STDOUT $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
        }else{
            if($outfile_gtcutoff){
                print OUTFILE $curr->header."\n".$curr->seq."\n+\n".$curr->qual."\n";
            }
        }
     
     }
}
close(OUTFILE);
exit;
