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
removeHeadNGz.pl

PURPOSE:
To keep reads lower or equal than --max_length <int> 

INPUT:
--infile <string>           : Sequence file

OUTPUT:
STDOUT                      : reads with length higher or equal to <min_length> AND lower or equal than <max_length>

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $max_length, $min_length, $outfile_gtcutoff);
my $verbose = 0;

GetOptions(
   'infile=s'           => \$infile,
   'max_length=i'       => \$max_length,
   'min_length=i'       => \$min_length,
   'outfile_gtcutoff=s' => \$outfile_gtcutoff,
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
print STDERR "[DEBUG] Processing $infile\n";
system("gunzip -c ".$infile." > $pipe &");

## MAIN
my $counter = 0;
#my $size =  -s $infile or die "$! : $infile".". Cannot calculate size of file.\n";
#open my $FASTQ, '<', $infile or die $!;  

#while( tell($FASTQ) < $size) {
#    my @lines = map scalar( <$FASTQ> ), 1 .. 4;
#    chomp @lines;
#    
#    my $header = $lines[0];  
#    my $seq = $lines[1];
#    my $qual = $lines[3];
#    my $first_char = substr($seq, 0, 1);
#    if($first_char eq "N"){
#        $counter++;
#        my $new_seq = substr($seq, 1);
#        my $new_qual = substr($qual, 1);
#        print STDOUT $header."\n".$new_seq."\n+\n".$new_qual."\n";
#    }
#}
my $ref_fastq_db = Iterator::FastqDb->new($pipe) or die("Unable to open Fasta file, $pipe\n");
while( my $curr = $ref_fastq_db->next_seq() ) {
    my $seq = $curr->seq;
    my $qual = $curr->qual;
    my $first_char = substr($seq, 0, 1);
    if($first_char eq "N"){
        $counter++;
        my $new_seq = substr($seq, 1);
        my $new_qual = substr($qual, 1);
        print STDOUT $curr->header."\n".$new_seq."\n+\n".$new_qual."\n";
    }
}
#close(OUTFILE);
print STDERR "[DEBUG] Removed leading N from $counter reads.\n";
exit;
