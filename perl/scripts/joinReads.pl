#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;
use Cwd 'abs_path';

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
joinReads.pl

PURPOSE:
Takes paired reads, R1 and R2 and generate what we refer to as "joined reads".
Which means R1 and R2 placed next to one another separated by 8Ns.

INPUT:
--infile_R1 <string>     : sequence file
--infile_R1 <string>     : sequence file
--Ns <int>               : Number of Ns to put between reads (default = 8).				

OUTPUT:
--outfile <stirng> : outfile where fasta (gz format) are written.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile_R1, $infile_R2, $Ns, $outfile);
my $verbose = 0;

GetOptions(
   'infile_R1=s'     => \$infile_R1,
   'infile_R2=s'     => \$infile_R2,
   'Ns=i'            => \$Ns,
   'outfile=s' => \$outfile,
   'verbose' 	      => \$verbose,
   'help'            => \$help
);
if ($help) { print $usage; exit; }

# reverse completment IUPAC
sub revcomp {
   my $dna = shift;

   # reverse the DNA sequence
   my $revcomp = reverse($dna);

   # complement the reversed DNA sequence
   $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
   return $revcomp;
}

my $tmpdir = File::Temp->newdir(
    "tmpdir-joinReads-XXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

my %hash;

$Ns = 8 unless($Ns);
die "--infile_R1 <string> missing...\n" unless($infile_R1);
die "--infile_R2 <string> missing...\n" unless($infile_R2);

print STDERR $infile_R1."\n";
print STDERR $infile_R2."\n";
 
# Declare output;
my $z = new IO::Compress::Gzip $outfile or die "gzip failed: $GzipError\n";

# make pipe
my $pipe1 = "$tmpdir/reads1.pipe";
my $pipe2 = "$tmpdir/reads2.pipe";
system("rm $pipe1");
system("rm $pipe2");
system("mkfifo $pipe1");
system("mkfifo $pipe2");

# gunzip to pipes
system("gunzip -c ".$infile_R1." > $pipe1 &");
system("gunzip -c ".$infile_R2." > $pipe2 &");

my $string = ""; for(my $i=0; $i<$Ns; $i++){ $string .= "N";}  

my $ref_fastq_db1 = Iterator::FastqDb->new($pipe1) or die("Unable to open Fastq file, $pipe1\n");
my $ref_fastq_db2 = Iterator::FastqDb->new($pipe2) or die("Unable to open Fastq file, $pipe2\n");
while( my $curr1 = $ref_fastq_db1->next_seq() ) {
   my $curr2 = $ref_fastq_db2->next_seq();
   
   $z->print(">".$curr1->base."\n".$curr1->seq.$string.revcomp($curr2->seq)."\n");
}
#system("rm $pipe1");
#system("rm $pipe2");
$z->close();

