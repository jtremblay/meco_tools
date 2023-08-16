#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Env qw(TMPDIR);
use File::Temp;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
generateGenbankToFastaHeader.pl

PURPOSE:

INPUT:
--indir <string> : Sequence file
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
   'indir=s'  => \$indir,
   'verbose'  => \$verbose,
   'help'     => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $tmpdir = File::Temp->newdir(
    "tmpDir-genbankToHeader-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if(substr($filename, -4) eq ".fna"){
            #my $pipe = "$tmpdir/reads.pipe";
            #system("mkfifo $pipe");

            # gunzip to pipes
            print STDERR "[DEBUG] processing $_\n";
            #system("gunzip -c ".$_." > $pipe &");

            my $filename2 = $filename;
            $filename2 =~ s/\.fna//;
            
            my $ref_fasta_db = Iterator::FastaDb->new($filename) or die("Unable to open Fasta file, $filename\n");
            while( my $curr = $ref_fasta_db->next_seq() ) {
                my ($header) = $curr->header =~ m/^>(\S+).*$/;
                print STDOUT $filename2."\t".$header."\n";
            }
            #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
            #system("rm $pipe");
		}
	}
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);


