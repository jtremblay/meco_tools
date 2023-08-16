#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Env qw(TMPDIR);
use File::Temp;
use File::Find;
use Cwd qw(cwd);

my $usage=<<'ENDHERE';
NAME:
correctSoilBacteriaContigsName.pl

PURPOSE:

INPUT:
--indir <string>  : Indir with sequence files

OUTPUT:
--outdir <string> : Outdir

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $outdir, $link);
my $verbose = 0;

GetOptions(
   'indir=s'   => \$indir,
   'outdir=s'  => \$outdir,
   'link=s'    => \$link,
   'verbose'   => \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN

$outdir = cwd()."/".$outdir;

my $tmpdir = File::Temp->newdir(
    "tmpDir-genbankToHeader-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);
my %hash;
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
            #system("gunzip -c ".$_." > $pipe &");

            my $filename2 = $filename;
            #$filename2 =~ s/\.fna//;
            $filename2 =~ s{^.*/}{};     # remove the leading path
            $filename2 =~ s{\.[^.]+$}{}; # remove the extension
            
            print STDERR "[DEBUG] processing $_   output: $outdir/$filename2.fna\n";
           
            open(OUT, ">"."$outdir/$filename2.fna") or die "Can't open $outdir/$filename2.fna  $!\n";
            my $ref_fasta_db = Iterator::FastaDb->new($filename) or die("Unable to open Fasta file, $filename\n");
            my $i = 1;
            while( my $curr = $ref_fasta_db->next_seq() ) {
                my ($header) = $curr->header =~ m/^>(\S+) .*$/;
                #_L006_R1_001_paired_trimmed_paired_contig
                #$header =~ s/_L\d\d\d_R\d_\d\d\d_paired_trimmed_paired_contig_//;
                #print OUT ">".$header."\n".$curr->seq."\n";
                print OUT ">".$filename2."_contig_id_".$i."\n".$curr->seq."\n";
                $hash{$filename2}{"contig_id_".$i} = $header;
                $i++;
            }
            close(OUT);
            #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
            #system("rm $pipe");
		}
	}
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

#finally, print link
open(OUT, ">".$link) or die "Cant' open $link\n";
foreach my $filename2 (keys %hash){
    foreach my $gene_id (keys %{ $hash{$filename2} }) {
        print OUT $filename2."\t".$gene_id."\t".$hash{$filename2}{$gene_id}."\n";
    }
}
close(OUT);
