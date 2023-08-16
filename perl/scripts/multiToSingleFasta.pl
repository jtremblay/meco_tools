#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use Cwd 'abs_path';
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
multiToSingleFasta.pl

PURPOSE:

INPUT:
--indir <string>     : Dir containing multifasta files
--extension <string> : default = fna
--recursive <string> : default = false

OUTPUT:
--outdir <string> : Dir where single fasta files will be created.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $indir, $outdir, $extension, $recursive);
my $verbose = 0;

GetOptions(
   'indir=s' 	 => \$indir,
   'outdir=s' 	 => \$outdir,
   'extension=s' => \$extension,
   'verbose' 	 => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

$indir = abs_path($indir);
$outdir = abs_path($outdir);

$extension = "fna" unless($extension);
$recursive = "false" unless($recursive);

print STDERR "indir: ".$indir."\n";
print STDERR "outdir: ".$outdir."\n";

## MAIN
sub eachFile{
   my $filename = $_;
   my $fullpath = $File::Find::name;
   #remember that File::Find changes your CWD, 
   #so you can call open with just $_
    
   print STDERR $fullpath."\n";
   print STDERR $filename."\n";
   if (-e $filename) { 
    
      if(substr($filename, -(length($extension)+1)) eq $extension){
         print STDERR "Parsing ".$fullpath." into single fasta file archive...\n";
         #system("pigz -p ".$num_threads." ".$fullpath);
         #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
         my $ref_fasta_db = Iterator::FastaDb->new($fullpath) or die("Unable to open Fasta file, $fullpath\n");
         my $outfile = "$outdir/$filename";
         open(OUT, ">".$outfile) or die "Can't open $outfile\n";
         my $header = 0;
         my $seq = "";
         while( my $curr = $ref_fasta_db->next_seq() ) {
             #if($header == 0){
             #  my $curr_header = $filename;
             #  $curr_header =~ s/\.$extension//g;
             #  print OUT ">".$curr_header."\n";
               $seq .= $curr->seq;
               #$header++; 
             #}else{
             #  $seq .= $curr->seq;
             #  $header++; 
             #}
         }

         my $i = 0;
         while($i < length($seq)){
             #print OUT ">".$filename."\n";
             print OUT substr($seq, $i, 80)."\n";
             $i = $i + 80;
         }
         close(OUT);
      }
   }
}

sub findNonRecursive(){
    my @filenames = glob($indir."/*.".$extension);
    foreach my $filename (@filenames){
        print STDERR $filename."\n";    
        if (-e $filename) {
            my $basename = basename($filename);
            next if($basename  eq "all_bins.fa");
            #print STDERR substr($basename, -(length($extension)+1))."\n";
            
           if(substr($basename, -(length($extension)+1)) eq ".".$extension){
                print STDERR "Parsing ".$filename." into single fasta file archive...\n";
                #system("pigz -p ".$num_threads." ".$fullpath);
                #$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
                my $ref_fasta_db = Iterator::FastaDb->new($filename) or die("Unable to open Fasta file, $filename\n");
                my $outfile = "$outdir/$basename";
                print STDERR $outfile."\n";
                open(OUT, ">".$outfile) or die "Can't open $outfile\n";
                my $header = 0;
                my $seq = "";
                while( my $curr = $ref_fasta_db->next_seq() ) {
                    #if($header == 0){
                    #  my $curr_header = $filename;
                    #  $curr_header =~ s/\.$extension//g;
                    #  print OUT ">".$curr_header."\n";
                    $seq .= $curr->seq;
                    #$header++; 
                    #}else{
                    #  $seq .= $curr->seq;
                    #  $header++; 
                    #}
                }

                my $i = 0;
                print OUT ">".$basename."\n";
                while($i < length($seq)){
                    print OUT substr($seq, $i, 80)."\n";
                    $i = $i + 80;
                }
                close(OUT);
            }
        }
    }
}

if($recursive eq "true"){
    find (\&eachFile, $indir);
}elsif($recursive eq "false"){
    findNonRecursive();
}
