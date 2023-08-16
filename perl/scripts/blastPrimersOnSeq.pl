#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::File;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Temp;
use Bio::Seq;
use Bio::Tools::IUPAC;
use Env qw/TMPDIR/;

my $usage=<<'ENDHERE';
NAME:
blastPrimersOnSeq.pl

PURPOSE:
To blast a primer sequence (short seq) against a reference sequence database. 
Generated output is a table of blast data that is used as input in a 
following script.

INPUT:
--primer_pair <fasta_file> : fasta file containing only one primer pair (two fasta sequences).
--ref_db <fasta_file>      : database containing fasta seq on which the primer seq is to be compared
--tmpdir <dir>             : Temporary directory to store intermediate files (optional).

OUTPUT:
--outfile <tab_file>       : file containing one best hit/ref_db seq.

NOTES:
This processing time of this script is quite long and should be used with the psub wrapper.

BUGS/LIMITATIONS:
        
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $primer_pair, $ref_db, $outfile);
my $verbose = 0;

my $tmpdir = $TMPDIR;

## SCRIPTS
GetOptions(
    'primer_pair=s' => \$primer_pair,
    'ref_db=s'      => \$ref_db,
    'tmpdir=s'      => \$tmpdir,
    'outfile=s'     => \$outfile,
    'verbose'       => \$verbose,
    'help'          => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE

die("--primer_pair (fasta format) file required\n") unless $primer_pair;
die("--ref_db  fasta reference database required\n") unless $ref_db;
die("--outfile outfile required\n") unless $outfile;

## MAIN

my $start = time;
$tmpdir = $tmpdir."/" if(substr($tmpdir, -1)) ne "/";

my $blast = `which blastn`;
chomp($blast);
die("blast executable does not exist at $blast\n") if($blast eq "");

open(OUT, ">".$outfile) or die "Can't open file ".$!."\n";
print OUT "fwd_subject\tfwd_aln_length\tfwd_mismatch\tfwd_qstart\tfwd_sstart\trev_subject\trev_aln_length\trev_mismatch\trev_qstart\trev_sstart\n";   

#OPEN AND VALIDATE PRIMER PAIR FILE
open(PRIMER, $primer_pair);
my $header_counter = 0;
while(<PRIMER>){
    $header_counter++ if(substr($_,0,1) eq ">");
    die "Fasta file contains more than one sequence or is not correctly formatted\n" if($header_counter>=3); 
}
close(PRIMER);

my $primer_seq = new Iterator::FastaDb($primer_pair) or die("Unable to open Fasta file, $primer_pair \n");
my $fwd_primer = "";
my $rev_primer = "";
$header_counter = 0;
my $name_fwd;
my $name_rev;
while( my $seq = $primer_seq->next_seq() ) {
    if($header_counter == 0){
        $seq->header() =~ />(\S+)/; 
        $name_fwd = $1;
        $fwd_primer = $seq->seq();
    }
    if($header_counter == 1){
        $seq->header() =~ />(\S+)/; 
        $name_rev = $1;
        $rev_primer = $seq->seq() 
    }
    $header_counter++;
}

my @fwd_primer = generate_primer_combination($fwd_primer);
my @rev_primer = generate_primer_combination($rev_primer);

my @fwd_primer_files = write_primer_files(\@fwd_primer, $name_fwd);
my @rev_primer_files = write_primer_files(\@rev_primer, $name_rev);

my $ref_fasta_db = new Iterator::FastaDb($ref_db) or die("Unable to open Fasta file, $ref_db\n");

while( my $ref_seq = $ref_fasta_db->next_seq() ) {
    blast_primers(\@fwd_primer_files, \@rev_primer_files ,$ref_seq);
}
close(OUT);

my $end = time;
my $total = $end - $start;
#print "Processing time: ".$total."s\n";
#cleanup();

##########
## SUBS

sub blast_primers{
    my @ref_array = ($_[0], $_[1]);
    my $current_seq = $_[2];    
    
    my $strand;
    my $prefix;
    my $counter = 0;
        
    my @array_fwd = ();
    my @array_rev = ();

    foreach(@ref_array){
        my @array = @{$_};
        my $lowest_evalue = 20; #Just put an unrealistic high evalue here.

        if($counter == 0){
            $strand = "plus";
            $prefix = "FWD";
        }else{
            $strand = "minus";
            $prefix = "REV";
        }
        
        #BLAST FWD PRIMERS AND THEN REV PRIMERS
        foreach(@array){
            my $line = "";
            #TEMP FILES===========================
            my $temp_file = File::Temp->new( 
                TEMPLATE => 'tempXXXXXXXXXX',
                DIR => $tmpdir,
                SUFFIX => $prefix
            );
            print $temp_file $current_seq->header()."\n".$current_seq->seq()."\n";
            #=====================================    
        
            my $command = $blast." -task blastn -query ".$_." -subject ".$temp_file." -strand ".$strand." -max_target_seqs 1 -word_size 4 -outfmt 6";
            
            #REDIRECT STDOUT TO VARIABLE
            open my $fh, $command." |" or die "cannot run command: $!";{
                  local $/;
                  $line = <$fh>;
            }
            close $fh;    
            #print $line;
            
            #Pass on to next seq if no hit.
            next if($line eq "");
    
            #PLACE VALUES IN A HASH AND ONLY KEEP HIT WITH BEST EVALUE.
            chomp($line);
            my @row = split(/\t/, $line);
            if($row[10] < $lowest_evalue and $row[10] != 0){
                $lowest_evalue = $row[10];
                
                if($counter == 0){
                    $array_fwd[0] = $row[1]; #subject
                     $array_fwd[1] = $row[3]; #align length
                     $array_fwd[2] = $row[4]; #mismatch
                     $array_fwd[3] = $row[8]; #q start
                     $array_fwd[4] = $row[9]; #s start
                }elsif($counter == 1){
                    $array_rev[0] = $row[1]; #subject
                     $array_rev[1] = $row[3]; #align length
                     $array_rev[2] = $row[4]; #mismatch
                     $array_rev[3] = $row[8]; #q start
                     $array_rev[4] = $row[9]; #s start
                }
            }        
        }
        $counter++;
    }
    print OUT $array_fwd[0]."\t".$array_fwd[1]."\t".$array_fwd[2]."\t".$array_fwd[3]."\t".$array_fwd[4]."\t".$array_rev[0]."\t".$array_rev[1]."\t".$array_rev[2]."\t".$array_rev[3]."\t".$array_rev[4]."\n" if(@array_fwd == 5 and @array_rev == 5);
        
}
#WRITE UNAMIGUOUS PRIMER SEQS TO A FILE
sub write_primer_files{
    my @array = @{$_[0]};
    my $prefix = $_[1];

    my @files = ();

    my $counter=1;
    foreach(@array){
        
        #TEMP_FILES===========================
        my $temp_file = File::Temp->new( 
            TEMPLATE => 'primerXXXXXXXXXX',
            DIR => $tmpdir,
            SUFFIX => $prefix."_".$counter
        );
        print $temp_file ">".$prefix."_".$counter."\n".$_."\n";
        push(@files, $temp_file);
        #=====================================    
        $counter++;
    }    
    return @files;
}

#GENERATE ALL POSSIBLE PRIMER COMBINATION
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

## FINAL CLEANUP
sub cleanup{
    system("rm -rf ".$tmpdir."\/*");
    mkdir $tmpdir unless -d $tmpdir;
}

exit;
