#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::File;
use Env qw/TMPDIR/;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Parallel::ForkManager;
use File::Temp;
#use Cache::FastMmap;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
qscoreSheets.pl

PURPOSE:
This script generates Qscore sheets using fastx_quality_stats. 
There are three options for the ouput. 

    1)One file with all qscores - each barcodes being separated by 
      a >>(header) and //(footer) string).
        >><ID>
        .
        .qscore rows
        .
        //

    2)Multiple separate files for each barcodes.

    3)One quality stats file regardless of barcodes.

INPUT:
--fastq <fastq>      :  Illumina fastq library in one file
--barcodes <fasta>   :  barcodes.
--phred <int>        :  specifiy if phred score is in phred+64 
                        (enter 64) or phred+33 (enter 33) 
--prefix <string>    :  prefix to output files
--suffix <string>    :  suffix string that will appear after barcode (optional)
--log <outfile>      :  log summary of barcode distrib (optional)
--tmp <DIR>          :  temp directory where intermediates files will be writtten. 
--num_threads <int>  :  Number of threads. Default=1.

OUTPUT:
--outfile <outfile>  :  outfile having all qscores
    OR
--outdir <outdir>    :  outdir having all qscores in separate files

NOTES:
Barcodes with a single base error will be corrected in the sequence header; 
all barcodes will be made uppercase

BUGS/LIMITATIONS:
Ambiguous nucleotides or multiple sequences per label are not supported

AUTHOR/SUPPORT:

Genomics and Microbiomes
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $barcodes_infile, $illumina_infile, $phred, $outfile, $prefix, $tmp, $outdir, $log, $suffix, $num_threads);
my $verbose = 0;
my $debug = 0;

## EXTERNAL SOFTWARE
my $barcodes_splitter_tool = `which barcodes.pl`;
my $qscore_tool = `which fastx_quality_stats`;

chomp $barcodes_splitter_tool;
chomp $qscore_tool;
die "fastx_quality_stats executable not found!\n" unless $qscore_tool;
die "barcodes binning tool not on path (barcodes.pl)\n" unless $barcodes_splitter_tool;

GetOptions(
  'fastq=s'       => \$illumina_infile,
  'barcodes=s'    => \$barcodes_infile,
  'prefix=s'      => \$prefix,
  'suffix=s'      => \$suffix,
  'outfile=s'     => \$outfile,
  'outdir=s'      => \$outdir,
  'phred=i'       => \$phred,
  'num_threads=i' => \$num_threads,
  'verbose'       => \$verbose,
  'log=s'         => \$log,
  'tmp=s'         => \$tmp,
  'help'          => \$help,
  'debug'         => \$debug
);
if ($help) { print $usage; exit; }

#VALIDATION
die("--phred provide a phred value (either 33 or 64)\n") unless $phred;
die("--fastq provide a fastq file\n") unless $illumina_infile;
die("provide either --outfile (one quality stats file for all indexes) or --outdir (each indexes have a separate quality stats file)\n") if $outfile and $outdir;
die("provide a barcode fasta file ir you want multiple file output\n") if($outdir && !$barcodes_infile);
$suffix = "_qscore.tab" unless $suffix;
$num_threads = 1 unless($num_threads);

die "No tempdir specified, please specify a \$TMPDIR variable in your environment.\n" if(!$TMPDIR);
$TMPDIR = $tmp if($tmp);
my $tmpdir = File::Temp->newdir(
	"tmpDirAmpliconTaggerQscoreXXXXXXX",
	DIR => $TMPDIR."/",
	CLEANUP => 0
);

my $tmpdir_fastq = File::Temp->newdir(
	"tmpDirAmpliconTaggerQscoreFastqXXXXXXX",
	DIR => $TMPDIR."/",
	CLEANUP => 0
);

my $tmpdir_qscores = File::Temp->newdir(
	"tmpDirAmpliconTaggerQscoreSheetsXXXXXXX",
	DIR => $TMPDIR."/",
	CLEANUP => 0
);

$tmpdir = $tmpdir."/";
$tmpdir_fastq = $tmpdir_fastq."/";
$tmpdir_qscores = $tmpdir_qscores."/";

my $separate;
if($outfile){
	$separate = 0; # One file with multiple barcodes
}elsif($outdir){
	$separate = 1; # Separate files with different barcodes
}

if($outfile && !$barcodes_infile){
	$separate = 2; # One quality score file regardless of barcodes (i.e. no barcodes).
}

## MAIN
print STDERR $tmpdir."\n" if($debug);
print STDERR $tmpdir_fastq."\n" if($debug);
print STDERR $tmpdir_qscores."\n" if($debug);
print STDERR "Separate value:\t".$separate."\n" if($debug);

if($outdir){
	$outdir = $outdir."/" if(substr($outdir,-1) ne "/");
}


# If one file with multiple barcodes.
if($separate == 0){
	open(OUT_Q, ">".$outfile) or die "Can't open file ".$outfile."\n";
	print STDERR "OUTQ opened\n" if($debug);
}

# If barcodes are to be handled.
if($separate != 2){

    # generate a hash table with barcode_header -> barcode_sequence
	my %hash = ();
	print STDERR "Placing barcodes values in hash\t".$barcodes_infile."\n" if($debug);
	my $db = Iterator::FastaDb->new($barcodes_infile) or die("Unable to open Fasta file, $barcodes_infile\n");
	while (my $seq=$db->next_seq) {
		$hash{$seq->header} = $seq->seq;
		print STDERR $seq->header()."\n".$seq->seq()."\n" if($debug); 
	}
	
	# Split library according to barcodes
	print STDERR "Splitting barcodes...\n using ".$barcodes_splitter_tool."\n" if $debug;
	if($log){
        my $cmd = $barcodes_splitter_tool." --infile ".$illumina_infile." --barcodes ".$barcodes_infile." --outdir ".$tmpdir_fastq." --log ".$log." --num_threads ".$num_threads;
        #my $cmd = $barcodes_splitter_tool." --infile_fastq ".$illumina_infile." --infile_barcodes ".$barcodes_infile." --outdir ".$tmpdir_fastq;
        print STDERR "Executing: ".$cmd."\n" if($debug);
		system($cmd);
        die "command failed: $!\n" if($? != 0);

	}else{
        my $cmd = $barcodes_splitter_tool." --infile ".$illumina_infile." --barcodes ".$barcodes_infile." --outdir ".$tmpdir_fastq." --log ".$tmpdir."/.barcodes_log.txt --num_threads ".$num_threads;
        #my $cmd = $barcodes_splitter_tool." --infile_fastq ".$illumina_infile." --infile_barcodes ".$barcodes_infile." --outdir ".$tmpdir_fastq;
        print STDERR "Executing: ".$cmd."\n" if($debug);
        system($cmd);
        die "command failed: $!\n" if($? != 0);
	}

	print STDERR $tmpdir_fastq."\n" if($debug);
		
	sleep(2);
	
	# Loop through all sequences and generate Qscore sheets.
	opendir(DIR, $tmpdir_fastq) or die "DIR ".$tmpdir_fastq." does not exists...\n";
	my @files_in_dir = readdir(DIR) or die "Can't read directory\n";
	closedir(DIR);
    
    my @files;
    foreach my $file (@files_in_dir){
        next if($file =~ /^\.$/); 
        next if($file =~ /^\.\.$/);
		# check if file don't exists or is empty
        #if current path is a file and non-empty
        if (-f $tmpdir_fastq.$file && -s $tmpdir_fastq.$file){
            push(@files, $file);
        }else{
            print STDERR $file." is empty or is not a file...\n";
        }
    }
		
	@files = sort(@files);
    print STDERR "Files to scan for Qscores:".join("\n", @files);


    system("mkdir -p ".$tmpdir_qscores);
	##======== Parallel::ForkManager starts here ========##
	my $pm = new Parallel::ForkManager($num_threads);
	my $index = 0;
	foreach(@files){
        
		$index++;
		my $pid = $pm->start($index) and next;
		print STDERR "Executing Fork process ".$index." of ".$num_threads." threads\n" if($debug);
		print STDERR "Processing file \t".$_."\n" if($debug);
		
		#$_ =~ m/(.*_.*).*/gc; # TODO find a way to get only the meaningful portion of the file name.
		my $id = $_;
        print STDERR ">>".$id."\n" if $debug;
		
		# Run qscore command, add output to a variable and finally print the content of that variable to a file.
		print STDERR "Calculate Qscores -> Tempdirqscoresubstring:\t".$tmpdir_qscores.$id."\n" if($debug);
		my $cmd = $qscore_tool." -Q ".$phred." -i ".$tmpdir_fastq.$_." -o ".$tmpdir_qscores.$id."_qscore.tab";
        print STDERR "Executing ".$cmd."\n";

        ##
	    my $attempts_left = 3;
        LOOP: {
            if(!eval { system($cmd); 1 }) {
                warn $@;
                if(--$attempts_left) {
                    warn "Retrying...\n";
                    redo;
                }else{
                    warn "Aborting.\n";
                    $pm->finish(1);
                }
            }
        }
        ##
		$pm->finish($index);
	}
	print STDERR (scalar localtime), " Waiting for some child process to finish.\n" if($debug);
	$pm->wait_all_children;
	##======== Parallel::ForkManager ends here ========##
	
    # Loop through all sequences and generate Qscore sheets.
    @files_in_dir = ();
	opendir(DIR, $tmpdir_qscores) or die "DIR ".$tmpdir_qscores." does not exists...\n";
	@files_in_dir = readdir(DIR) or die "Can't read directory\n";
	closedir(DIR);
    
    my @Qfiles;
	my @IDs;
    foreach my $file (@files_in_dir){
        next if($file =~ /^\.$/); 
        next if($file =~ /^\.\.$/);
		# check if file don't exists or is empty
        #if current path is a file and non-empty
        if (-f $tmpdir_qscores.$file && -s $tmpdir_qscores.$file){
            push(@Qfiles, $file);
            my $id = $file;
            $id =~ s{.*/}{};      # removes path  
            $id =~ s{\.[^.]+$}{}; # removes extension
            push(@IDs, $id);
        }else{
            print STDERR $file." is empty or is not a file...\n";
        }
    }
	
	#CREATE FILES FOR MULTIPLE OUTPUT OPTION
	foreach my $file (@Qfiles){
        my $content_qscores = "";
		my $id = shift(@IDs);
		print STDERR $id."\t".$file."\n" if($debug);
    
		my $an_outfile = "";
		my $a_fh;
		if($separate == 1){
		    $an_outfile = $outdir.$prefix."_".$id.$suffix;
		    $a_fh=IO::File->new(">$an_outfile") or die($!);
		}
		open FH, '<'.$tmpdir_qscores.$file or die "Can't open file ".$file."\n";
	
		$content_qscores .= ">>".$id."\n";	
		my @lines = <FH>;
		foreach my $line (@lines){
			chomp($line);
			print $a_fh $line."\n" if($separate == 1);
			$content_qscores .= $line."\n" if($separate == 0);
		}
		$content_qscores .= "//\n";
		close(FH);
		close $a_fh if($separate == 1);
		print OUT_Q $content_qscores if($separate == 0);
	}
	close(OUT_Q);

# If only one file regardless of barcode information. Run only one 
# fastx_quality_stats job for the whole undemultiplexed file.
}elsif($separate == 2){
	system $qscore_tool." -Q ".$phred." -i ".$illumina_infile." -o ".$outfile;
    die "command failed: $!\n" if($? != 0);
}

## REMOVE TEMP FILES. Cannot place this inside END subroutine
## Because END is being called for each $pm->finish function...
unless($debug){
    system("rm ".$tmpdir." -rf");
    system("rm ".$tmpdir_qscores." -rf");
    system("rm ".$tmpdir_fastq." -rf");
}
exit;

## REMOVE TEMP FILES
sub END{
}
