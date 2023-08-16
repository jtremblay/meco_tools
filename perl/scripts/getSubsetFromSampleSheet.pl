#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use List::Util qw(sum);
use Getopt::Long;
use File::Which;
use File::Temp;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use Statistics::Descriptive;
use SampleSheet;
use File::Basename;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
getSubsetFromSampleSheet.pl

PURPOSE:
Sample a random number of reads for library specified in sample sheet. The purpose of this
is to test a pipeline with a reduced dataset of reads.

INPUT:
--indir <string>               : If data to be assembled is from our GQIC internal
                                 libraries. The wrapper will look for the DIR 
                                 <INDIR>/raw_reads/
--sampleSheet <string>         : GQIC formatted sample sheet.
--num_threads <int>            : Number of threads
--n <int>                      : Number of reads you want to subsample to get
                                 paired-end assembly stats.
--QStatPlots                   : If quality stats plots are to be generated.
                                 Plot will be done from subset reads only.
--paired
OR
--single      

OUTPUT:
--outdir <string>              :  output directory

NOTES:

BUGS/LIMITATIONS:
The flash binary has to be in your $PATH.

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com
ENDHERE

## OPTIONS
my ($help, $paired, $single,$num_threads, $outdir, $qscores, $start_at, $n, $sampleSheet, $indir);
my $verbose = 0;

GetOptions(
	'indir=s'		=> \$indir,
	'sampleSheet=s'	=> \$sampleSheet,
	'outdir=s'		=> \$outdir,
	'num_threads=i' => \$num_threads,
	'n=i'			=> \$n,
	'QstatPlots'	=> \$qscores,
	'paired'		=> \$paired,
	'single'		=> \$single,

	'start_at=i'	=> \$start_at,
    'verbose' => \$verbose,
    'help' => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--outdir arg missing...\n") unless($outdir);
#die("--num_threads arg required\n") unless($num_threads);
#die("no \$TMPDIR variable defined\n") unless($TMPDIR);
die("--outdir does not exists...\n") unless(-d $outdir);

#my $tmpdir = File::Temp->newdir(
#    "tmpDirFlashWrapperXXXXXXX",
#    DIR => $TMPDIR."/",
#    CLEANUP => 0 # 1-Delete after execution. 0-Do not delete.
#);


## MAIN
$n = 5000 unless($n);
$start_at = 0 unless($start_at);
$num_threads = 1 unless($num_threads);
my $tmpdir = $TMPDIR;
mkdir $tmpdir unless -d $tmpdir;
mkdir $outdir."/raw_reads_subset" unless -d $outdir."/raw_reads_subset";
mkdir $outdir."/Qscores" unless -d $outdir."/Qscores";

my %hash1;
my %hash2;
my %hashAssembled;
my %hashLog;

if($indir && $sampleSheet){
	my $rHoAoH_sampleInfo  = SampleSheet::parseSampleSheetAsHash($sampleSheet);
	
	# Loop through each samples in sample sheet
	my $rawReadDir    = $indir."/raw_reads";
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
		my $R1_gz;
		my $R2_gz;
		my $R1_name;
		my $R2_name;
		my $R_subset_dir;

		for my $rH_laneInfo (@$rAoH_sampleLanes) { #this rH_laneInfo contains the complete line info from the sample sheet for this sample.
			$R1_gz = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'} if($single);
		 	$R2_gz = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read2File'} if($paired);
			$R1_name = $rH_laneInfo->{'read1File'} if($single);
		 	$R2_name = $rH_laneInfo->{'read2File'} if($paired);	
			$R_subset_dir = $outdir."/raw_reads_subset/".$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
			mkdir $R_subset_dir unless -d $R_subset_dir;
			mkdir $outdir."/raw_reads_subset/".$sampleName unless -d $outdir."/raw_reads_subset/".$sampleName;
			mkdir $outdir."/raw_reads_subset/".$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} unless -d $outdir."/raw_reads_subset/".$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
		}

		# Uncompress
		$R1_name =~ s/\.gz// if($single);
		$R2_name =~ s/\.gz// if($paired);
		my $R1_fastq = $tmpdir."/".$R1_name if($single);
		my $R2_fastq = $tmpdir."/".$R2_name if($paired);
		system("gunzip -c ".$R1_gz." > ".$R1_fastq) if($single and $start_at <= 0);
		$? != 0 ? die "command failed: $!\n" : print STDERR $R1_gz." decompressed... into ".$R1_fastq."\n" if($verbose);
		system("gunzip -c ".$R2_gz." > ".$R2_fastq) if($paired and $start_at <= 0);
		$? != 0 ? die "command failed: $!\n" : print STDERR $R2_gz." decompressed... into ".$R2_fastq."\n" if($verbose);

		print STDERR $R1_fastq."\n" if($verbose && $single);
		print STDERR $R2_fastq."\n" if($verbose && $paired);

		# GetSubset
		getSubsetSingle($R1_fastq, $R_subset_dir."/".$R1_name) if($single);
		getSubsetPaired($R1_fastq, $R2_fastq, $R_subset_dir."/".$R1_name, $R_subset_dir."/".$R2_name) if($paired);

		# Compute Qstats and store qstat sheet in a hash for later retrival when building qscore plots.
		$hash1{'fastq'}{$R1_name} = $R_subset_dir."/".$R1_name if($single);
		$hash2{'fastq'}{$R2_name} = $R_subset_dir."/".$R2_name if($paired);
	}
}
my ($qscores1, $qscores2, $qscoresAss) = generateQscoreSheet(\%hash1, \%hash2, \%hashAssembled) if($qscores);	
generateQscorePlots($qscores1, $qscores2, $qscoresAss) if($qscores);

exit;

## SUBROUTINES

sub generateQscorePlots{
	my $qs1 = shift;
	my $qs2 = shift;
	my $qsA = shift;

	# Generate plots for R1 and R2
	my $cmd = "~/build/tools/generateQscorePlots.pl";
	$cmd .= " --infile_1 ".$qs1;
	$cmd .= " --infile_2 ".$qs2;
	$cmd .= " --name R1_R2";
	$cmd .= " --pdf ".$outdir."/R1R2Qplots.pdf";
	$cmd .= " --display 1";
	$cmd .= " --paired";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDERR "R1R2 Qplots generated...\n" if($verbose);

	# Generate plots for assembled reads.
	$cmd = "~/build/tools/generateQscorePlots.pl";
	$cmd .= " --infile_1 ".$qsA;
	$cmd .= " --name Assembled_reads";
	$cmd .= " --pdf ".$outdir."/assembledReadsQplots.pdf";
	$cmd .= " --display 1";
	$cmd .= " --single";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDERR "Extended fragments Qplots generated...\n" if($verbose);
}

sub generateQscoreSheet{

	my( $href1, $href2, $hrefAssembled) = @_;
	my @prefixes =("R1", "R2", "Assembled");
	
	my @AoH = ($href1, $href2, $hrefAssembled);
	my @Qscores;	

	# Get file names (with full path).
	for my $href ( @AoH ) {
		my %hash = %$href;	
		my $prefix = shift @prefixes;
		my $files = "";
		foreach my $key (keys %hash) {
			foreach my $key2 (keys %{ $hash{$key} }) {
				$files .= " --fastq ".$hash{'fastq'}{$key2};
			}
			my $cmd = "~/build/tools/generateQscoreSheets.pl ";
			$cmd .= $files;
			$cmd .= " --outfile ".$outdir."/Qscores/qscore_".$prefix.".txt";
			$cmd .= " --phred 33";
			$cmd .= " --num_threads ".$num_threads;
			system($cmd);
			$? != 0 ? die "command failed: $!\n" : print STDERR "Reads subset generated...\n" if($verbose);
			push(@Qscores, $outdir."/Qscores/qscore_".$prefix.".txt");
	    }
	}
	return @Qscores;	
}

sub getSubsetPaired{

	my $infile_1 = shift;
	my $infile_2 = shift;
	my $outfile_1 = shift;
	my $outfile_2 = shift;

	my $R1Prefix = basename($infile_1);
	my $R2Prefix = basename($infile_2);
	$R1Prefix =~ s/\.pair1\.fastq//;	
	$R2Prefix =~ s/\.pair2\.fastq//;
	my $ASSPrefix = $R1Prefix;
	$ASSPrefix =~ s/\.pair1\.fastq//;

	# Get subset of reads (random).
	my $cmd = "~/build/tools/getSubset.pl ";
	$cmd .= " --infile_1 ".$infile_1;
	$cmd .= " --infile_2 ".$infile_2;
	$cmd .= " --outfile_1 ".$outfile_1;
	$cmd .= " --outfile_2 ".$outfile_2;
	$cmd .= " --n ".$n;
	system($cmd) if($start_at <= 0);
	$? != 0 ? die "command failed: $!\n" : print STDERR "Reads subset generated...\n" if($verbose);	
}

sub getSubsetSingle{

	my $infile = shift;
	my $outfile = shift;

	# Get subset of reads (random).
	my $cmd = "~/build/tools/getSubset.pl ";
	$cmd .= " --infile ".$infile;
	$cmd .= " --outfile ".$outfile;
	$cmd .= " --n ".$n;
	print $cmd."\n";
	system($cmd) if($start_at <= 1);
	$? != 0 ? die "command failed: $!\n" : print STDERR "Reads subset generated...\n" if($verbose);	
}

# Remove temp files
sub END{
	#system("rm ".$tmpdir." -rf");
}

