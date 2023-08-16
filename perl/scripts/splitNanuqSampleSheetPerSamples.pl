#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use Cwd;
use File::Path qw(make_path);
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
splitNanuqSampleSheetPerSamples.pl

PURPOSE:
Takes a sample sheet in nanuq format and create new sample sheets (1 per sample).
Also write a cat script to obtain one fastq per file. Generates a readset as well
(to use with a pipeline).

INPUT:
--infile <string>  : nanuq sample sheet. HISEQ READ SET SHEET, not sample sheet...
        
OUTPUT:
--readset <string> : File into which will be written read sets
                     (mandatory for pipeline usage).
--makeSL <string>  : Directory in which to make symlinks.

STDOUT:            : Will write cat commands to STD output.

NOTES:
Make a folder called raw_data and put all raw reads fetched from GQ in there.
Then make a folder called raw_reads.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $readset, $makeSL);
my $verbose = 0;

GetOptions(
    'infile=s'  => \$infile,
    'readset=s' => \$readset,
    'makeSL=s'  => \$makeSL,
    'verbose'   => \$verbose,
    'help'      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "missing --infile\n" unless($infile);
die "missing --readset\n" unless($readset);
die "missing --makeSL\n" unless($makeSL);

my $raw_reads_prefix = "raw_data/";
my $dir = getcwd;
$dir .= "/csvs";
make_path($dir);

## MAIN
my %hash;
my %catHash;
my $header;
open(IN, "<".$infile);
while(<IN>){
   chomp;
   if($. == 1){
      $header = $_;  
   }
   my @row = split(/,/, $_);
   
   my $z=0;
   foreach my $el (@row){
      print STDERR "$z: $el\n";
      $z++;
   }


   next if $row[0] =~ m/"Name"/;
   $hash{$row[0]} .= $_."\n";
   my $file = $row[27];
   $file =~ s/"//g;
   my $name = $row[0];
   $name =~ s/"//g;
   #print STDERR "file: $file\n";
   if(-e $raw_reads_prefix."/".$file."_R1.fastq.gz"){
      $catHash{$name}{R1}             .= " ".$raw_reads_prefix.$file."_R1.fastq.gz";
      ($catHash{$name}{libraryBarcode} = $row[1]) =~ s/"//g; 
      ($catHash{$name}{run}            = $row[20]) =~ s/"//g;
      ($catHash{$name}{lane}           = $row[22]) =~ s/"//g;
      ($catHash{$name}{adaptor1}       = $row[17]) =~ s/"//g;
      ($catHash{$name}{adaptor2}       = $row[18]) =~ s/"//g;
      ($catHash{$name}{qualityOffset}  = $row[21]) =~ s/"//g;
      ($catHash{$name}{readset}        = $row[27]) =~ s/"//g;
      #$catHash{$name}{bam}            = 
  }
  if(-e $raw_reads_prefix."/".$file."_R2.fastq.gz"){
      $catHash{$name}{R2}             .= " ".$raw_reads_prefix.$file."_R2.fastq.gz";
   }
}
close(IN);

print STDERR Dumper(\%catHash);


for my $key (keys %hash){
  my $filename = $dir."/".$key;
  $filename =~ s/"//g;

  my $string = $hash{$key};
  chomp($string);

  open(OUT, ">".$filename.".csv");
  print OUT $header."\n";
  print OUT $string."\n";
  close(OUT);
}

# Then loop through cat hash to print to STDOUT
open(READSET, ">".$readset) or die "Can't open $readset\n";
print READSET "Sample\tReadset\tRunType\tBAM\tFASTQ1\tFASTQ2\tLibrary\tRun\tLane\tAdaptor1\tAdaptor2\tQualityOffset\tBED\n";

system("mkdir -p ".$makeSL);

for my $key (keys %catHash){
   print READSET $key."\t".$key."_LANE_X\tPAIRED_END\t"; 
  
   my $R1_string = $catHash{$key}{R1};
   #$R1_string =~ s/^\s+|\s+$//g;
   my $R2_string = $catHash{$key}{R2};
   #$R2_string =~ s/^\s+|\s+$//g;

   print STDOUT "cat ".$R1_string." > ./raw_reads/".$key."_R1.fastq.gz\n";
   print STDOUT "cat ".$R2_string." > ./raw_reads/".$key."_R2.fastq.gz\n";
   
   #system("ln -sf ".$catHash{$key}{R1}." ".$makeSL."/".$key."_R1.fastq.gz");
   #system("ln -sf ".$catHash{$key}{R2}." ".$makeSL."/".$key."_R2.fastq.gz");
   
   print READSET $key.".bam\t$raw_reads_prefix".$key."_R1.fastq.gz\t$raw_reads_prefix".$key."_R2.fastq.gz\t".$catHash{$key}{libraryBarcode}."\t".$catHash{$key}{run}."\t".$catHash{$key}{lane}."\t".$catHash{$key}{adaptor1}."\t".$catHash{$key}{adaptor2}."\t".$catHash{$key}{qualityOffset}."\t".$key.".bed\n";

   if($makeSL){
      symlink(abs_path($key."_R1.fastq.gz"), $makeSL."/".$key."_R1.fastq.gz"); 
      symlink(abs_path($key."_R2.fastq.gz"), $makeSL."/".$key."_R2.fastq.gz"); 
   }
   
}
close(READSET);

exit;
