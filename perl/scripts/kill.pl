#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
kill.pl

PURPOSE:

INPUT:
--all      : Kill all jobs
--blocked  : Kill blocked jobs
--active   : Kill active jobs
--eligible : Kill eligible jobs
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $all, $active,$blocked, $eligible);
my $verbose = 0;

GetOptions(
    'all' 	 => \$all,
	'blocked'  => \$blocked,
	'eligible' => \$eligible,
	'active'   => \$active,
    'verbose' 	 => \$verbose,
    'help' 		 => \$help
);
if ($help) { print $usage; exit; }

my $jobs = `showq -n -v -u jtrembla`;
my @jobs = split(/\n/, $jobs);

my @active;
my @eligible;
my @blocked;
my $flag;
foreach (@jobs){
	#print STDOUT $_."\n";
	if($_ =~ m/active jobs/){
		$flag = "A";
		next;
	}
	if($_ =~ m/eligible jobs/){
		$flag = "E";
		next;
	}
	if($_ =~ m/blocked jobs/){
		$flag = "B";
		next;
	}

	my $value;
	if($_ =~ m/^(\d+)\//){
		$value=$1;
		
		push(@active, $1) if($flag eq "A");
		push(@eligible, $1) if($flag eq "E");
		push(@blocked, $1) if($flag eq "B");
	}
}

if($blocked){ print STDOUT $_."\n" foreach(@blocked) };
if($active){ print STDOUT $_."\n" foreach(@active) };
if($eligible){ print STDOUT $_."\n" foreach(@eligible) };

if($all){
	print STDOUT $_."\n" foreach(@blocked);
	print STDOUT $_."\n" foreach(@active);
	print STDOUT $_."\n" foreach(@eligible);
}

exit;
