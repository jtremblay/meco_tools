#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
multilineBiomTable.pl

PURPOSE:

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
 - Biomonitoring
Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
   'infile=s' 	=> \$infile,
   'verbose' 	=> \$verbose,
   'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#open(IN, "<".$infile) or die "Can't open $infile\n";
#$/ = \1024;
#$/ = "}{";
#while(<IN>){
#   #print "\n" if $closing and /^{/;
#   #undef $closing;
#   #$_ =~ s/\{/\{\n/g;
#   #$_ =~ s/\}/\}\n/g;
#   $_ =~ s/\{/\{\n/g;
#   print STDOUT $_;
#   #$closing = 1 if /}$/;
#}

my $lastchunk = '';
my $buffer = '';
open(FILE, $infile) or die "Can't open `$infile': $!";
while (sysread FILE, $buffer, 1024) {
    my $matching = $lastchunk . $buffer;
    #if ($matching =~ /\{/) {
    #    print STDOUT "{\n";
    #}
    #if ($matching =~ /\}/) {
    #    print STDOUT "}\n";
    #}
    $matching =~ s/\{/\{\n/g;
    $matching =~ s/\}/\}\n/g;
    $matching =~ s/\[\[/\[\n\[/g;
    $matching =~ s/\]\]/\]\n\]/g;
    print STDOUT $matching;
    $lastchunk = $buffer;
}
close FILE;

#perl -i~ -e ' $/ = \1024;
#              while (<>) {
#                  print "\n" if $closing and /^{/;
#                  undef $closing;
#                  s/{/{\n}/g;
#                  s/}/}\n}/g;
#                  print;
#                  $closing = 1 if /}$/;
#              } ' 
#
