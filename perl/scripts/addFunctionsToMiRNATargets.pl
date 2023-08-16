#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
addFunctionsToMiRNATargets.pl

PURPOSE:

INPUT:
--infile <string>   : miRNA coordinate file (output of parseMiRNATargets.pl) 
--KO <string>       : output of kofamscan for each gene of each reference genome.
--COG <string>      : output of RPSBLAST COG for each gene of each reference genome.
--taxonomy <string> : Taxonomy file. 
--gff <string>      : merged gff file of all reference genome.


OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:

Julien Tremblay - jtremblay514@gmail.com

ENDHERE

## OPTIONS
my ($help, $infile, $gff, $KO, $COG, $taxonomy);
my $verbose = 0;

GetOptions(
   'infile=s'   => \$infile,
   'gff=s'      => \$gff,
   'KO=s'       => \$KO,
   'COG=s'      => \$COG,
   'taxonomy=s' => \$taxonomy,
   'verbose'    => \$verbose,
   'help'       => \$help
);
if ($help) { print $usage; exit; }

## SUBROUTINES
# Returns index of x in arr if present, else -1
sub binary_search{
    my ($arr, $low, $high, $x) = @_;
    my @arr = @{$arr};

    # because array is sorted.
    print STDERR "low: ".$low."   high:".$high."\n";

    if($high >= $low){
 
        my $mid = ($high + $low) / 2;
 
        # If element is present at the middle itself
        if($arr[$mid] == $x){
            return $mid;
 
        # If element is smaller than mid, then it can only
        # be present in left subarray
        }elsif($arr[$mid] > $x){
            return binary_search(\@arr, $low, $mid - 1, $x);
 
        # Else the element can only be present in right subarray
        }else{
            return binary_search(\@arr, $mid + 1, $high, $x);
        }
 
    }else{
        # Element is not present in the array
        return $high;
        #return -1;
    }
}

sub getNearest{
    my ($needle, $haystack) = @_;
    my @haystack = @{$haystack};
    if(!$haystack) {
        die('empty haystack');
    }
    my $bestDistance = 999999999999999999999999999;
    my $keep;
    my $keep_index = 0;
    my $index = 0;
    foreach my $value (@haystack) {
        if ($value == $needle) {
            #return $needle;
            return $index;
        }
        my $distance = abs($value - $needle);
        if ($distance < $bestDistance) {
            $bestDistance = $distance;
            $keep = $value;
            $keep_index = $index;
        }
        $index++;
    }
    return $keep_index;# ?? $value;
}

## MAIN

##my @array = (1,3,5,8,11,12,14,15,16,17,18,22,39,40);
##my $high = $array[$#array];
##my $low = $array[0];
###my $index = binary_search(\@array, $low, $high, 16);
##my $index = getNearest(16, \@array);
##print STDERR "found index (binary search): ".sprintf("%.0f", $index)."\n";
##exit;

my %hash;
# parse gff file.
open(IN, "<".$gff) or die "Can't open $gff\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    while(my $line = <IN>){
        next if($line =~ m/^#/);
        chomp;
        my @row = split(/\t/, $line);
        next if($row[2] ne "CDS");
        
        # find coords.
        my $contig = $row[0];
        my $start = $row[3];
        my $end = $row[4];
        my $strand = $row[6];

        # find gene id
        my $row8 = $row[8];
        my @row8 = split(/;/, $row8);
        my $gene = $row8[0];
        $gene =~ s/ID\=//;
        if($gene =~ m/\d+_(\d+)/){
            $gene = $1;
        }else{
            die("Could not parse gene id...\n");
        }
        $gene = $contig."_".$gene;

        $hash{$contig}{$gene}{start} = $start;
        $hash{$contig}{$gene}{end} = $end;
        $hash{$contig}{$gene}{strand} = $strand;
    }
}
close(IN);
#print STDERR Dumper(\%hash);
#exit;
# Parse taxonomy file
my %hash_taxonomy;
open(IN, "<".$taxonomy) or die "Can't open $taxonomy\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    while(my $line = <IN>){
        next if($line =~ m/^#/);
        chomp;
        my @row = split(/\t/, $line);
        
        my $genome_id = $row[0];
        my $s = $row[1];
        my $p = $row[5];
        my $c = $row[6];
        my $o = $row[7];
        my $f = $row[8];
        my $g = $row[9];

        $hash_taxonomy{$genome_id}{s} = $s;
        $hash_taxonomy{$genome_id}{p} = $p;
        $hash_taxonomy{$genome_id}{c} = $c;
        $hash_taxonomy{$genome_id}{o} = $o;
        $hash_taxonomy{$genome_id}{f} = $f;
        $hash_taxonomy{$genome_id}{g} = $g;
    }
}
close(IN);
#print STDERR Dumper(\%hash_taxonomy);

my %hash_cog;
open(IN, "<".$COG) or die "Can't open $COG\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    while(my $line = <IN>){
        next if($line =~ m/^#/);
        chomp;
        my @row = split(/\t/, $line);
        my $gene_id = $row[0];
        #if($gene_id =~ m/contig_id_(\d+_\d+)$/){
        #    $gene_id = $1;
        #}else{
        #    die("could not parse gene id from contig id...\n");
        #}

        my @terms = split(/,/, $row[12]);
        
        $hash_cog{$gene_id}{cog_access} = $terms[0];;
        my $term1 = $terms[1];
        $term1 =~ s/^\s+//;
        $hash_cog{$gene_id}{cog_name} = $term1;;
        my @desc = @terms[2 .. (@terms - 1)];
        my $desc = join(",", @desc);
        $desc =~ s/^\s+//;
        my $cog1;
        my $cog2;
        if($desc =~ m/^(.*)\[/){
           $cog1 = $1;
           $cog1 =~ s/\s+$//; 
           $hash_cog{$gene_id}{cog_function} = $cog1;
        }
        if($desc =~ m/\[(.*)\]/){
           $cog2 = $1; 
           $cog2 =~ s/\s+$//; 
           $hash_cog{$gene_id}{cog_category} = $cog2;
        }
    }
}
close(IN);

my %hash_kegg;
open(IN, "<".$KO) or die "Can't open $KO\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);

    while(my $line = <IN>){
        next if($line =~ m/^#/);
        chomp;
        my @row = split(/\t/, $line);
        my $gene_id = $row[0];
        #if($gene_id =~ m/contig_id_(\d+_\d+)$/){
        #    $gene_id = $1;
        #}else{
        #    die("could not parse gene id from contig id...\n");
        #}
        
        my $KO = $row[1];
        my $gene_name = $row[2];
        my $KO_desc = $row[4];
        my $pathway_id = $row[5];
        my $pathway_name = $row[6];
        my $module_id = $row[7];
        my $module_name = $row[8];

        $hash_kegg{$gene_id}{KO} = $KO;
        $hash_kegg{$gene_id}{gene_name} = $gene_name;
        $hash_kegg{$gene_id}{KO_desc} = $KO_desc;
        $hash_kegg{$gene_id}{pathway_id} = $pathway_id;
        $hash_kegg{$gene_id}{pathway_name} = $pathway_name;
        $hash_kegg{$gene_id}{module_id} = $module_id;
        $hash_kegg{$gene_id}{module_name} = $module_name;
    }
}
close(IN);

# lastly, find on which gene corresponds each miRNA coordinate hits.
# Also include the position from each 3', 5' target from its orf.
# Put the miRNA table data into a hash.
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    if($. == 1){ # print header
        print STDOUT $_."\ttarget_gene\tcategory\tdistance_from_5p\tdistance_from_3p\tdistance_from_left\tdistance_from_right\t";
        print STDOUT "phylum\tclass\torder\tfamily\tgenus\tspecies\tcog_accession\tcog_name\tcog_function\tcog_category\t";
        print STDOUT "KO\tgene_name\tKO_desc\tpathway_id\tpathway_name\tmodule_id\tmodule_name\n";
        next;
    }
    my @row = split(/\t/, $_);
   
    my $mirna = $row[0];
    my $contig = $row[1];
    my $mirna_start = $row[7];
    my $mirna_end = $row[8];

    if(exists $hash{$contig}){
        # First, get taxonomy.
        my ($genome_id) = $contig =~ m/^(\d+)_.*/;
        my $s = $hash_taxonomy{$genome_id}{s};
        my $p = $hash_taxonomy{$genome_id}{p};
        my $c = $hash_taxonomy{$genome_id}{c};
        my $o = $hash_taxonomy{$genome_id}{o};
        my $f = $hash_taxonomy{$genome_id}{f};
        my $g = $hash_taxonomy{$genome_id}{g};

        # we have a candidate contig. Lets put the coords of genes into an array.
        my @array_start;
        my @array_end;
        my @array_gene_id;

        my @array;
        foreach my $gene_id (keys %{ $hash{$contig} }) {
            my $start = $hash{$contig}{$gene_id}{start};
            my $end = $hash{$contig}{$gene_id}{end};
            push(@array, {gene_id => $gene_id, start => $start, end => $end});
        }

        # Then sort the array by the gff start coordinates
        foreach my $gene_id ( sort { $a->{start} <=> $b->{start} } @array ) {
            push(@array_gene_id, $gene_id->{'gene_id'});
            push(@array_start,   $gene_id->{'start'});
            push(@array_end,     $gene_id->{'end'});
        }

        # return index
        my $index_mirna_start_vs_gff_start = getNearest($mirna_start, \@array_start);
        my $index_mirna_start_vs_gff_end   = getNearest($mirna_start, \@array_end);
        my $index_mirna_end_vs_gff_start   = getNearest($mirna_end, \@array_start);
        my $index_mirna_end_vs_gff_end     = getNearest($mirna_end, \@array_end);

        # get value from index.
        my $value_mirna_start_vs_gff_start = $array_start[$index_mirna_start_vs_gff_start]; 
        my $value_mirna_start_vs_gff_end   = $array_end[$index_mirna_start_vs_gff_end];
        my $value_mirna_end_vs_gff_start   = $array_start[$index_mirna_end_vs_gff_start];
        if(!defined $value_mirna_end_vs_gff_start){
            print STDERR $_."\n";
            print STDERR $index_mirna_end_vs_gff_start."\n";
            exit;
        }
        my $value_mirna_end_vs_gff_end     = $array_end[$index_mirna_end_vs_gff_end];
        
        # Once we have the indexes we can then fetch the values and can 
        # determine the shortest distances.
        # Once we find the shortest distance, we can then find the gene id using that index position.
        my $shortest_dist_mirna_start = 999999999999999999;
        my $chosen_gene_mirna_start;
        if(abs($mirna_start - $value_mirna_start_vs_gff_start) < $shortest_dist_mirna_start){
            $shortest_dist_mirna_start = abs($mirna_start - $value_mirna_start_vs_gff_start);
            $chosen_gene_mirna_start = $array_gene_id[$index_mirna_start_vs_gff_start];
        }
        if(abs($mirna_start - $value_mirna_start_vs_gff_end) < $shortest_dist_mirna_start){
            $shortest_dist_mirna_start = abs($mirna_start - $value_mirna_start_vs_gff_end);
            $chosen_gene_mirna_start = $array_gene_id[$index_mirna_start_vs_gff_end];
        }
        
        my $shortest_dist_mirna_end = 999999999999999999;
        my $chosen_gene_mirna_end;
        if(abs($mirna_end - $value_mirna_end_vs_gff_start) < $shortest_dist_mirna_end){
            $shortest_dist_mirna_end = abs($mirna_end - $value_mirna_end_vs_gff_start);
            $chosen_gene_mirna_end = $array_gene_id[$index_mirna_end_vs_gff_start];
        }
        if(abs($mirna_end - $value_mirna_end_vs_gff_end) < $shortest_dist_mirna_end){
            $shortest_dist_mirna_end = abs($mirna_end - $value_mirna_end_vs_gff_end);
            $chosen_gene_mirna_end = $array_gene_id[$index_mirna_end_vs_gff_end];
        }

        # At this point we have our chosen gene. Which category is the interaction mirna:target gene?
        #Determine where exactly is located the miRNA target on the genome:
        #Is it 
        #5', 
        #5' overlapping with CDS, 
        #inside CDS, 
        #3' overlapping with CDS 
        #3'
        #or between 2 genes.
        my $chosen_gene;
        my $category = "undef";
        my $distance_from_5p;
        my $distance_from_3p;
        my $distance_from_left;
        my $distance_from_right;
        if($chosen_gene_mirna_end eq $chosen_gene_mirna_start){
            $chosen_gene = $chosen_gene_mirna_start;
            
            #functions
            my $cog_access; my $cog_name; my $cog_function; my $cog_category;
            if(exists $hash_cog{$chosen_gene}){
                $cog_access = $hash_cog{$chosen_gene}{cog_access};
                $cog_name = $hash_cog{$chosen_gene}{cog_name};
                $cog_function = $hash_cog{$chosen_gene}{cog_function};
                $cog_category = $hash_cog{$chosen_gene}{cog_category};
            }else{
                $cog_access = "undef";
                $cog_name = "undef";
                $cog_function = "undef";
                $cog_category = "undef";
            }

            my $KO; my $KO_desc; my $gene_name; my $pathway_id; my $pathway_name; my $module_id; my $module_name;
            if(exists $hash_kegg{$chosen_gene}){
                $KO = $hash_kegg{$chosen_gene}{KO};
                $KO_desc = $hash_kegg{$chosen_gene}{KO_desc};
                $gene_name = $hash_kegg{$chosen_gene}{gene_name};
                $pathway_id = $hash_kegg{$chosen_gene}{pathway_id};
                $pathway_name = $hash_kegg{$chosen_gene}{pathway_name};
                $module_id = $hash_kegg{$chosen_gene}{module_id};
                $module_name = $hash_kegg{$chosen_gene}{module_name};
            }else{
                $KO = "undef";
                $gene_name = "undef";
                $pathway_id = "undef";
                $pathway_name = "undef";
                $KO_desc = "undef";
                $module_id = "undef";
                $module_name = "undef";
            }

            #print STDERR Dumper($hash{$contig}{$chosen_gene});

            # 5' or 3' or 5' overlapping or 3'overlapping?
            if($hash{$contig}{$chosen_gene}{strand} eq "+"){
                # inside CDS
                #                              start             (+)               end
                # -------------------------------5'=================================>3'---------------------------
                #                                        start          end
                if($mirna_start >= $hash{$contig}{$chosen_gene}{start} && $mirna_end <= $hash{$contig}{$chosen_gene}{end}){
                    $category = "inside_CDS";
                    $distance_from_5p = $mirna_start - $hash{$contig}{$chosen_gene}{start};
                    $distance_from_3p = $hash{$contig}{$chosen_gene}{end} - $mirna_end;

                # is it in  5' nc region?
                #                              start             (+)                end
                # -------------------------------5'=================================>3'---------------------------
                #           start          end
                }elsif($mirna_start < $hash{$contig}{$chosen_gene}{start} && $mirna_end <= $hash{$contig}{$chosen_gene}{start}){
                    $category = "5p";
                    $distance_from_5p = $hash{$contig}{$chosen_gene}{start} - $mirna_end; # always the shortest distance between gene coord and mirna target position. In this case the mirna_end is the closest to the gene.
                    $distance_from_3p = $hash{$contig}{$chosen_gene}{end} - $mirna_end;

                # is it in  3' nc region?
                #                              start             (+)                end
                # -------------------------------5'=================================>3'---------------------------
                #                                                                         start          end
                }elsif($mirna_start >= $hash{$contig}{$chosen_gene}{end} && $mirna_end > $hash{$contig}{$chosen_gene}{end}){
                    $category = "3p";
                    $distance_from_5p = $mirna_start - $hash{$contig}{$chosen_gene}{start};
                    $distance_from_3p = $mirna_start - $hash{$contig}{$chosen_gene}{end};
                                
                #if not, is it overlapping in 5' and nc region?
                #                              start             (+)                end
                # -------------------------------5'=================================>3'---------------------------
                #                        start          end
                }elsif($mirna_start < $hash{$contig}{$chosen_gene}{start} && $mirna_end > $hash{$contig}{$chosen_gene}{start}){
                    $category = "5p_overlapping";
                    $distance_from_5p = "NA";
                    $distance_from_3p = "NA";
                
                #if not, is it overlapping in 3' and nc region?
                #                              start             (+)                end
                # -------------------------------5'=================================>3'---------------------------
                #                                                           start          end
                }elsif($mirna_start < $hash{$contig}{$chosen_gene}{end} && $mirna_end > $hash{$contig}{$chosen_gene}{end}){
                    $category = "3p_overlapping";
                    $distance_from_5p = "NA";
                    $distance_from_3p = "NA";
                }else{
                    die("Uncaught exception:\n$chosen_gene\n$_\n");
                }
                

            }elsif($hash{$contig}{$chosen_gene}{strand} eq "-"){
                # inside CDS
                #                                start            (-)              end
                # -------------------------------3'<================================5'---------------------------
                #                                        start          end
                if($mirna_start >= $hash{$contig}{$chosen_gene}{start} && $mirna_end <= $hash{$contig}{$chosen_gene}{end}){
                    $category = "inside_CDS";
                    $distance_from_5p = $hash{$contig}{$chosen_gene}{end} - $mirna_end;
                    $distance_from_3p = $mirna_start - $hash{$contig}{$chosen_gene}{start};
            
                # 5' nc region
                #                                start            (-)              end
                # -------------------------------3'<================================5'---------------------------
                #                                                                           start          end
                }elsif($mirna_start >= $hash{$contig}{$chosen_gene}{end} && $mirna_end > $hash{$contig}{$chosen_gene}{end}){
                    $category = "5p";
                    $distance_from_5p = $mirna_start - $hash{$contig}{$chosen_gene}{end};
                    $distance_from_3p = $mirna_start - $hash{$contig}{$chosen_gene}{start};
                
                # 3' nc region
                #                                start            (-)              end
                # -------------------------------3'<================================5'---------------------------
                #         start          end
                }elsif($mirna_start < $hash{$contig}{$chosen_gene}{start} && $mirna_end <= $hash{$contig}{$chosen_gene}{start}){
                    $category = "3p";
                    $distance_from_5p = $hash{$contig}{$chosen_gene}{end} - $mirna_end;
                    $distance_from_3p = $hash{$contig}{$chosen_gene}{start} - $mirna_end;
                                
                #if not, is it overlapping in 5' and nc?
                #                                start            (-)              end
                # -------------------------------3'<================================5'---------------------------
                #                                                          start          end
                }elsif($mirna_start < $hash{$contig}{$chosen_gene}{end} && $mirna_end > $hash{$contig}{$chosen_gene}{end}){
                    $category = "5p_overlapping";
                    $distance_from_5p = "NA";
                    $distance_from_3p = "NA";
                
                #if not, is it overlapping in 3' and nc?
                #                                start            (-)              end
                # -------------------------------3'<================================5'---------------------------
                #                        start          end
                }elsif($mirna_start < $hash{$contig}{$chosen_gene}{start} && $mirna_end > $hash{$contig}{$chosen_gene}{start}){
                    $category = "3p_overlapping";
                    $distance_from_5p = "NA";
                    $distance_from_3p = "NA";
                }else{
                    print STDERR Dumper($hash{$contig}{$chosen_gene});
                    print STDERR "mirna_start: ".$mirna_start."\n";
                    print STDERR "mirna_end: ".$mirna_end."\n";
                    die("Exception:\n$chosen_gene\n$_\n");
                }
            }
            $distance_from_left = "NA";
            $distance_from_right = "NA";

            print STDOUT $_."\t".$chosen_gene."\t".$category."\t".$distance_from_5p."\t".$distance_from_3p."\t".$distance_from_left."\t".$distance_from_right."\t";
            print STDOUT $p."\t".$c."\t".$o."\t".$f."\t".$g."\t".$s."\t".$cog_access."\t".$cog_name."\t".$cog_function."\t".$cog_category."\t";
            print STDOUT $KO."\t".$KO_desc."\t".$gene_name."\t".$pathway_id."\t".$pathway_name."\t".$module_id."\t".$module_name."\n";
        }else{
            #If chosen_gene_mirna_end and chosen_gene_mirna_start disagree, it means 
            #that we probably have a microrna localized between two genes at relatively the same distance of each gene...

            #functions
            my $cog_access_left; my $cog_name_left; my $cog_function_left; my $cog_category_left;
            if(exists $hash_cog{$chosen_gene_mirna_start}){
                $cog_access_left = $hash_cog{$chosen_gene_mirna_start}{cog_access};
                $cog_name_left = $hash_cog{$chosen_gene_mirna_start}{cog_name};
                $cog_function_left = $hash_cog{$chosen_gene_mirna_start}{cog_function};
                $cog_category_left = $hash_cog{$chosen_gene_mirna_start}{cog_category};
            }else{
                $cog_access_left = "undef";
                $cog_name_left = "undef";
                $cog_function_left = "undef";
                $cog_category_left = "undef";
            }
            
            my $cog_access_right; my $cog_name_right; my $cog_function_right; my $cog_category_right;
            if(exists $hash_cog{$chosen_gene_mirna_end}){
                $cog_access_right = $hash_cog{$chosen_gene_mirna_end}{cog_access};
                $cog_name_right = $hash_cog{$chosen_gene_mirna_end}{cog_name};
                $cog_function_right = $hash_cog{$chosen_gene_mirna_end}{cog_function};
                $cog_category_right = $hash_cog{$chosen_gene_mirna_end}{cog_category};
            }else{
                $cog_access_right = "undef";
                $cog_name_right = "undef";
                $cog_function_right = "undef";
                $cog_category_right = "undef";
            }
            
            my $KO_left; my $KO_desc_left; my $gene_name_left; my $pathway_id_left; my $pathway_name_left; my $module_id_left; my $module_name_left;
            if(exists $hash_kegg{$chosen_gene_mirna_start}){
                $KO_left = $hash_kegg{$chosen_gene_mirna_start}{KO};
                $KO_desc_left = $hash_kegg{$chosen_gene_mirna_start}{KO_desc};
                $gene_name_left = $hash_kegg{$chosen_gene_mirna_start}{gene_name};
                $pathway_id_left = $hash_kegg{$chosen_gene_mirna_start}{pathway_id};
                $pathway_name_left = $hash_kegg{$chosen_gene_mirna_start}{pathway_name};
                $module_id_left = $hash_kegg{$chosen_gene_mirna_start}{module_id};
                $module_name_left = $hash_kegg{$chosen_gene_mirna_start}{module_name};
            }else{
                $KO_left = "undef";
                $KO_desc_left = "undef";
                $gene_name_left = "undef";
                $pathway_id_left = "undef";
                $pathway_name_left = "undef";
                $module_id_left = "undef";
                $module_name_left = "undef";
            }
            
            my $KO_right; my $KO_desc_right; my $gene_name_right; my $pathway_id_right; my $pathway_name_right; my $module_id_right; my $module_name_right;
            if(exists $hash_kegg{$chosen_gene_mirna_end}){
                $KO_right = $hash_kegg{$chosen_gene_mirna_end}{KO};
                $KO_desc_right = $hash_kegg{$chosen_gene_mirna_end}{KO_desc};
                $gene_name_right = $hash_kegg{$chosen_gene_mirna_end}{gene_name};
                $pathway_id_right = $hash_kegg{$chosen_gene_mirna_end}{pathway_id};
                $pathway_name_right = $hash_kegg{$chosen_gene_mirna_end}{pathway_name};
                $module_id_right = $hash_kegg{$chosen_gene_mirna_end}{module_id};
                $module_name_right = $hash_kegg{$chosen_gene_mirna_end}{module_name};
            }else{
                $KO_right = "undef";
                $KO_desc_right = "undef";
                $gene_name_right = "undef";
                $pathway_id_right = "undef";
                $pathway_name_right = "undef";
                $module_id_right = "undef";
                $module_name_right = "undef";
            }
            
            my $left_strand = $hash{$contig}{$chosen_gene_mirna_start}{strand};
            my $right_strand = $hash{$contig}{$chosen_gene_mirna_end}{strand};
            $category = "between_2_genes_".$left_strand."_".$right_strand;
            $distance_from_left = $hash{$contig}{$chosen_gene_mirna_start}{start} - $mirna_end;
            $distance_from_right = $hash{$contig}{$chosen_gene_mirna_end}{end} - $mirna_end ;
            $distance_from_5p = "NA";
            $distance_from_3p = "NA";
        
            print STDOUT $_."\t".$chosen_gene_mirna_start.";".$chosen_gene_mirna_end."\t".$category."\t".$distance_from_5p."\t".$distance_from_3p."\t".$distance_from_left."\t".$distance_from_right."\t";
            print STDOUT $p."\t".$c."\t".$o."\t".$f."\t".$g."\t".$s."\t".$cog_access_left."|".$cog_access_right."\t".$cog_name_left."|".$cog_name_right."\t".$cog_function_left."|".$cog_function_right."\t".$cog_category_left."|".$cog_category_right."\t";
            print STDOUT $KO_left."|".$KO_right."\t".$KO_desc_left."|".$KO_desc_right."\t".$gene_name_left."|".$gene_name_right."\t".$pathway_id_left."|".$pathway_id_right."\t".$pathway_name_left."|".$pathway_name_right."\t".$module_id_left."|".$module_id_right."\t".$module_name_left."|".$module_name_right."\n";
        }
    }
}
close(IN);
