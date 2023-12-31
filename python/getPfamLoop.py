#!/usr/bin/env python
 
"""Takes input in rps-blast format against PFAM database. (parsed output that is...)
and an abundance (gene count) table and output a PFAM abundance table.
Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import re
import csv
import fnmatch
import multiprocessing
from collections import defaultdict
import operator
import numpy as np

def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--infile-pfam', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-t', '--num-threads', help="Number of threads", type=int)
    #parser.add_argument('-g', '--infile-gene-abundance', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-i', '--indir', help="Input directory where are your DDA.tsv files.")
    parser.add_argument('-p', '--prefix', help="Prefix of your expression files (DDA or DEG).")
   

    args = parser.parse_args(arguments)
    infile_pfam = os.path.abspath(args.infile_pfam.name)
    indir = os.path.join(args.indir)
    num_threads = args.num_threads
    prefix = str(args.prefix)

    sys.stderr.write("infile pfam: " + infile_pfam + "\n")

    # First find each DDA.tsv files.
    matches = []
    for root, dirnames, filenames in os.walk(indir):
        for filename in fnmatch.filter(filenames, '*' + prefix + '_normalized_significant.tsv'):
            matches.append(os.path.join(root, filename))

    sys.stderr.write("\n".join(matches))

    # Do this loop in multiprocessing.
    #q = Queue(maxsize=0)
    #get_all_files(q, matches, indir, infile_pfam)
    #num_threads = 11
    p = multiprocessing.Pool(num_threads)
    for match in matches:
        p.apply_async(do_work, [match, indir, infile_pfam])
    
    #p.apply_async(do_work, [matches[1], indir, infile_pfam])
    
    p.close()
    p.join()

    
def do_work(match, indir, infile_pfam):    
    #param = q.get()
    #match = param[0]
    #indir = param[1]
    #infile_pfam = param[2]
    
    #while True:
        # Do work....
    outfile = os.path.splitext(match)[0] + "_PFAM.tsv"
    OUT = open(outfile, 'w')
    infile_gene_abundance = match
    
    sys.stderr.write("\n===BEGIN===\n" + match + "\n")
    sys.stderr.write("\nwriting to outfile: " + outfile + "\n")
    
    fhand = open(infile_gene_abundance, "r")
    gene_abundance_dict = defaultdict(list)
    
    # First, go through gene abundance list and store each value in a dict of lists.
    j = 0;
    header = ""
    for line in fhand:
        if j == 0:
            line = line.rstrip('\n')
            header = line
        else:    
            line = line.rstrip('\n')
            row = line.split("\t")
            gene_id = row[0]
            value_list = row[1:len(row)]
            gene_abundance_dict[gene_id] = value_list
        
        j = j + 1
    fhand.close() 
    # Next parse blastp table.
    #print gene_abundance_dict
    
    hash_ko = defaultdict(list)
    
    OUT.write(header + "\n")
    
    my_pfam_file = open(infile_pfam, "rU")
    #reader=csv.DictReader(fhand, delimiter="\t")

    #k=0
    # skipi first 3 lines (header)
    #for _ in xrange(3):
    #    next(my_pfam_file)

    seen = dict()
    for line in my_pfam_file:
        line = line.rstrip('\n')
        #:sys.stderr.write(line + "\n")

        if(line[:1] == '#'):
            sys.stderr.write(line[:1] + "\n")
        else:
            #sys.stderr.write("FIRST LINE NOT EQUAL TO #: " + line[:1] + "\n")
            row = re.split('\s+', line)
            gene_id     = row[2]
            pfam_id     = row[1]
            pfam_id2    = row[0]
            descs       = row[18:len(row)]
            desc        = ' '.join(descs)
    
            if(gene_id in seen):
                pass
            else:
                #sys.stderr.write("desc: " + desc + "\n")
                #sys.stderr.write("gene_id: " + gene_id + "\n")
                #sys.stderr.write("pfam_id: " + pfam_id + "\n")
                #sys.stderr.write("pfam_id2: " + pfam_id2 + "\n")
                hash_ko[pfam_id2 + "==" + pfam_id + "==" + desc].append(gene_id)
                seen[gene_id] = gene_id
    
    sys.stderr.write("complete processing hmmscan file...\n")
    my_pfam_file.close()
    gene_abundance_dict2 = gene_abundance_dict.copy()
    #once blastp file has been parsed, loop through its hash_ko.


    for ko in hash_ko.iterkeys():
        # get gene id
        gene_ids = map(str, hash_ko[ko])
        # then add all rpkm values in gene_abundance_dict
        list_of_lists = []
    
        i = 0;
        curr_list = []
        for gene_id in gene_ids:
            if gene_id in gene_abundance_dict.keys():
                sys.stderr.write("gene_id" +  gene_id + "\n")
                curr_list = map(float, gene_abundance_dict[gene_id])
                #del gene_abundance_dict[gene_id] # delete element found, because in the end, we'll print all elements not found.
                list_of_lists.append(curr_list)
                
                i = i + 1
             #else: # If does not exists. Still incude in gene list but with dummy KO (i.e. NA) values.
             #   curr_list = map(float, gene_abundance_dict[gene_id])
             #   curr_list = 
    
        if i > 1:
            curr_sum = [sum(x) for x in zip(*list_of_lists)]
            #curr_sum = []
            #k = 0
            #for this_list in list_of_lists:
            #    if k == 0:
            #        curr_sum = this_list
            #        #print curr_sum
            #    else:
            #        curr_sum = np.sum([curr_sum, this_list], axis=0)
            #        #print curr_sum
            #    k = k + 1
    
            OUT.write(ko + "\t")
            OUT.write("\t".join(map(str, curr_sum)) + "\n")
            #print "\t".join(map(str, curr_sum))
            #print "\t".join(str(curr_sum)[0])
        else:
            if not curr_list:
                x=1
            else:
                OUT.write(ko + "\t")
                OUT.write("\t".join(map(str, curr_list)) + "\n")
                #print "\t".join(map(str, curr_list))
    
    # Finally print gene_ids not found (not ided by blastp).
    for gene_id in gene_abundance_dict2:
        curr_value = map(float, gene_abundance_dict[gene_id])
        ## Still not sure if we should print values without annotation... For now we won't...
        #print gene_id + "_NULL_KEGG\t" + "\t".join(map(str, curr_value))
    
    OUT.close
    #q.task_done()

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


