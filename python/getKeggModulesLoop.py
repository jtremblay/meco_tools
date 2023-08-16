#!/usr/bin/env python
 
"""Takes input in blastp against kegg database. (parsed output that is...)
and an abundance (gene count) table and output a Kxxxx abundance table.
Julien Tremblay - 
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
    parser.add_argument('-b', '--infile-blastp', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-t', '--num-threads', help="Number of threads", type=int)
    parser.add_argument('-i', '--indir', help="Input directory where are your DDA.tsv files.")
    parser.add_argument('-p', '--prefix', help="Prefix of your expression files (DDA or DEG).")
   

    args = parser.parse_args(arguments)
    infile_blastp = os.path.abspath(args.infile_blastp.name)
    indir = os.path.join(args.indir)
    num_threads = args.num_threads
    prefix = str(args.prefix)

    # First find each DDA.tsv files.
    matches = []
    for root, dirnames, filenames in os.walk(indir):
        for filename in fnmatch.filter(filenames, '*' + prefix + '_normalized_significant.tsv'):
            matches.append(os.path.join(root, filename))

    sys.stderr.write("\n".join(matches))

    # Do this loop in multiprocessing.
    #q = Queue(maxsize=0)
    #get_all_files(q, matches, indir, infile_blastp)
    #num_threads = 11
    p = multiprocessing.Pool(num_threads)
    for match in matches:
        p.apply_async(do_work, [match, indir, infile_blastp])
    
    p.close()
    p.join()

    
def do_work(match, indir, infile_blastp):    
    #param = q.get()
    #match = param[0]
    #indir = param[1]
    #infile_blastp = param[2]
    
    #while True:
        # Do work....
    outfile = os.path.splitext(match)[0] + "_keggModules.tsv"
    OUT = open(outfile, 'w')
    infile_gene_abundance = match
    
    sys.stderr.write("===BEGIN===\n" + match + "\n")
    
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
    
    # Next parse blastp table.
    #print gene_abundance_dict
    fhand=open(infile_blastp, "r")
    reader=csv.DictReader(fhand, delimiter="\t")
    
    hash_ko = defaultdict(list)
    
    OUT.write(header + "\n")
    
    for row in reader:
        gene_id            = row["#query"]
        kegg_gene_id       = row["kegg_gene_id"]
        ko_id              = row["KO_id"]
        entry              = row["ENTRY"]
        definition         = row["DEFINITION"]
        pathway_id_field   = row["PATHWAY_ID"]
        pathway_desc_field = row["PATHWAY_DESC"]
        module_id_field    = row["MODULE_ID"]
        module_desc_field  = row["MODULE_DESC"]
    
        # split MODULE_ID AND MODULE_DESC
        module_ids = module_id_field.split("==")
        module_descs = module_desc_field.split("==")
    
        for module_id in module_ids:
            module_desc = module_descs.pop(0)
            if module_id != "NULL":
                hash_ko[module_id + "--" + module_desc].append(gene_id)
                #print ko_id + "\t" + definition + "\t" + gene_id + "\n"
    
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
                #print gene_id
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


