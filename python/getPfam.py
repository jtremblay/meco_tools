#!/usr/bin/env python
 
"""Takes input in hmmsearch of genes vs pfam-A db. (parsed output that is...)
and an abundance (gene count) table and output a PFAM domain abundance table.
Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import re
import csv
import itertools
from collections import defaultdict
import operator
import numpy as np
import pprint

def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--infile-pfam', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-g', '--infile-gene-abundance', help="Input file", type=argparse.FileType('r'))
    
    args = parser.parse_args(arguments)
    infile_pfam = os.path.abspath(args.infile_pfam.name)
    infile_gene_abundance = os.path.abspath(args.infile_gene_abundance.name)
    sys.stderr.write(infile_pfam + "\n")
    sys.stderr.write(infile_gene_abundance + "\n")
     
    sys.stderr.write(str(args))

    # First parse blastp table.
    #print gene_abundance_dict
    fhand=open(infile_pfam, "r")
    #reader=csv.DictReader(fhand, delimiter="\t")
    reader=csv.reader(fhand, delimiter=" ")
    #gene_abundance_dict2 = gene_abundance_dict.copy()

    hash_pfam = defaultdict(list)
    hash_genes = defaultdict(list)
    p = re.compile("(PF\d+\.\d+)")

    for row in reader:
        if "#" not in row:
        
            gene_id  = row[2]
            pfam     = row[1]
            #pfam     = p.search(pfam)
            #pfam     = pfam.group(1)
            #sys.stderr.write("pfam:" + pfam + "\n")
            #sys.stderr.write("pfam:" + pfam + "\n")

            #hash_pfam[ko_id + "--" + definition].append(gene_id)
            hash_pfam[pfam].append(gene_id)
            hash_genes[gene_id] = ""
        #sys.stderr.write(ko_id + "\t" + definition + "\t" + gene_id + "\n")

    #pprint.pprint(hash_pfam)
    #exit(0)

    sys.stderr.write("Completed hashing parsed PFAM table\n")
    
    # Second, go through gene abundance list and store each value in a dict of lists.
    fhand = open(infile_gene_abundance, "r")
    gene_abundance_dict = defaultdict(list)

    j = 0;
    header = ""
    for line in fhand:
        if j == 0:
            line = line.rstrip('\n')
            row = line.split("\t")
            row[0] = "PFAM"
            header = "\t".join(row);
            #header = line
        else:    
            line = line.rstrip('\n')
            row = line.split("\t")
            gene_id = row[0]
            if gene_id in hash_genes:
                value_list = row[1:len(row)]
                gene_abundance_dict[gene_id] = value_list
        
        j = j + 1
   
    #sys.stderr.write("header: " + header + "\n")
    #pprint.pprint(gene_abundance_dict)
    #exit(0)
    
    sys.stderr.write("Completed hashing gene abundance table\n")
    
    print(header)

    #once blastp file has been parsed, loop through its hash_pfam.
    #for k in iter(d.keys())
    #for ko in hash_pfam.iterkeys():
    for pfam in iter(hash_pfam.keys()):
        # get gene id
        gene_ids = map(str, hash_pfam[pfam])
        #pprint.pprint(hash_pfam[pfam])
        # then add all values in gene_abundance_dict
        #sys.stderr.write("gene_ids: " + str(gene_ids) + "\n")
        list_of_lists = []

        i = 0;
        curr_list = []

        # First get all genes associated with that curr K orthology
        for gene_id in gene_ids:
            #if gene_id in gene_abundance_dict.keys():
            #print gene_id
            curr_list = list(map(float, gene_abundance_dict[gene_id]))
            #sys.stderr.write("curr_list inside : " + str(gene_abundance_dict[gene_id]) + "\n")
            # make sure that curr list is not empty. Faster to check empty list than to check for existence of key in huge hash.
            if len(list(curr_list)) != 0:
                list_of_lists.append(curr_list)
            
            del gene_abundance_dict[gene_id] # delete element found, because in the end, we'll print all elements not found.
            i = i + 1
             #else: # If does not exists. Still incude in gene list but with dummy KO (i.e. NA) values.
             #   curr_list = map(float, gene_abundance_dict[gene_id])
             #   curr_list = 
        #sys.stderr.write("curr_list: " + str(list_of_lists) + "\n")

        if i > 1:
            curr_sum = [sum(x) for x in zip(*list_of_lists)]

            #sys.stderr.write("curr_sum: " + str(curr_sum) + "\n")
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

            sys.stdout.write(pfam + "\t")
            print("\t".join(map(str, curr_sum)))
            #print "\t".join(str(curr_sum)[0])
        else:
            if not curr_list:
                x=1
            else:
                sys.stdout.write(pfam + "\t")
                print("\t".join(map(str, curr_list)))

    # Finally print gene_ids not found (not ided by blastp).
    #for gene_id in gene_abundance_dict2:
    #    curr_value = map(float, gene_abundance_dict[gene_id])
    #    ## Still not sure if we should print values without annotation... For now we won't...
    #    #print gene_id + "_NULL_KEGG\t" + "\t".join(map(str, curr_value))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


