#!/usr/bin/env python
 
"""Takes input in blastp against kegg database. (parsed output that is...)
and an abundance (gene count) table and output a KEGG-module abundance table.
Julien Tremblay - 
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

def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--infile-blastp', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-g', '--infile-gene-abundance', help="Input file", type=argparse.FileType('r'))
    
    args = parser.parse_args(arguments)
    infile_blastp = os.path.abspath(args.infile_blastp.name)
    infile_gene_abundance = os.path.abspath(args.infile_gene_abundance.name)
    sys.stderr.write(infile_blastp)
    sys.stderr.write(infile_gene_abundance)
     
    sys.stderr.write(str(args))

    # Parse gene abundance
    fhand=open(infile_blastp, "r")
    reader=csv.DictReader(fhand, delimiter="\t")

    hash_ko = defaultdict(list)
    hash_genes = defaultdict(list)
    
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
        pathway_ids = pathway_id_field.split("==")
        pathway_descs = pathway_desc_field.split("==")
        for pathway_id in pathway_ids:
            pathway_desc = pathway_descs.pop(0)
            #hash_ko[pathway_id + "--" + pathway_desc].append(gene_id)
            hash_ko[pathway_id].append(gene_id)
            hash_genes[gene_id] = ""

    # Parse gene abundance
    fhand = open(infile_gene_abundance, "r")
    gene_abundance_dict = defaultdict(list)
   

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
            if gene_id in hash_genes:
                value_list = row[1:len(row)]
                gene_abundance_dict[gene_id] = value_list
        
        j = j + 1
    
    print header
        
    # Then iterate and add values.
    for ko in hash_ko.iterkeys():
        # get gene id
        gene_ids = map(str, hash_ko[ko])
        # then add all rpkm values in gene_abundance_dict
        list_of_lists = []

        i = 0;
        curr_list = []
        for gene_id in gene_ids:
            #if gene_id in gene_abundance_dict.keys():
            #print gene_id
            curr_list = map(float, gene_abundance_dict[gene_id])
            # make sure that curr list is not empty. Faster to check empty list than to check for existence of key in huge hash.
            if len(curr_list) != 0:
                list_of_lists.append(curr_list)
            
            i = i + 1

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

            sys.stdout.write(ko + "\t")
            print "\t".join(map(str, curr_sum))
            #print "\t".join(str(curr_sum)[0])
        else:
            if not curr_list:
                x=1
            else:
                sys.stdout.write(ko + "\t")
                print "\t".join(map(str, curr_list))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


