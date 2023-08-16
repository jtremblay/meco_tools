#!/usr/bin/env python
 
"""A simple python script template.

Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import re
import pickle
import h5py
import numpy as np

  
def vprint(obj):
    sys.stderr.write("[DEBUG] " + str(obj) + "\n")

def main(arguments):
    h5 = h5py.File('file.hdf5', 'r')
    ann = h5.get('annotations')
    abun = h5.get('gene_abundance')
    
    with open('index_annotations.pkl', 'rb') as f:
        index_annotations = pickle.load(f)
    
    with open('index_gene_abundance.pkl', 'rb') as f:
        index_gene_abundance = pickle.load(f)
    
    with open('index_kegg.pkl', 'rb') as f:
        index_kegg = pickle.load(f)
    
    curr_gene_id = "gene_id_3"
    index_ann = index_annotations[curr_gene_id]
    vprint("index_ann:" + str(index_ann))
    curr_ann = ann[index_ann]
    vprint(curr_ann)
    
    index_abun = index_gene_abundance[curr_gene_id]
    curr_abun = abun[index_abun]
    vprint(curr_abun)

    curr_kegg = "K00003"
    curr_kegg_genes = index_kegg[curr_kegg]
    vprint(curr_kegg_genes)
    curr_kegg_genes_indexes_abun = list( map(index_gene_abundance.get, curr_kegg_genes) )
    curr_kegg_genes_indexes_abun.sort()
    vprint(curr_kegg_genes_indexes_abun)
    #vprint(abun[curr_kegg_indexes_abun])
    curr_abun = abun[curr_kegg_genes_indexes_abun]
    print(curr_abun)
    #print(curr_kegg_indexes_abun)
    #curr_ = index_annotations[curr_
    


     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


