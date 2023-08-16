#!/usr/bin/env python
 
"""A simple python script template.

Julien Tremblay - 
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

    #gene_1 = dict(gene_id="gene_id_1", KO="K12345", symbol="aaa")
    #gene_2 = dict(gene_id="gene_id_2", KO="K12346", symbol="bbb")
    #gene_3 = dict(gene_id="gene_id_3", KO="K12347", symbol="ccc")

    #list_of_genes = [gene_1, gene_2, gene_3]
    #print(list_of_genes)
    #print(len(list_of_genes))
    #print(dict(list_of_genes))

    """Annotations
    """

    with h5py.File('file.hdf5','w') as h5w:
        dt = h5py.special_dtype(vlen=str)
        feature_names = np.array(
            [
                ['gene_id_1', 'K00001', 'COG0111'],
                ['gene_id_2', 'K00002', 'COG0222'],
                ['gene_id_4', 'K00001', 'COG0111'],
                ['gene_id_3', 'K00003', 'COG4441']
            ], 
            dtype=dt
        )
        h5w.create_dataset('annotations', data=feature_names, maxshape=(None, 3), compression="gzip", compression_opts=9)

        # Then add new data
        vprint(h5w['annotations'].shape[0])
        feature_names2 = np.array(
            [
                ['gene_id_5', 'K00002', 'COG0222'],
                ['gene_id_6', 'K00004', 'COG0555'],
                ['gene_id_7', 'K00003', 'COG4441']
            ], 
            dtype=dt
        )
        

        h5w['annotations'].resize((h5w['annotations'].shape[0] + 3), axis=0)
        h5w['annotations'][-feature_names2.shape[0]:] = feature_names2

        
        """Gene abundance
        """
        gene_abundance1 = np.array(
            [
                [33, 54, 66],
                [2, 3, 0],
                [11, 12, 15],
                [1, 0, 1],
            ],
            dtype='int32'
        )
        h5w.create_dataset('gene_abundance', data=gene_abundance1, maxshape=(None, 3), compression="gzip", compression_opts=9)
        gene_abundance2 = np.array(
            [
                [112, 145, 150],
                [50, 52, 52],
                [200, 201, 199]
            ],
            dtype='int32'
        )
        h5w['gene_abundance'].resize((h5w['gene_abundance'].shape[0] + 3), axis=0)
        h5w['gene_abundance'][-gene_abundance2.shape[0]:] = gene_abundance2

        """indexes
           gene_id_n = row[x] in the h5 'annotations' dataset
        """ 
        index_annotations = dict(
            gene_id_1 = 0, 
            gene_id_2 = 1, 
            gene_id_4 = 2, 
            gene_id_3 = 3, 
            gene_id_5 = 4, 
            gene_id_6 = 5, 
            gene_id_7 = 6
        )
        with open('index_annotations.pkl', 'wb') as f:
            pickle.dump(index_annotations, f)

        index_gene_abundance = dict(
            gene_id_1 = 6, 
            gene_id_2 = 5, 
            gene_id_4 = 4, 
            gene_id_3 = 3, 
            gene_id_5 = 2, 
            gene_id_6 = 0, 
            gene_id_7 = 1
        )
        with open('index_gene_abundance.pkl', 'wb') as f:
            pickle.dump(index_gene_abundance, f)

        index_kegg = dict(
            K00001 = ['gene_id_1', 'gene_id_4'],
            K00002 = ['gene_id_2', 'gene_id_5'],
            K00003 = ['gene_id_3', 'gene_id_7'],
            K00004 = ['gene_id_6']
        )
        with open('index_kegg.pkl', 'wb') as f:
            pickle.dump(index_kegg, f)
        

     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


