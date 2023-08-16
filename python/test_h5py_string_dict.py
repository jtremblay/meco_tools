#!/usr/bin/env python
 
"""A simple python script template.

Julien Tremblay - 
"""
 
import os
import sys
import argparse
import re
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

    with h5py.File('file.hdf5','w') as h5w:
        dt = h5py.special_dtype(vlen=str)
        feature_names = dict(
            gene_id1 = ['aaa', 'bbb', 'cccccccccccccccccccccccccccc'],
            gene_id2 = ['bbb', '123', 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'],
            gene_id3 = ['cc', '12345', 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'],
            gene_id4 = ['dddb', '124566', 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee']
        ) 
        h5w.create_dataset('annotations', data=feature_names, compression="gzip", compression_opts=9)
        #h5w.create_dataset('annotations', data=feature_names, maxshape=(None, 3), compression="gzip", compression_opts=9)
        exit(0);

        # Then add new data
        vprint(h5w['annotations'].shape[0])
        feature_names2 = np.array(
            [
                ['JJJ', 'bcd', 'cccccccccccc8990-cc'],
                ['HH', '1255555', 'eeeeee877eeeeeeeee'],
                ['VVV', '124563333', '658eeeeeeeeeeee']
            ], 
            dtype=dt
        )



        h5w['annotations'].resize((h5w['annotations'].shape[0] + 3), axis=0)
        vprint(h5w['annotations'][-feature_names2.shape[0]:])
        h5w['annotations'][-feature_names2.shape[0]:] = feature_names2


        ## Then query file:
        vprint(h5w.keys())
        with h5py.File('file.hdf5', "r") as f:
        # List all groups
            print("Keys: %s" % f.keys())
            a_group_key = list(f.keys())[0]

             # Get the data
            data = list(f[a_group_key])
            vprint(data[0][0])


     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


